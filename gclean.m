function gclean(cfg, subj)
%GCLEAN use German's toolbox to clean the data
% Read the data from seldata (it assumes that the data is continuous), low
% pass filter, reject bad channels and eye, ecg, emg activity with ICA.
% It returns a fieldtrip structure made of trials of different length only
% containing good data. Bad channels are interpolated.
%
% CFG
%  .data: path of /data1/projects/PROJ/subjects/
%  .nick: NICK in /data1/projects/PROJ/subjects/0001/MOD/NICK/
%  .mod: modality, MOD in /data1/projects/PROJ/subjects/0001/MOD/NICK/
%  .endname: includes preprocessing steps (e.g. '_seldata')
%
%  .log: name of the file and directory to save log
%
%  .sens.file: file with EEG sensors. It can be sfp or mat
%  .sens.dist: distance between sensors to consider them neighbors (in the units of cfg.sens.file)
%
%  .step: all the analysis step (for cfg.clear)
%  .clear: cell with the name of preprocessing steps to delete
%
%  .gtool.fsample: manually specify the frequency (very easily bug-prone, but in this way it does not read "data" all the time)
%  .gtool.saveall: false
%  .gtool.verbose: true
%  .gtool.lpfreqn = [.5 / (cfg.gtool.fsample/2)]; % normalized by half of the sampling frequency!
%  .gtool.bad_samples.MADs = 5;
%  .gtool.bad_channels.MADs = 8;
%  .gtool.eog.correction = 50;
%  .gtool.emg.correction = 30;
%  .gtool.pwl.pca: if not empty, number of PC to work with in JADE
%
% IN
%  data in /data1/projects/PROJ/subjects/SUBJ/MOD/NICK/
% 
% OUT
%  data, after ICA rejection for eyes-blinks and movement, rejection and
%  interpolation of bad electrodes, rejection of bad epochs
%  It appends '_gclean' at the end of the filename
%
% Part of EVENTBASED preprocessing
% see also SELDATA, GCLEAN, REDEF

%---------------------------%
%-start log
output = sprintf('%s (%04d) began at %s on %s\n', ...
  mfilename, subj, datestr(now, 'HH:MM:SS'), datestr(now, 'dd-mmm-yy'));
tic_t = tic;
%---------------------------%

%---------------------------%
%-dir and files
ddir = sprintf('%s%04d/%s/%s/', cfg.data, subj, cfg.mod, cfg.nick); % data dir
allfile = dir([ddir '*' cfg.endname '.mat']); % files matching a preprocessing
%---------------------------%

%---------------------------%
%-elec or grad
haselec = false;
if isfield(cfg.sens, 'file') && ~isempty(cfg.sens.file)
  haselec = true;
  sens = ft_read_sens(cfg.sens.file);
  sens.label = upper(sens.label);
  
  cfg1 = [];
  cfg1.elec = sens;
  cfg1.method = 'distance';
  cfg1.neighbourdist = cfg.sens.dist;
  neigh = ft_prepare_neighbours(cfg1);
end

if strcmp(cfg.mod, 'meg')
  hasgrad = true;
else
  hasgrad = false;
end
%---------------------------%

%------------------------------------%
%-loop over files
import pset.node.*

for i = 1:numel(allfile)
  
  %-----------------%
  %-load the data
  outtmp = sprintf('%s\n', allfile(i).name);
  output = [output outtmp];
  %-----------------%
  
  %--------------------------%
  %-gtoolbox parameters
  
  %------%
  %-PCA
  if isfield(cfg.gtool, 'pwl') && ...
      isfield(cfg.gtool.pwl, 'pca')
    opt.pcapwl = spt.pca('MaxDimOut', cfg.gtool.pwl.pca);
  else
    opt.pcapwl = [];
  end
  %------%
  
  %------%
  %-JADE
  opt.bsspwl          = spt.jade;
  opt.bsseog          = spt.jade;
  opt.bssecg          = spt.jade;
  opt.bssemg          = spt.jade;
  opt.bssrs           = spt.jade;
  %------%
  
  %------%
  %-filter
  if ~isempty(cfg.gtool.lpfreqn)
    filterNode1 = ...
      pset.node.bpfilt('fp', [cfg.gtool.lpfreqn 1], ...
      'save', cfg.gtool.saveall);
  else
    filterNode1 = [];
  end
  %------%
  %--------------------------%
  
  %--------------------------%
  %-prepare pipeline  
  cleanpipe = pipeline(...
     pset.node.physioset_import(...
    'Importer', pset.import.fieldtrip, ...
    'Verbose', cfg.gtool.verbose), ...
    ...
    pset.node.center(...
    'Save', cfg.gtool.saveall, ...
    'Verbose', cfg.gtool.verbose), ...
    ...
    pset.node.detrend(...
    'Save', cfg.gtool.saveall, ...
    'Verbose', cfg.gtool.verbose), ...
    ...
    filterNode1, ...
    ...
    pset.node.bad_samples(...
    'MADs', cfg.gtool.bad_samples.MADs, ...
    'Save', cfg.gtool.saveall, ...
    'Verbose', cfg.gtool.verbose), ...
    ...
    pset.node.bad_channels(...
    'MADs', cfg.gtool.bad_channels.MADs * [1 1], ... % bc there are two filters
    'Filter', {[], filter.hpfilt('fc', .5)}, ... % you can specify multiple filters
    'Save', cfg.gtool.saveall, ...
    'Verbose', cfg.gtool.verbose), ...
    ...
    pset.node.bss_regression.pwl(cfg.gtool.fsample, ...
    'PCA', opt.pcapwl, ...
    'BSS', opt.bsspwl, ...
    'Save', cfg.gtool.saveall, ...
    'Verbose', cfg.gtool.verbose), ...
    ...
    pset.node.bss_regression.eog(cfg.gtool.fsample, ...
    'BSS', opt.bssecg, ...
    'Correction', cfg.gtool.eog.correction, ...
    'Save', cfg.gtool.saveall, ...
    'Verbose', cfg.gtool.verbose), ...
    ...
    'OGE', false, ...
    'Save', false);

  %     pset.node.bss_regression.ecg(cfg.gtool.fsample, ...
  %     'BSS', opt.bssecg, ...
  %     'Save', cfg.gtool.saveall, ...
  %     'Verbose', cfg.gtool.verbose), ...
  %     ...

  %     pset.node.bss_regression.emg(cfg.gtool.fsample, ...
  %     'BSS', opt.bssecg, ...
  %     'Correction', cfg.gtool.emg.correction, ...
  %     'Save', cfg.gtool.saveall, ...
  %     'Verbose', cfg.gtool.verbose), ...
  
  %--------------------------%
  
  %--------------------------%
  %-run pipeline
  gdata = process(cleanpipe, [ddir allfile(i).name]);
  %--------------------------%
  
  %--------------------------%
  %-convert and prepare cfg
  data = fieldtrip(gdata);
  %--------------------------%
  
  %--------------------------%
  %-artifacts  
  %-----------------%
  %-reject trials
  if ~isempty(gdata.BadSampleIdx)
    badsmp = gdata.BadSampleIdx;
    bound = find(diff(badsmp)~=1);
    bound = [1 bound+1; bound numel(badsmp)]';
    artmat = badsmp(bound);
    artmat = artmat + data.sampleinfo(1);
  else
    artmat = [];
  end
  
  cfg1 = [];
  cfg1.artfctdef.gclean.artifact = artmat;
  cfg1.artfctdef.reject = 'partial';
  [data] = ft_rejectartifact(cfg1, data);
  %-----------------%
  
  %-----------------%
  %-repair channels
  cfg2 = [];
  cfg2.badchannel = data.label(gdata.BadChanIdx);
  if haselec
    cfg2.neighbours = neigh;
  elseif hasgrad
    cfg1 = [];
    cfg1.grad = grad;
    cfg1.method = 'distance';
    cfg1.neighbourdist = cfg.sens.dist;
    neigh = ft_prepare_neighbours(cfg1);
    cfg2.neighbours = neigh;
  end
  cfg2.feedback = 'none';
  
  [~, filename] = fileparts(allfile(i).name);
  cfg2.outputfile = [ddir filename '_' mfilename];
  
  ft_channelrepair(cfg2, data)
  %-----------------%
  %--------------------------%
  
  %--------------------------%
  %-save event
  load([ddir allfile(i).name], 'event')
  save(cfg2.outputfile, 'event', '-append')
  %--------------------------%
  
  %--------------------------%
  %-output
  outtmp = sprintf('%s ', data.label{gdata.BadChanIdx});
  output = sprintf('%s   Bad channels: %s \n', output, outtmp);
  
  output = sprintf('%s   N bad samples:% 4.f\n', ...
    output, numel(gdata.BadSampleIdx));
  %--------------------------%
  
  %-----------------%
  %-clear
  previousstep = cfg.step{find(strcmp(cfg.step, mfilename))-1};
  if any(strcmp(cfg.clear, previousstep))
    delete([ddir allfile(i).name])
  end
  %-----------------%
  
end
%------------------------------------%

%---------------------------%
%-end log
toc_t = toc(tic_t);
outtmp = sprintf('%s (%04d) ended at %s on %s after %s\n\n', ...
  mfilename, subj, datestr(now, 'HH:MM:SS'), datestr(now, 'dd-mmm-yy'), ...
  datestr( datenum(0, 0, 0, 0, 0, toc_t), 'HH:MM:SS'));
output = [output outtmp];

%-----------------%
fprintf(output)
fid = fopen([cfg.log '.txt'], 'a');
fwrite(fid, output);
fclose(fid);
%-----------------%
%---------------------------%
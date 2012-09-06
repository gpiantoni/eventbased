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
%  .gclean.fsample: manually specify the frequency (very easily bug-prone, but in this way it does not read "data" all the time)
%  .gclean.saveall: false
%  .gclean.verbose: true
%  .gclean.lpfreqn = [.3 / (cfg.gclean.fsample/2)]; % normalized by half of the sampling frequency!
%  .gclean.bad_samples.MADs = 5;
%  .gclean.bad_channels.MADs = 8;
%  .gclean.bad_samples.Percentile = [25 75];
%  .gclean.eog.correction = 50;
%  .gclean.emg.correction = 30;
%  .gclean.pwl.pca: if not empty, number of PC to work with in JADE
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
import pset.node.*;

for i = 1:numel(allfile)
  
  %-----------------%
  %-load the data
  outtmp = sprintf('%s\n', allfile(i).name);
  output = [output outtmp];
  %-----------------%
  
  
  %--------------------------%
  %-German toolbox parameters
  %------%
  %-PCA
  if isfield(cfg.gclean, 'pwl') && ...
      isfield(cfg.gclean.pwl, 'pca')
    opt.pcapwl = spt.pca('MaxDimOut', cfg.gclean.pwl.pca);
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
  if ~isempty(cfg.gclean.lpfreqn)
    filterNode1 = ...
      bpfilt('fp', [cfg.gclean.lpfreqn 1], ...
      'save', cfg.gclean.saveall);
  else
    filterNode1 = [];
  end
  %------%
  %--------------------------%
  
  %--------------------------%
  %  Construct the cleaning pipeline
  cleanpipe = pipeline(...
    ...
    ...  % Data importer
    physioset_import(...
    'Importer', pset.import.fieldtrip, ...
    'Verbose', cfg.gclean.verbose), ...
    ...
    ... % Remove data mean
    center(...
    'Save', cfg.gclean.saveall, ...
    'Verbose', cfg.gclean.verbose), ...
    ...
    ... % Detrend
    detrend(...
    'Save', cfg.gclean.saveall, ...
    'Verbose', cfg.gclean.verbose), ...
    ...
    ... % Bad channel rejection: using raw data and HP-filtered data
    bad_channels(...
    'MADs', cfg.gclean.bad_channels.MADs, ... % maybe [1 1] bc there are two filters
    'Filter', {[], filter.hpfilt('fc', 0.5)}, ... % you can specify multiple filters
    'Save', cfg.gclean.saveall, ...
    'Verbose', cfg.gclean.verbose), ...
    ...
    ... % Bad sample rejection
    bad_samples(...
    'MADs', cfg.gclean.bad_samples.MADs, ...
    'Percentile', cfg.gclean.bad_samples.Percentile, ...
    'Save', cfg.gclean.saveall, ...
    'Verbose', cfg.gclean.verbose), ...
    ...
    ... % HP filter
    filterNode1, ...
    ...
    ... % PWL removal
    bss_regression.pwl(cfg.gclean.fsample, ...
    'PCA', opt.pcapwl, ...
    'BSS', opt.bsspwl, ...
    'Chopper', chopper.dummy, ...
    'Save', cfg.gclean.saveall, ...
    'Verbose', cfg.gclean.verbose), ...
    ...
    ... % ECG removal
    bss_regression.ecg(cfg.gclean.fsample, ...
    'BSS', opt.bssecg, ...
    'Chopper', chopper.dummy, ...
    'Save', cfg.gclean.saveall, ...
    'Verbose', cfg.gclean.verbose), ...
    ...
    ... % EOG removal
    bss_regression.eog(cfg.gclean.fsample, ...
    'BSS', opt.bsseog, ...
    'Correction', cfg.gclean.eog.correction, ...
    'Chopper', chopper.dummy, ...
    'Save', cfg.gclean.saveall, ...
    'Verbose', cfg.gclean.verbose), ...
    ...
    ... % EMG removal
    bss_regression.emg(cfg.gclean.fsample, ...
    'BSS', opt.bssecg, ...
    'Correction', cfg.gclean.emg.correction, ...
    'Save', cfg.gclean.saveall, ...
    'Verbose', cfg.gclean.verbose), ...
    ...
    ... % Extra arguments to pipeline constructor
    'OGE', false, ...
    'Save', false);
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
  if any(gdata.BadSample)
    badsmp = find(gdata.BadSample);
    bound = find(diff(badsmp)~=1);
    bound = [1 bound+1; bound numel(badsmp)]';
    artmat = badsmp(bound);
    artmat = artmat + data.sampleinfo(1);
  else
    artmat = [];
    badsmp = 0;
  end
  
  cfg1 = [];
  cfg1.artfctdef.gclean.artifact = artmat;
  cfg1.artfctdef.reject = 'partial';
  [data] = ft_rejectartifact(cfg1, data);
  %-----------------%
  
  %-----------------%
  %-repair channels
  cfg2 = [];
  cfg2.badchannel = data.label(gdata.BadChan);
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
  outtmp = sprintf('%s ', data.label{gdata.BadChan});
  output = sprintf('%s   Bad channels: %s \n', output, outtmp);
  
  output = sprintf('%s   N bad samples: %0.3f percent\n', ...
    output, numel(badsmp)/length(gdata));
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
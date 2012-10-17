function gclean(info, opt, subj)
%GCLEAN use German's toolbox to clean the data
% Read the data from seldata (it assumes that the data is continuous), low
% pass filter, reject bad channels and eye, ecg, emg activity with ICA.
% It returns a fieldtrip structure made of trials of different length only
% containing good data. Bad channels are interpolated.
%
% INFO
%  .data: path of /data1/projects/PROJ/subjects/
%  .nick: NICK in /data1/projects/PROJ/subjects/0001/MOD/NICK/
%  .mod: modality, MOD in /data1/projects/PROJ/subjects/0001/MOD/NICK/
%
%  .log: name of the file and directory to save log
%
%  .sens.file: file with EEG sensors. It can be sfp or mat
%  .sens.dist: distance between sensors to consider them neighbors (in the units of cfg.sens.file)
%
% CFG.OPT
%  .fsample: manually specify the frequency (very easily bug-prone, but in this way it does not read "data" all the time)
%  .lpfreqn = [.3 / (opt.fsample/2)]; % normalized by half of the sampling frequency!
%  .bad_samples.MADs = 5;
%  .bad_channels.MADs = 8;
%  .bad_samples.Percentile = [25 75];
%  .eog.correction = 50;
%  .emg.correction = 30;
%  .pwl.pca: if not empty, number of PC to work with in JADE
%
% IN
%  data in /data1/projects/PROJ/subjects/SUBJ/MOD/NICK/
%
% OUT
%  data, after ICA rejection for eyes-blinks and movement, rejection and
%  interpolation of bad electrodes, rejection of bad epochs
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
ddir = sprintf('%s%04d/%s/%s/', info.data, subj, info.mod, info.nick); % data dir
allfile = dir([ddir '*_A.mat']); % files matching a preprocessing

prepr_name = '_B'; % preprocessing name to append
%---------------------------%

%---------------------------%
%-elec or grad
haselec = false;
if isfield(info.sens, 'file') && ~isempty(info.sens.file)
  haselec = true;
  sens = ft_read_sens(info.sens.file);
  sens.label = upper(sens.label);
  
  cfg = [];
  cfg.elec = sens;
  cfg.method = 'distance';
  cfg.neighbourdist = info.sens.dist;
  neigh = ft_prepare_neighbours(cfg);
end

if strcmp(info.mod, 'meg')
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
  if isfield(opt, 'pwl') && ...
      isfield(opt.pwl, 'pca')
    opt.pcapwl = spt.pca('MaxDimOut', opt.pwl.pca);
  else
    opt.pcapwl = [];
  end
  %------%
  
  %------%
  %-JADE
  opt.bsspwl = spt.jade;
  opt.bsseog = spt.jade;
  opt.bssecg = spt.jade;
  opt.bssemg = spt.jade;
  opt.bssrs  = spt.jade;
  %------%
  %--------------------------%
  
  %--------------------------%
  %  Construct the cleaning pipeline
  cleanpipe = pipeline(...
    ...
    ...  % Data importer
    physioset_import(...
    'Importer', pset.import.fieldtrip, ...
    'verbose', true), ...
    ...
    center('verbose', true), ... % Remove data mean
    ...
    detrend('verbose', true), ... % Detrend
    ...
    bad_channels(... % Bad channel rejection: using raw data and HP-filtered data
    'MADs', opt.bad_channels.MADs, ...
    'Filter', {[], filter.hpfilt('fc', 0.5)}, ... % you can specify multiple filters
    'verbose', true), ... 
    ...
    bad_samples(... % Bad sample rejection
    'MADs', opt.bad_samples.MADs, ...
    'Percentile', opt.bad_samples.Percentile, ...
    'verbose', true), ...
    ...
    bpfilt('fp', [opt.hpfreq 1], ...
    'verbose', true), ... % highpass filter
    ...
    bss_regression.pwl(opt.fsample, ... % PWL removal
    'PCA', opt.pcapwl, ...
    'BSS', opt.bsspwl, ...
    'Chopper', chopper.dummy, ...
    'verbose', true), ...
    ...
    ... % ECG removal
    bss_regression.ecg(opt.fsample, ...
    'BSS', opt.bssecg, ...
    'Chopper', chopper.dummy, ...
    'verbose', true), ...
    ...
    ... % EOG removal
    bss_regression.eog(opt.fsample, ...
    'BSS', opt.bsseog, ...
    'Correction', opt.eog.correction, ...
    'Chopper', chopper.dummy, ...
    'verbose', true), ...
    ...
    ... % EMG removal
    bss_regression.emg(opt.fsample, ...
    'BSS', opt.bssecg, ...
    'Correction', opt.emg.correction, ...
    'verbose', true), ...
    ...
    ... % Extra arguments to pipeline constructor
    'OGE', true, ...
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
  
  cfg = [];
  cfg.artfctdef.gclean.artifact = artmat;
  cfg.artfctdef.reject = 'partial';
  [data] = ft_rejectartifact(cfg, data);
  %-----------------%
  
  %-----------------%
  %-repair channels
  cfg = [];
  cfg.badchannel = data.label(gdata.BadChan);
  if haselec
    cfg.neighbours = neigh;
  elseif hasgrad
    tmpcfg = [];
    tmpcfg.grad = grad;
    tmpcfg.method = 'distance';
    tmpcfg.neighbourdist = info.sens.dist;
    neigh = ft_prepare_neighbours(tmpcfg);
    cfg.neighbours = neigh;
  end
  cfg.feedback = 'none';
  
  [~, filename] = fileparts(allfile(i).name);
  cfg.outputfile = [ddir filename prepr_name];
  
  ft_channelrepair(cfg, data)
  %-----------------%
  %--------------------------%
  
  %--------------------------%
  %-save event
  load([ddir allfile(i).name], 'event')
  save(cfg.outputfile, 'event', '-append')
  %--------------------------%
  
  %--------------------------%
  %-output
  outtmp = sprintf('%s ', data.label{gdata.BadChan});
  output = sprintf('%s   Bad channels: %s \n', output, outtmp);
  
  output = sprintf('%s   N bad samples: %0.3f percent\n', ...
    output, numel(badsmp)/length(gdata));
  %--------------------------%
  
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
fid = fopen([info.log '.txt'], 'a');
fwrite(fid, output);
fclose(fid);
%-----------------%
%---------------------------%
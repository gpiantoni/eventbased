function gclean(cfg, subj)
%PREPROCESS DATA USING GDEV TOOLBOX
% use German's toolbox to clean the data

%TODO: make it parallel over files (I think the gtool can already do it,
%only fixing the fields in data should not be parallel, but it should be
%very fast

mversion = 6;
%06 12/02/08 added emg and multiple filter for channel rejection
%05 12/02/02 do interpolation, but keep track of bad channels and bad samples
%04 12/02/02 using gtoolbox now (it keeps all the fields)
%03 12/02/01 reject trials only if there are artifacts
%02 12/02/01 fixed bug using i twice for two loops
%01 12/01/31 created

%-----------------%
%-input
if nargin == 1
  subj = cfg.subj;
end
%-----------------%

%---------------------------%
%-start log
output = sprintf('(p%02.f) %s (v%02.f) started at %s on %s\n', ...
  subj, mfilename,  mversion, datestr(now, 'HH:MM:SS'), datestr(now, 'dd-mmm-yy'));
tic_t = tic;
%---------------------------%

%---------------------------%
%-dir and files
ddir = sprintf('%s%04.f/%s/%s/', cfg.data, subj, cfg.mod, cfg.cond); % data
allfile = dir([ddir '*' cfg.endname '.mat']); % files matching a preprocessing

import aar.node.pipeline

sens = ft_read_sens(cfg.sens.file);
sens.label = upper(sens.label);

cfg1 = [];
cfg1.elec = sens;
cfg1.method = 'distance';
cfg1.neighbourdist = cfg.sens.dist;
neigh = ft_prepare_neighbours(cfg1);
%---------------------------%

%------------------------------------%
%-loop over files
for i = 1:numel(allfile) % this can run in parallel
  
  %-----------------%
  %-load the data
  outtmp = sprintf('%s\n', allfile(i).name);
  output = [output outtmp];
  %-----------------%
  
  %--------------------------%
  %-gtoolbox parameters
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
      aar.node.bpfilt('fp', [cfg.gtool.lpfreqn 1], ...
      'save', cfg.gtool.saveall);
  else
    filterNode1 = [];
  end
  %------%
  %--------------------------%
  
  %--------------------------%
  %-prepare pipeline  
  cleanpipe = pipeline('Node', ...
    {...
    aar.node.eegset_import(...
    'Importer', pset.import.fieldtrip, ...
    'Verbose', cfg.gtool.verbose), ...
    ...
    aar.node.center(...
    'Save', cfg.gtool.saveall, ...
    'Verbose', cfg.gtool.verbose), ...
    ...
    aar.node.detrend(...
    'Save', cfg.gtool.saveall, ...
    'Verbose', cfg.gtool.verbose), ...
    ...
    filterNode1, ...
    ...
    aar.node.bad_samples(...
    'MADs', cfg.gtool.bad_samples.MADs, ...
    'Percentile', cfg.gtool.bad_samples.perc, ...
    'Save', cfg.gtool.saveall, ...
    'Verbose', cfg.gtool.verbose), ...
    ...
    aar.node.bad_channels(...
    'MADs', cfg.gtool.bad_channels.MADs * [1 1], ... % bc there are two filters
    'Filter', {[], filter.hpfilt('fc', .5)}, ... % you can specify multiple filters
    'Save', cfg.gtool.saveall, ...
    'Verbose', cfg.gtool.verbose), ...
    ...
    aar.node.bss_regression.pwl(cfg.gtool.fsample, ...
    'BSS', opt.bsspwl, ...
    'Save', cfg.gtool.saveall, ...
    'Verbose', cfg.gtool.verbose), ...
    ...
    aar.node.bss_regression.ecg(cfg.gtool.fsample, ...
    'BSS', opt.bssecg, ...
    'Save', cfg.gtool.saveall, ...
    'Verbose', cfg.gtool.verbose), ...
    ...
    aar.node.bss_regression.eog(cfg.gtool.fsample, ...
    'BSS', opt.bssecg, ...
    'Correction', cfg.gtool.eog.correction, ...
    'Save', cfg.gtool.saveall, ...
    'Verbose', cfg.gtool.verbose), ...
    ...
    aar.node.bss_regression.emg(cfg.gtool.fsample, ...
    'BSS', opt.bssecg, ...
    'Correction', cfg.gtool.emg.correction, ...
    'Save', cfg.gtool.saveall, ...
    'Verbose', cfg.gtool.verbose), ...
    }, ...
    'OGE', cfg.gtool.oge, ...
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
  if ~isempty(gdata.BadSampleIdx)
    badsmp = gdata.BadSampleIdx;
    bound = find(diff(badsmp)~=1);
    bound = [1 bound+1; bound numel(badsmp)]';
    artmat = badsmp(bound);
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
  cfg2.neighbours = neigh;
  cfg2.feedback = 'none';
  
  [~, filename] = fileparts(allfile(i).name);
  cfg2.outputfile = [ddir filename '_' mfilename];
  
  ft_channelrepair(cfg2, data)
  %-----------------%
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
  if any(strcmp(mfilename, cfg.step(cfg.clear+1)))
    delete([ddir allfile(i).name])
  end
  %-----------------%
  
end
%------------------------------------%

%---------------------------%
%-end log
toc_t = toc(tic_t);
outtmp = sprintf('(p%02.f) %s (v%02.f) ended at %s on %s after %s\n\n', ...
  subj, mfilename, mversion, datestr(now, 'HH:MM:SS'), datestr(now, 'dd-mmm-yy'), ...
  datestr( datenum(0, 0, 0, 0, 0, toc_t), 'HH:MM:SS'));
output = [output outtmp];

%-----------------%
fprintf(output)
fid = fopen([cfg.log '.txt'], 'a');
fwrite(fid, output);
fclose(fid);
%-----------------%
%---------------------------%
function preproc(cfg, subj)
%PREPROC do preprocessing (including rerefering) on all data
%
% CFG
%  .data: name of projects/PROJNAME/subjects/
%  .mod: name of the modality used in recordings and projects
%  .cond: name to be used in projects/PROJNAME/subjects/0001/MOD/CONDNAME/
%  .endname: includes previous steps '_seldata_gclean'
%  .log: name of the file and directory with analysis log
%
%  .sens.file: file with EEG sensors. It can be sfp or mat
%  .sens.dist: distance between sensors to consider them neighbors (in the units of cfg.sens.file)
%
%  .step: all the analysis step (for cfg.clear)
%  .clear: index of cfg.step to remove from subject directory
%
%  .preproc: anything that goes into ft_preprocessing(cfg.preproc, data)
%  Suggested:
%     .preproc.reref = 'yes';
%     .preproc.refchannel = 'all';
%     .preproc.implicit = [];
%     .preproc.hpfilter = 'yes';
%     .preproc.hpfreq = 0.5;
%     .preproc.hpfiltord = 4;
%
%  .csd.do: do current source density ('yes' or 'no')
%  .csd.method: which method for CSD ('finite' or 'spline' or 'hjorth')
%
% Part of EVENTBASED preprocessing
% see also SELDATA, GCLEAN, PREPROC, REDEF

%---------------------------%
%-start log
output = sprintf('(p%02.f) %s started at %s on %s\n', ...
  subj, mfilename,  datestr(now, 'HH:MM:SS'), datestr(now, 'dd-mmm-yy'));
tic_t = tic;
%---------------------------%

%---------------------------%
%-dir and files
ddir = sprintf('%s%04.f/%s/%s/', cfg.data, subj, cfg.mod, cfg.cond); % data
allfile = dir([ddir '*' cfg.endname '.mat']); % files matching a preprocessing

sens = ft_read_sens(cfg.sens.file);
sens.label = upper(sens.label);
%---------------------------%

%---------------------------%
%-loop over files
for i = 1:numel(allfile)

  %-----------------%
  %-files
  inputfile = [ddir allfile(i).name];
  [~, filename] = fileparts(allfile(i).name);
  outputfile = [ddir filename '_' mfilename];
  %-----------------%

  %-----------------%
  %-prepr
  load(inputfile); % ft_preprocessing does not really use inputfile
  cfg1 = cfg.preproc;
  cfg1.feedback = 'none';
  cfg1.outputfile = outputfile;
  ft_preprocessing(cfg1, data);
  save(outputfile, 'event', '-append')
  %-----------------%

  %-----------------%
  %-scalp current density
  if strcmpi(cfg.csd.do, 'yes')
    cfg1 = [];
    cfg1.method = cfg.csd.method;
    cfg1.elec = sens;
    cfg1.feedback = 'none';
    cfg1.inputfile = outputfile; % it rewrites the same file
    cfg1.outputfile = outputfile;
    ft_scalpcurrentdensity(cfg1);
  end
  %-----------------%
  
  %-----------------%
  %-clear
  if any(strcmp(mfilename, cfg.step(cfg.clear+1)))
    delete([ddir allfile(i).name])
  end
  %-----------------%

end
%---------------------------%

%---------------------------%
%-end log
toc_t = toc(tic_t);
outtmp = sprintf('(p%02.f) %s ended at %s on %s after %s\n\n', ...
  subj, mfilename, datestr(now, 'HH:MM:SS'), datestr(now, 'dd-mmm-yy'), ...
  datestr( datenum(0, 0, 0, 0, 0, toc_t), 'HH:MM:SS'));
output = [output outtmp];

%-----------------%
fprintf(output)
fid = fopen([cfg.log '.txt'], 'a');
fwrite(fid, output);
fclose(fid);
%-----------------%
%---------------------------%
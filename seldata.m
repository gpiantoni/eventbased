function seldata(cfg, subj)
%SELDATA get data from recordings and put them in subject directory
% it recreates the "cfg.cond" folder for each subject
%
% CFG
%  .recs: name of recordings/RECNAME/subjects/
%  .mod: name of the modality used in recordings and projects
%  .rcnd: specific name of the condition of interest in the raw recording folder
%  .data: name of projects/PROJNAME/subjects/
%  .cond: name to be used in projects/PROJNAME/subjects/0001/MOD/CONDNAME/
%  .log: name of the file and directory with analysis log
% 
%  .sens.file: file with EEG sensors. It can be sfp or mat.
% 
%  .seldata.trialfun: name of the trialfun used to read the data before preprocessing, the function should be in PROJNAME_private/
%  .seldata.selchan: channels to read. It can be a vector or a cell of strings with the elec names on file (Micromed elec names are '  1' '  2'  '  3'
%  .seldata.label: if not empty, labels of electrodes to rename (same length as seldata.selchan)
%
% You need to write your own function to read the data and the events.
% Call the function something like "trialfun_XXX" and use as:
%   [trl, event] = trialfun_XXX(cfg)
% where trl is 1x3 vector (as in ft_definetrial) and event is the structure
% which can be used later on in redef.m to prepare the actual trials.
% It's better if you prepare only one big trials. gclean will clean the
% whole trial and prefers continuous data. You can create smaller trials
% later, during redef.m
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
rdir = sprintf('%s%04.f/%s/%s/', cfg.recs, subj, cfg.mod, 'raw'); % recording raw
ddir = sprintf('%s%04.f/%s/%s/', cfg.data, subj, cfg.mod, cfg.cond); % data
if isdir(ddir); rmdir(ddir, 's'); end
mkdir(ddir)

%-----------------%
%-read sensors (can be sfp or mat)
sens = ft_read_sens(cfg.sens.file);
sens.label = upper(sens.label); % <- EGI labels are uppercase, but the elec file is lowercase
%-----------------%
%---------------------------%

%---------------------------%
%-find raw data
allfile = dir([rdir '*' cfg.rcnd '*']);

for i = 1:numel(allfile)
  
  dataset = [rdir allfile(i).name];
  
  %-----------------%
  %-definetrials
  cfg1 = [];
  cfg1.trialfun = cfg.seldata.trialfun;
  cfg1.dataset = dataset;
  
  cfg2 = ft_definetrial(cfg1);
  
  %--------%
  %-ignore files with no events
  if all(cfg2.trl(1,1:3) == [0 0 0])
    continue
  end
  %--------%
  %-----------------%
    
  %-----------------%
  %-preprocessing
  cfg2.feedback = 'off';
  if iscell(cfg.seldata.selchan)
    cfg2.channel = cfg.seldata.selchan{subj};
  else
    cfg2.channel = cfg.seldata.selchan;
  end
  
  data = ft_preprocessing(cfg2);
  %-----------------%
  
  %-----------------%
  %-fix channels
  if ~isempty(cfg.seldata.label)
    data.label = cfg.seldata.label;
  end
  data.elec = sens;
  %-----------------%
  
  %-----------------%
  %-save data
  [~, filename] = fileparts(allfile(i).name);
  savename = [cfg.proj '_' filename '_' mfilename]; % <-- add proj name
  save([ddir savename], 'data');
  clear data
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
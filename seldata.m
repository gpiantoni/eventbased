function seldata(info, opt, subj)
%SELDATA get data from recordings and put them in subject directory
% it recreates the "info.nick" folder for each subject
%
% INFO
%  .recs: path of /data1/projects/PROJ/recordings/REC/subjects/
%  .mod: modality, MOD in /data1/projects/PROJ/subjects/SUBJ/MOD/NICK/
%
%  .data: name of projects/PROJ/subjects/
%  .nick: NICK in /data1/projects/PROJ/subjects/SUBJ/MOD/NICK/
%  .log: name of the file and directory to save log
% 
%  .sens.file: file with EEG sensors. It can be sfp or mat. It's included
%              in data struct. If empty, it does not read the sensors.
%
% CFG.OPT
%  .rcnd: specific name of the condition of interest in the raw recording folder
% 
%  .trialfun: name of the trialfun used to read the data, see below.
%                     The function should be in NICK_private/
%  .trialopt: options to pass to trialfun_XXX
%  .selchan: channels to read. It can be a vector or a cell of
%                    strings with the elec names on file (Micromed elec
%                    names are '  1' '  2'  '  3') 
%  .label: if not empty, labels of electrodes to rename (same
%                  length as cfg.opt.selchan) 
%
% IN
%  raw data in any format Fieldtrip can read in the recording folder:
%  /data1/projects/PROJ/recordings/REC/subjects/SUBJ/MOD/RAW/
% 
% OUT
%  data in FieldTrip format in /data1/projects/PROJ/subjects/SUBJ/MOD/NICK/
%
% You need to write your own function to read the data and the events.
% Call the function something like "trialfun_XXX" and use as:
%   [trl, event] = trialfun_XXX(opt)
% where trl is 1x3 vector (as in ft_definetrial) and event is the structure
% which can be used later on in redef.m to prepare the actual trials.
% It's better if you prepare only one big trial. gclean will clean the
% whole trial and prefers continuous data. You can create smaller trials
% later, during redef.m
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
rdir = sprintf('%s%04d/%s/%s/', info.recs, subj, info.mod, 'raw'); % recording raw
ddir = sprintf('%s%04d/%s/%s/', info.data, subj, info.mod, info.nick); % data dir
if isdir(ddir); rmdir(ddir, 's'); end
mkdir(ddir)

%-----------------%
%-read sensors (can be sfp or mat)
hassens = false;
if isfield(info.sens, 'file') && ~isempty(info.sens.file)
  hassens = true;
  sens = ft_read_sens(info.sens.file);
  sens.label = upper(sens.label); % <- EGI labels are uppercase, but the elec file is lowercase
end
%-----------------%

prepr_name = '_A'; % preprocessing name to append
%---------------------------%

%---------------------------%
%-find raw data
allfile = dir([rdir '*' opt.rcnd '*']);

for i = 1:numel(allfile)
  
  dataset = [rdir allfile(i).name];
  
  %-----------------%
  %-definetrials
  if isfield(opt, 'trialcfg')
    cfg = opt.trialcfg;
  else
    cfg = [];
  end
  cfg.trialfun = opt.trialfun;
  cfg.dataset = dataset;
  
  cfg = ft_definetrial(cfg);
  
  %--------%
  %-ignore files with no events
  if all(cfg.trl(1,1:3) == [0 0 0])
    continue
  end
  %--------%
  %-----------------%
    
  %-----------------%
  %-preprocessing
  cfg.feedback = 'off';
  if iscell(opt.selchan)
    cfg.channel = opt.selchan{subj};
  else
    cfg.channel = opt.selchan;
  end
  
  cfg.continuous = 'yes'; % necessary for MEG data over trials
  data = ft_preprocessing(cfg);
  event = ft_findcfg(data.cfg, 'event');
  
  if ischar(event) && ...
      strcmp(event, 'empty - this was cleared by checkconfig')
    event = ft_read_event(event);
  end
  %-----------------%
  
  %-----------------%
  %-fix channels
  if ~isempty(opt.label)
    data.label = opt.label;
  end
  
  if hassens
    data.elec = sens;
  end
  %-----------------%
  
  %-----------------%
  %-save data
  [~, filename] = fileparts(allfile(i).name);
  savename = [info.nick '_' filename prepr_name]; % <-- add nick name
  save([ddir savename], 'data', 'event');
  clear data
  %-----------------%
 
end
%---------------------------%

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
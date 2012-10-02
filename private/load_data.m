function [data badchan] = load_data(info, subj, cond)
%LOAD_DATA load data, and optionally get bad channels
% Use as:
%   [data badchan] = load_data(info, subj, cond)
%
% INFO
%  .data: path of /data1/projects/PROJ/subjects/
%  .rec: REC in /data1/projects/PROJ/recordings/REC/
%  .nick: NICK in /data1/projects/PROJ/subjects/SUBJCODE/MOD/NICK/
%  .mod: modality, MOD in /data1/projects/PROJ/subjects/SUBJCODE/MOD/NICK/
%
% SUBJ
%   number indicating the subject number
%
% COND
%   a string with the name used to read the data. The file name is
%   structured as: REC_NICK_SUBJ_MOD_COND_endname
%   COND should not have leading or trailing underscores
%   So, you need to specify COND with astericks only if necessary
%
% DATA
%   data in the specified condition
% 
% BADCHANNEL
%   if gclean has been run, it'll return the labels of the channels
%   labeled as bad
% 
% Part of EVENTBASED/PRIVATE

%-----------------%
%-input and output for each condition
ddir = sprintf('%s%04d/%s/%s/', info.data, subj, info.mod, info.nick); % data
beginname = sprintf('%s_%s_%04d_%s_', info.nick, info.rec, subj, info.mod); % beginning of datafile

prepr_name = '_A_B_C';
dnames = dir([ddir beginname cond prepr_name '.mat']); % files matching a preprocessing
%-----------------%

%-----------------%
%-concatenate only if you have more datasets
if numel(dnames) > 1
  spcell = @(name) sprintf('%s%s', ddir, name);
  allname = cellfun(spcell, {dnames.name}, 'uni', 0);
  
  cfg = [];
  cfg.inputfile = allname;
  data = ft_appenddata(cfg);
  
elseif numel(dnames) == 1
  load([ddir dnames(1).name], 'data')
  
else
  data = [];
  
end
%-----------------%

%-----------------%
%-if empty
if isempty(data)
  badchan = [];
  return
end
%-----------------%

%-----------------%
%-find bad channels (bad channel if bad in at least one dataset)
if ~isfield(data.cfg, 'previous') || ~iscell(data.cfg.previous) % not appenddata
  badchan = ft_findcfg(data.cfg, 'badchannel');
  
else
  
  badchan = {};
  for i = 1:numel(data.cfg.previous)
    badtemp = ft_findcfg(data.cfg.previous{i}, 'badchannel');
    badchan = union(badchan, badtemp);
  end
  
end
%-----------------%
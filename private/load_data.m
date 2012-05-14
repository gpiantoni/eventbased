function [data badchan] = load_data(cfg, subj, cond)
%LOAD_DATA load data, and optionally get bad channels
% Use as:
%   [data badchan] = load_data(cfg, subj, cond)
%
% CFG
%  .data: path of /data1/projects/PROJ/subjects/
%  .rec: REC in /data1/projects/PROJ/recordings/REC/
%  .nick: NICK in /data1/projects/PROJ/subjects/SUBJCODE/MOD/NICK/
%  .mod: modality, MOD in /data1/projects/PROJ/subjects/SUBJCODE/MOD/NICK/
%  .endname: includes preprocessing steps (e.g. '_seldata_gclean_preproc_redef')
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
%   labelled as bad
% 
% Part of EVENTBASED/PRIVATE

%-----------------%
%-input and output for each condition
ddir = sprintf('%s%04d/%s/%s/', cfg.data, subj, cfg.mod, cfg.nick); % data
beginname = sprintf('%s_%s_%04d_%s_', cfg.nick, cfg.rec, subj, cfg.mod); % beginning of datafile

dnames = dir([ddir beginname cond cfg.endname '.mat']); % files matching a preprocessing
%-----------------%

%-----------------%
%-concatenate only if you have more datasets
if numel(dnames) > 1
  spcell = @(name) sprintf('%s%s', ddir, name);
  allname = cellfun(spcell, {dnames.name}, 'uni', 0);
  
  cfg1 = [];
  cfg1.inputfile = allname;
  data = ft_appenddata(cfg1);
  
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
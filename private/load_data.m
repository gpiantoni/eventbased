function [data badchan] = load_data(cfg, subj, cond)
%LOAD_DATA load data, and optionally get bad channels
% Use as:
%   [data badchan] = load_data(cfg, subj, cond)
%
% CFG
%  .data: path of /data1/projects/PROJNAME/subjects/
%  .rec: RECNAME in /data1/projects/PROJNAME/recordings/RECNAME/
%  .nick: NICKNAME in /data1/projects/PROJNAME/subjects/0001/MOD/NICKNAME/
%  .mod: modality, MOD in /data1/projects/PROJNAME/subjects/0001/MOD/NICKNAME/
%  .endname: includes preprocessing steps (e.g. '_seldata_gclean_preproc_redef')
%
% SUBJ
%   number indicating the subject number
%
% COND
%   a string with the name used to read the data. The file name is
%   structured as: RECNAME_NICKNAME_SUBJCODE_MOD_COND_endname
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
beginname = sprintf('%s_%s_%04d_%s_', cfg.rec, cfg.nick, subj, cfg.mod); % beginning of datafile

allfile = dir([ddir beginname cond cfg.endname '.mat']); % files matching a preprocessing
%-----------------%

%-----------------%
%-concatenate only if you have more datasets
if numel(allfile) > 1
  spcell = @(name) sprintf('%s%s', ddir, name);
  allname = cellfun(spcell, {allfile.name}, 'uni', 0);
  
  cfg1 = [];
  cfg1.inputfile = allname;
  data = ft_appenddata(cfg1);
  
elseif numel(allfile) == 1
  load([ddir allfile(1).name], 'data')
  
else
  data = [];
  
end
%-----------------%

%-----------------%
%-find bad channels (bad channel if bad in at least one dataset)
if ~iscell(data.cfg.previous) % not appenddata
  badchan = ft_findcfg(data.cfg, 'badchannel');
  
else
  
  badchan = {};
  for i = 1:numel(data.cfg.previous)
    badtemp = ft_findcfg(data.cfg.previous{i}, 'badchannel');
    badchan = union(badchan, badtemp);
  end
  
end
%-----------------%
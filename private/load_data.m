function [data badchan] = load_data(cfg, subj, cond)
%LOAD_DATA load data, and optionally get bad channels
% TODO: a better way to pass data to read
% 
% Part of EVENTBASED/PRIVATE

%-----------------%
%-input and output for each condition
ddir = sprintf('%s%04.f/%s/%s/', cfg.data, subj, cfg.mod, cfg.cond); % data

allfile = dir([ddir cond cfg.endname '.mat']); % files matching a preprocessing
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
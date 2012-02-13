function pow_subj(cfg, subj)
%POW_SUBJ check power changes at single-subject level

mversion = 5;
%05 12/02/02 renamed to pow_subj
%04 12/01/12 no csd at this point
%03 11/12/09 baseline can be empty
%02 11/11/20 fix wild card in DIR
%01 11/09/27 created from gosdconn

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
%---------------------------%

%-------------------------------------%
%-loop over conditions
for k = 1:numel(cfg.test)
  
  %-----------------%
  %-input and output for each condition
  allfile = dir([ddir '*' cfg.test{k} cfg.endname '.mat']); % files matching a preprocessing
  if isempty(allfile)
    continue
  end
  
  condname = regexprep(cfg.test{k}, '*', '');
  outputfile = sprintf('pow_%02.f_%s', subj, condname);
  %-----------------%
  
  %-----------------%
  %-concatenate only if you have more datasets
  if numel(allfile) > 1
    spcell = @(name) sprintf('%s%s', ddir, name);
    allname = cellfun(spcell, {allfile.name}, 'uni', 0);
    
    cfg1 = [];
    cfg1.inputfile = allname;
    data = ft_appenddata(cfg1);
    
  else
    load([ddir allfile(1).name], 'data')
    
  end
  %-----------------%

  %-----------------%
  %-calculate power
  cfg2 = cfg.pow;
  cfg2.feedback = 'none';
  data = ft_freqanalysis(cfg2, data);
  data.time = cfg2.toi;
  %-----------------%
  
  %-----------------%
  %-fix baseline
  if ~isempty(cfg.pow.bl.baseline)
    cfg3 = cfg.pow.bl;
    freq = ft_freqbaseline(cfg3, data);
    
  else
    freq = data;
    
  end
  %-----------------%
  save([cfg.dpow outputfile], 'freq')
  
end
%-------------------------------------%

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

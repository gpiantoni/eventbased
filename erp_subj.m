function erp_subj(cfg, subj)
%ERP_SUBJ create subject-specific erp

%---------------------------%
%-start log
output = sprintf('(p%02.f) %s started at %s on %s\n', ...
  subj, mfilename,  datestr(now, 'HH:MM:SS'), datestr(now, 'dd-mmm-yy'));
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
  outputfile = sprintf('erp_%02.f_%s', subj, condname);
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
  %-timelock analysis
  cfg2 = cfg.erp;
  cfg2.feedback = 'none';
  cfg2.preproc.feedback = 'none';
  cfg2.outputfile = [cfg.derp outputfile];
  ft_timelockanalysis(cfg2, data);
  %-----------------%
  
end
%-------------------------------------%

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

function erp_subj(cfg, subj)
%ERP_SUBJ create subject-specific erp
%
% CFG
%  .data: name of projects/PROJNAME/subjects/
%  .mod: name of the modality used in recordings and projects
%  .cond: name to be used in projects/PROJNAME/subjects/0001/MOD/CONDNAME/
%  .endname: includes previous steps '_seldata_gclean_preproc_redef'
%  .test: a cell with the condition defined by redef. This function will loop over cfg.test
%  .log: name of the file and directory with analysis log
%  .derp: directory to save ERP data
%
%  .erp: a structure with cfg to pass to ft_timelockanalysis
%
% OUT
%  [cfg.derp 'erp_001_TEST']: timelock analysis for single-subject
%
% Part of EVENTBASED single-subject
% see also ERP_SUBJ, ERP_GRAND, ERPSOURCE_SUBJ, ERPSOURCE_GRAND, 
% POW_SUBJ, POW_GRAND, POWSOURCE_SUBJ, POWSOURCE_GRAND, 
% POWCORR_SUBJ, POWCORR_GRAND,
% CONN_SUBJ, CONN_GRAND, CONN_STAT

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
  allfile = dir([ddir cfg.test{k} cfg.endname '.mat']); % files matching a preprocessing

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
    
  elseif numel(allfile) == 1
    load([ddir allfile(1).name], 'data')
    
  else
    output = sprintf('%sCould not find any file in %s for test %s\n', ...
      output, ddir, cfg.test{k});
    continue
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

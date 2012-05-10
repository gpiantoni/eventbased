function erp_subj(cfg, subj)
%ERP_SUBJ create subject-specific erp
%
% CFG
%  .data: path of /data1/projects/PROJNAME/subjects/
%  .rec: RECNAME in /data1/projects/PROJNAME/recordings/RECNAME/
%  .nick: NICKNAME in /data1/projects/PROJNAME/subjects/0001/MOD/NICKNAME/
%  .mod: modality, MOD in /data1/projects/PROJNAME/subjects/0001/MOD/NICKNAME/
%  .endname: includes preprocessing steps (e.g. '_seldata_gclean_preproc_redef')
%
%  .log: name of the file and directory to save log
%  .derp: directory for ERP data
%  .erp.cond: cell with conditions (e.g. {'*cond1' '*cond2'})'
%
%  .erp: a structure with cfg to pass to ft_timelockanalysis
%
% IN:
%  data in /PROJNAME/subjects/0001/MOD/NICKNAME/
% 
% OUT
%  [cfg.derp 'erp_SUBJCODE_CONDNAME']: timelock analysis for single-subject
%
% Part of EVENTBASED single-subject
% see also ERP_SUBJ, ERP_GRAND, ERPSOURCE_SUBJ, ERPSOURCE_GRAND, 
% POW_SUBJ, POW_GRAND, POWSOURCE_SUBJ, POWSOURCE_GRAND, 
% POWCORR_SUBJ, POWCORR_GRAND, POWSTAT_SUBJ, POWSTAT_GRAND, 
% CONN_SUBJ, CONN_GRAND, CONN_STAT

%---------------------------%
%-start log
output = sprintf('(p%02.f) %s started at %s on %s\n', ...
  subj, mfilename,  datestr(now, 'HH:MM:SS'), datestr(now, 'dd-mmm-yy'));
tic_t = tic;
%---------------------------%

%-------------------------------------%
%-loop over conditions
for k = 1:numel(cfg.erp.cond)
  cond     = cfg.erp.cond{k};
  condname = regexprep(cond, '*', '');
  
  %---------------------------%
  %-read data
  [data] = load_data(cfg, subj, cond);
  if isempty(data)
    output = sprintf('%sCould not find any file for condition %s\n', ...
      output, cond);
    continue
  end
  
  outputfile = sprintf('erp_%04d_%s', subj, condname);
  %---------------------------%
  
  %---------------------------%
  %-timelock analysis
  cfg2 = cfg.erp;
  cfg2.feedback = 'none';
  cfg2.preproc.feedback = 'none';
  cfg2.outputfile = [cfg.derp outputfile];
  ft_timelockanalysis(cfg2, data);
  %---------------------------%
  
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

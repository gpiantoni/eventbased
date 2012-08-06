function erp_subj(cfg, subj)
%ERP_SUBJ create subject-specific erp
%
% CFG
%  .data: path of /data1/projects/PROJ/subjects/
%  .rec: REC in /data1/projects/PROJ/recordings/REC/
%  .nick: NICK in /data1/projects/PROJ/subjects/0001/MOD/NICK/
%  .mod: modality, MOD in /data1/projects/PROJ/subjects/0001/MOD/NICK/
%  .endname: includes preprocessing steps (e.g. '_seldata_gclean_redef')
%
%  .log: name of the file and directory to save log
%  .derp: directory with ERP data
%  .erp.cond: cell with conditions (e.g. {'*cond1' '*cond2'})'
%  .erp.source: read virtual electrode data (logical)
%
%  .erp: a structure with cfg to pass to ft_timelockanalysis
%
% IN:
%  data in /PROJ/subjects/SUBJ/MOD/NICK/
% OR if cfg.erp.source
%  source in cfg.dsou from SOURCE_SUBJ
% 
% OUT
%  [cfg.derp 'erp_SUBJ_COND'] 'erp_s': timelock analysis for single-subject
%
% Part of EVENTBASED single-subject
% see also ERP_SUBJ, ERP_GRAND, 
% ERPSOURCE_SUBJ, ERPSOURCE_GRAND, ERPSTAT_SUBJ, ERPSTAT_GRAND,
% POW_SUBJ, POW_GRAND, POWCORR_SUBJ, POWCORR_GRAND,
% POWSOURCE_SUBJ, POWSOURCE_GRAND, POWSTAT_SUBJ, POWSTAT_GRAND,
% SOURCE_SUBJ, CONN_SUBJ, CONN_GRAND, CONN_STAT

%---------------------------%
%-start log
output = sprintf('%s (%04d) began at %s on %s\n', ...
  mfilename, subj, datestr(now, 'HH:MM:SS'), datestr(now, 'dd-mmm-yy'));
tic_t = tic;
%---------------------------%

%-------------------------------------%
%-loop over conditions
for k = 1:numel(cfg.erp.cond)
  cond     = cfg.erp.cond{k};
  condname = regexprep(cond, '*', '');
  
  %---------------------------%
  %-read data
  if ~isfield(cfg.erp, 'source') || cfg.erp.source
    [data] = load_data(cfg, subj, cond);
  else
    [data] = load_source(cfg, subj, cond);
  end
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
  erp_s = ft_timelockanalysis(cfg2, data);
  %---------------------------%
  
  save([cfg.derp outputfile], 'erp_s')
  
end
%-------------------------------------%

%---------------------------%
%-end log
toc_t = toc(tic_t);
outtmp = sprintf('%s (%04d) ended at %s on %s after %s\n\n', ...
  mfilename, subj, datestr(now, 'HH:MM:SS'), datestr(now, 'dd-mmm-yy'), ...
  datestr( datenum(0, 0, 0, 0, 0, toc_t), 'HH:MM:SS'));
output = [output outtmp];

%-----------------%
fprintf(output)
fid = fopen([cfg.log '.txt'], 'a');
fwrite(fid, output);
fclose(fid);
%-----------------%
%---------------------------%

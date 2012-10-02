function erp_subj(info, opt, subj)
%ERP_SUBJ timelock-analysis for each subject
%
% INFO
%  .derp: directory with ERP data
%  .log: name of the file and directory to save log
%
% CFG.OPT
%  .source: read virtual electrode data (logical)
%  .cond*: cell with conditions (e.g. {'*cond1' '*cond2'})'
%  .erp*: a structure with cfg to pass to ft_timelockanalysis
%
% IN
%  LOAD_DATA: data in /PROJ/subjects/SUBJ/MOD/NICK/
% OR if cfg.opt.source
%  LOAD_SOURCE: source in info.dsou after SOURCE_SUBJ
% 
% OUT
%  [info.derp 'erp_SUBJ_COND'] 'erp_s': timelock analysis for single-subject
%
% * indicates obligatory parameter
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
for k = 1:numel(opt.cond)
  cond = opt.cond{k};
  condname = regexprep(cond, '*', '');
  
  %---------------------------%
  %-read data
  if ~isfield(opt, 'source') || ~info.source
    [data] = load_data(info, subj, cond);
  else
    [data] = load_source(info, subj, cond);
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
  cfg = opt.erp;
  cfg.feedback = 'none';
  cfg.preproc.feedback = 'none';
  erp_s = ft_timelockanalysis(cfg, data);
  %---------------------------%
  
  save([info.derp outputfile], 'erp_s')
  
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
fid = fopen([info.log '.txt'], 'a');
fwrite(fid, output);
fclose(fid);
%-----------------%
%---------------------------%

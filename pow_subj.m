function pow_subj(cfg, subj)
%POW_SUBJ create subject-specific pow
%
% CFG
%  .data: path of /data1/projects/PROJNAME/subjects/
%  .rec: RECNAME in /data1/projects/PROJNAME/recordings/RECNAME/
%  .nick: NICKNAME in /data1/projects/PROJNAME/subjects/0001/MOD/NICKNAME/
%  .mod: modality, MOD in /data1/projects/PROJNAME/subjects/0001/MOD/NICKNAME/
%  .endname: includes preprocessing steps (e.g. '_seldata_gclean_preproc_redef')
%
%  .log: name of the file and directory to save log
%  .dpow: directory for POW data
%  .pow.cond: cell with conditions (e.g. {'*cond1' '*cond2'})'
%
%  .pow: a structure with cfg to pass to ft_freqanalysis
%
% Baseline correction is applied in POW_GRAND
%
% IN:
%  data in /PROJNAME/subjects/0001/MOD/NICKNAME/
%
% OUT
%  [cfg.dpow 'pow_SUBJCODE_COND']: power analysis for single-subject
%
% Part of EVENTBASED single-subject
% see also ERP_SUBJ, ERP_GRAND, ERPSOURCE_SUBJ, ERPSOURCE_GRAND, 
% POW_SUBJ, POW_GRAND, POWSOURCE_SUBJ, POWSOURCE_GRAND, 
% POWCORR_SUBJ, POWCORR_GRAND, POWSTAT_SUBJ, POWSTAT_GRAND, 
% CONN_SUBJ, CONN_GRAND, CONN_STAT

%---------------------------%
%-start log
output = sprintf('(p%04d) %s started at %s on %s\n', ...
  subj, mfilename,  datestr(now, 'HH:MM:SS'), datestr(now, 'dd-mmm-yy'));
tic_t = tic;
%---------------------------%

%-------------------------------------%
%-loop over conditions
for k = 1:numel(cfg.pow.cond)
  cond     = cfg.pow.cond{k};
  condname = regexprep(cond, '*', '');
  
  %---------------------------%
  %-read data
  [data] = load_data(cfg, subj, cond);
  if isempty(data)
    output = sprintf('%sCould not find any file for condition %s\n', ...
      output, cond);
    continue
  end
  
  outputfile = sprintf('pow_%04d_%s', subj, condname);
  %---------------------------%
  
  %---------------------------%
  %-calculate power
  cfg2 = cfg.pow;
  cfg2.feedback = 'etf';
  
  if cfg.pow.outliers
    cfg2.keeptrials = 'yes';
  end
  
  freq = ft_freqanalysis(cfg2, data);
  
  if isfield(cfg.pow, 'toi')
    freq.time = cfg.pow.toi;
  end
  %---------------------------%
  
  save([cfg.dpow outputfile], 'freq')
  
end
%-------------------------------------%

%---------------------------%
%-end log
toc_t = toc(tic_t);
outtmp = sprintf('(p%04d) %s ended at %s on %s after %s\n\n', ...
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

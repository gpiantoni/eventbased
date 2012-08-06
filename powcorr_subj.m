function powcorr_subj(cfg, subj)
%POWCORR_SUBJ correlate power with trialinfo
%
% CFG
%  .data: path of /data1/projects/PROJ/subjects/
%  .rec: REC in /data1/projects/PROJ/recordings/REC/
%  .nick: NICK in /data1/projects/PROJ/subjects/0001/MOD/NICK/
%  .mod: modality, MOD in /data1/projects/PROJ/subjects/0001/MOD/NICK/
%  .endname: includes preprocessing steps (e.g. '_seldata_gclean_redef')
%
%  .log: name of the file and directory to save log
%  .dpow: directory for POW data
%  .powcorr.cond: cell with conditions (e.g. {'*cond1' '*cond2'})'
%  .powcorr.source: read virtual electrode data (logical)
%
%  .powcorr: a structure with cfg to pass to ft_freqanalysis
%  .powcorr.bl.baseline: two scalars with baseline windows (if empty, no baseline)
%
%  .powcorr.info: column of trialinfo to use for the correlation
%  .powcorr.log: logical (take the log of power, strongly advised)
%
% IN:
%  data in /PROJ/subjects/SUBJ/MOD/NICK/
% OR if cfg.erp.source
%  source in cfg.dsou from SOURCE_SUBJ
%
% OUT
%  [cfg.dpow 'powcorr_SUBJ_COND'] 'powcorr_s': power correlation for single-subject
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
for k = 1:numel(cfg.powcorr.cond)
  cond     = cfg.powcorr.cond{k};
  condname = regexprep(cond, '*', '');
  
  %---------------------------%
  %-read data
  if ~isfield(cfg.powcorr, 'source') || cfg.powcorr.source
    [data] = load_data(cfg, subj, cond);
  else
    [data] = load_source(cfg, subj, cond);
  end
  if isempty(data)
    output = sprintf('%sCould not find any file for condition %s\n', ...
      output, cond);
    continue
  end
  
  outputfile = sprintf('powcorr_%04d_%s', subj, condname);
  %---------------------------%
  
  %---------------------------%
  %-calculate power
  cfg2 = cfg.powcorr;
  cfg2.feedback = 'none';
  cfg2.keeptrials = 'yes';
  data = ft_freqanalysis(cfg2, data);
  
  if isfield(cfg.powcorr, 'toi')
    data.time = cfg.powcorr.toi;
  end
  %---------------------------%
  
  %---------------------------%
  %-fix baseline
  if isfield(cfg.powcorr, 'bl') && ~isempty(cfg.powcorr.bl.baseline)
    cfg3 = cfg.powcorr.bl;
    powcorr_s = ft_freqbaseline(cfg3, data);
    
  else
    powcorr_s = data;
    
  end
  %---------------------------%
  
  %---------------------------%
  %-regression at each point
  if cfg.powcorr.log
    powcorr_s.powspctrm = log(powcorr_s.powspctrm);
  end
  
  [s1 s2 s3 s4] = size(powcorr_s.powspctrm);
  
  powspctrm = nan(s2,s3,s4);
  warning('off', 'MATLAB:rankDeficientMatrix') % NaN in \
  
  for i2 = 1:s2
    for i3 = 1:s3
      for i4 = 1:s4
        regr = [ones(size(powcorr_s.trialinfo,1),1) powcorr_s.powspctrm(:,i2,i3,i4)];
        beta = regr \ powcorr_s.trialinfo(:, cfg.powcorr.info);
        powspctrm(i2,i3,i4) = beta(2);
      end
    end
  end
  warning('on', 'MATLAB:rankDeficientMatrix')
  %---------------------------%
  
  %---------------------------%
  %-restructure freq for multiple subjects
  powcorr_s = rmfield(powcorr_s, 'trialinfo');
  powcorr_s.dimord = powcorr_s.dimord(5:end);
  powcorr_s.powspctrm = powspctrm;
  %---------------------------%
  
  save([cfg.dpow outputfile], 'powcorr_s')
  
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

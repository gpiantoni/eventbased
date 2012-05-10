function powcorr_subj(cfg, subj)
%POWCORR_SUBJ correlate power with trialinfo
%
% CFG
%  .data: name of projects/PROJNAME/subjects/
%  .mod: name of the modality used in recordings and projects
%  .cond: name to be used in projects/PROJNAME/subjects/0001/MOD/CONDNAME/
%  .endname: includes previous steps '_seldata_gclean_preproc_redef'
%  .log: name of the file and directory with analysis log
%  .test: a cell with the condition defined by redef. This function will loop over cfg.test
%  .dpow: directory to save ERP data
%
%  .powcorr: a structure with cfg to pass to ft_freqanalysis
%  .powcorr.bl.baseline: two scalars with baseline windows (if empty, no baseline)
%
%  .powcorr.info: column of trialinfo to use for the correlation
%  .powcorr.log: logical (take the log of power, strongly advised)
%
% OUT
%  [cfg.dpow 'powcorr_001_TEST']: power analysis for single-subject
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

%-------------------------------------%
%-loop over conditions
for k = 1:numel(cfg.powcorr.cond) % DOC: CFG.POWCORR.COND
  cond     = cfg.powcorr.cond{k};
  condname = regexprep(cond, '*', '');
  
  %---------------------------%
  %-read data
  [data] = load_data(cfg, subj, cond);
  if isempty(data)
    output = sprintf('%sCould not find any file for condition %s\n', ...
      output, cond);
    continue
  end
  
  outputfile = sprintf('powcorr_%02.f_%s', subj, condname);
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
  if ~isempty(cfg.powcorr.bl.baseline)
    cfg3 = cfg.powcorr.bl;
    freq = ft_freqbaseline(cfg3, data);
    
  else
    freq = data;
    
  end
  %---------------------------%
  
  %---------------------------%
  %-regression at each point
  if cfg.powcorr.log
    freq.powspctrm = log(freq.powspctrm);
  end
  
  [s1 s2 s3 s4] = size(freq.powspctrm);
  
  powspctrm = nan(s2,s3,s4);
  warning('off', 'MATLAB:rankDeficientMatrix') % NaN in \
  
  for i2 = 1:s2
    for i3 = 1:s3
      for i4 = 1:s4
        regr = [ones(size(freq.trialinfo,1),1) freq.powspctrm(:,i2,i3,i4)];
        beta = regr \ freq.trialinfo(:, cfg.powcorr.info);
        powspctrm(i2,i3,i4) = beta(2);
      end
    end
  end
  warning('on', 'MATLAB:rankDeficientMatrix')
  %---------------------------%
  
  %---------------------------%
  %-restructure freq for multiple subjects
  freq = rmfield(freq, 'trialinfo');
  freq.dimord = freq.dimord(5:end);
  freq.powspctrm = powspctrm;
  %---------------------------%
  
  save([cfg.dpow outputfile], 'freq')
  
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

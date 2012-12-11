function powcorr_subj(info, opt, subj)
%POWCORR_SUBJ power-trial correlation for each subject
% 
% INFO
%  .dpow: directory for POW data
%  .log: name of the file and directory to save log
%
% CFG.OPT
%  .source: read virtual electrode data (logical)
%  .cond*: cell with conditions (e.g. {'*cond1' '*cond2'})
%  .pow*: a structure with cfg to pass to ft_freqanalysis
%
%  .powcorr*: column of trialinfo to use for the correlation
%  .powlog: logical (take the log of power, strongly advised)
%
%  Baseline correction at the single-trial level:
%  .bl: if empty, no baseline. Otherwise:
%  .bl.baseline: two scalars with baseline windows
%  .bl.baselinetype: type of baseline ('relchange')
%
% IN
%  LOAD_DATA: data in /PROJ/subjects/SUBJ/MOD/NICK/
% OR if cfg.opt.source
%  LOAD_SOURCE: source in info.dsou after SOURCE_SUBJ
%
% OUT
%  [info.dpow 'powcorr_SUBJ_COND'] 'powcorr_s': power correlation for single-subject
%
% * indicates obligatory parameter
%
% Part of EVENTBASED single-subject
% see also ERP_SUBJ, ERP_GRAND, 
% ERPSOURCE_SUBJ, ERPSOURCE_GRAND, ERPSTAT_SUBJ, ERPSTAT_GRAND,
% POW_SUBJ, POW_GRAND, POW_GRP, POWCORR_SUBJ, POWCORR_GRAND,
% POWSOURCE_SUBJ, POWSOURCE_GRAND, POWSTAT_SUBJ, POWSTAT_GRAND,
% SOURCE_SUBJ, CONN_SUBJ, CONN_GRAND, CONN_STAT

%---------------------------%
%-start log
output = sprintf('%s (%04d) began at %s on %s\n', ...
  mfilename, subj, datestr(now, 'HH:MM:SS'), datestr(now, 'dd-mmm-yy'));
tic_t = tic;
%---------------------------%

if ~isfield(opt, 'powlog'); opt.powlog = false; end

%-------------------------------------%
%-loop over conditions
for k = 1:numel(opt.cond)
  cond = opt.cond{k};
  condname = regexprep(cond, '*', '');
  
  %---------------------------%
  %-read data
  if ~isfield(opt, 'source') || ~opt.source
    [data] = load_data(info, subj, cond);
  else
    [data] = load_source(info, subj, cond);
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
  cfg = opt.pow;
  cfg.feedback = 'none';
  cfg.keeptrials = 'yes';
  powcorr_s = ft_freqanalysis(cfg, data);
  
  if isfield(opt.pow, 'toi')
    powcorr_s.time = opt.pow.toi;
  end
  
  %-----------------%
  %-when no time info (mtmfft), create empty time
  if ~isfield(powcorr_s, 'time')
    powcorr_s.time = 0;
    powcorr_s.dimord = [powcorr_s.dimord '_time'];
  end
  %-----------------%
  %---------------------------%
  
  %---------------------------%
  %-baseline
  if isfield(opt, 'bl') && ~isempty(opt.bl)
    cfg = opt.bl;
    powcorr_s = ft_freqbaseline(cfg, powcorr_s);
  end
  %---------------------------%
  
  %---------------------------%
  %-regression at each point
  if opt.powlog
    powcorr_s.powspctrm = log(powcorr_s.powspctrm);
  end
  
  [s1 s2 s3 s4] = size(powcorr_s.powspctrm);
  
  powspctrm = nan(s2,s3,s4);
  warning('off', 'MATLAB:rankDeficientMatrix') % NaN in \
  
  for i2 = 1:s2
    for i3 = 1:s3
      for i4 = 1:s4
        regr = [ones(size(powcorr_s.trialinfo,1),1) powcorr_s.trialinfo(:, opt.powcorr)];
        beta = regr \ powcorr_s.powspctrm(:,i2,i3,i4);
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
  
  save([info.dpow outputfile], 'powcorr_s')
  
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

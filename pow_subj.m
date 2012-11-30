function pow_subj(info, opt, subj)
%POW_SUBJ power-analysis for each subject
%
% INFO
%  .dpow: directory for POW data
%  .log: name of the file and directory to save log
%
% CFG.OPT
%  .source: read virtual electrode data (logical)
%  .cond*: cell with conditions (e.g. {'*cond1' '*cond2'})
%  .pow*: a structure with cfg to pass to ft_freqanalysis
%  .planar: planar transformation, MEG-only (logical)
%
%  Baseline correction at the single-trial level:
%  .bl: if empty, no baseline. Otherwise:
%  .bl.baseline: two scalars with baseline windows
%  .bl.baselinetype: type of baseline ('relchange')
%  .bl.log: take the log BEFORE taking baseline (data becomes more normal)
%
% IN
%  LOAD_DATA: data in /PROJ/subjects/SUBJ/MOD/NICK/
% OR if cfg.opt.source
%  LOAD_SOURCE: source in info.dsou after SOURCE_SUBJ
%
% OUT
%  [info.dpow 'pow_SUBJ_COND'] 'pow_s': power analysis for single-subject
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

if ~isfield(opt, 'planar'); opt.planar = []; end
if ~isfield(opt, 'bl'); opt.bl = []; end

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
  
  outputfile = sprintf('pow_%04d_%s', subj, condname);
  %---------------------------%
  
  %---------------------------%
  %-calculate power
  %-----------------%
  %-planar
  if isfield(data, 'grad') && opt.planar
    
    cfg = [];
    cfg.grad = data.grad;
    cfg.method = 'distance';
    cfg.neighbourdist = info.sens.dist;
    nbor = ft_prepare_neighbours(cfg);
    
    cfg = [];
    cfg.neighbours = nbor;
    data = ft_megplanar(cfg, data);
    
  end
  %-----------------%
  
  %-----------------%
  %-power
  cfg = opt.pow;
  cfg.feedback = 'none';
  cfg.keeptrials = 'yes';
  pow_s = ft_freqanalysis(cfg, data);
  
  if isfield(opt.pow, 'toi')
    pow_s.time = opt.pow.toi;
  end
  %-----------------%
  
  %-----------------%
  %-when no time info (mtmfft), create empty time
  if ~isfield(pow_s, 'time')
    pow_s.time = 0;
    pow_s.dimord = [pow_s.dimord '_time'];
  end
  %-----------------%
  
  %-----------------%
  %-baseline
  if isfield(opt, 'bl') && ~isempty(opt.bl)

    if isfield(opt.bl, 'log')
      pow_s.powspctrm = log(pow_s.powspctrm);
      pow_s.powspctrm(isinf(pow_s.powspctrm)) = 0;
    end
    
    cfg = opt.bl;
    pow_s = ft_freqbaseline(cfg, pow_s);

  end
  %-----------------%
  
  %-----------------%
  %-average
  cfg = [];
  pow_s = ft_freqdescriptives(cfg, pow_s);
  %-----------------%
  
  %-----------------%
  %-planar
  if isfield(data, 'grad') && opt.planar
    cfg = [];
    pow_s = ft_combineplanar(cfg, pow_s);
  end
  %-----------------%
  %---------------------------%
  
  save([info.dpow outputfile], 'pow_s')
  
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

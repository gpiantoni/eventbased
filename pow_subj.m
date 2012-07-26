function pow_subj(cfg, subj)
%POW_SUBJ create subject-specific pow
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
%  .pow.cond: cell with conditions (e.g. {'*cond1' '*cond2'})'
%
%  .pow: a structure with cfg to pass to ft_freqanalysis
%
%  .pow.planar: planar transformation, MEG-only (logical)
%
% Baseline correction is applied in POW_GRAND
%
% IN:
%  data in /PROJ/subjects/SUBJ/MOD/NICK/
%
% OUT
%  [cfg.dpow 'pow_SUBJ_COND'] 'pow_s': power analysis for single-subject
%
% Part of EVENTBASED single-subject
% see also ERP_SUBJ, ERP_GRAND, 
% ERPSOURCE_SUBJ, ERPSOURCE_GRAND, ERPSTAT_SUBJ, ERPSTAT_GRAND,
% POW_SUBJ, POW_GRAND, POWCORR_SUBJ, POWCORR_GRAND,
% POWSOURCE_SUBJ, POWSOURCE_GRAND, POWSTAT_SUBJ, POWSTAT_GRAND,
% CONN_SUBJ, CONN_GRAND, CONN_STAT

%---------------------------%
%-start log
output = sprintf('%s (%04d) began at %s on %s\n', ...
  mfilename, subj, datestr(now, 'HH:MM:SS'), datestr(now, 'dd-mmm-yy'));
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
  %-------%
  %-planar
  if isfield(data, 'grad') && cfg.pow.planar
    
    tmpcfg = [];
    tmpcfg.grad = data.grad;
    tmpcfg.method = 'distance';
    tmpcfg.neighbourdist = cfg.sens.dist;
    nbor = ft_prepare_neighbours(tmpcfg);
    
    tmpcfg = [];
    tmpcfg.neighbours = nbor;
    data = ft_megplanar(tmpcfg, data);
    
  end
  %-------%
  
  cfg2 = cfg.pow;
  cfg2.feedback = 'etf';
  pow_s = ft_freqanalysis(cfg2, data);
  
  if isfield(cfg.pow, 'toi')
    pow_s.time = cfg.pow.toi;
  end
  
  %-------%
  %-when no time info (mtmfft), create empty time
  if ~isfield(pow_s, 'time')
    pow_s.time = 0;
    pow_s.dimord = [pow_s.dimord '_time'];
  end
  %-------%
  
  %-------%
  %-planar
  if isfield(data, 'grad') && cfg.pow.planar
    tmpcfg = [];
    pow_s = ft_combineplanar(tmpcfg, pow_s);
  end
  %-------%
  %---------------------------%
  
  save([cfg.dpow outputfile], 'pow_s')
  
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

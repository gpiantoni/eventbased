function erpsource_subj(info, opt, subj)
%ERPSOURCE_SUBJ: sources of ERP using LCMV for each subject
%
% INFO
%  .log: name of the file and directory to save log
%  .derp: directory for ERP data
%
% CFG.OPT
%  .cond*: cell with conditions (e.g. {'*cond1' '*cond2'})
%  .erp: a structure with cfg to pass to ft_timelockanalysis 
%  .keepfilter: keep filters or not, keep them only if you plan to use erpstat or conn analyses (logical)

%
%  .erpsource.refcomp: cell of string(s) of the comparison whose peaks
%                     will be localized (one of the cells of cfg.gerp.comp)
%

%  .erpsource.bline: one number in s, the center of the covariance window
%                    of the baseline (the window length depends on erp_peak) 
%                    If empty, no baseline.
%
%  .erpsource.lcmv: options that will be passed to beamformer. Examples:
%     .lambda: regularization parameter of beamformer ('25%')
%     .powmethod: power method of beamformer ('trace' or 'lambda1')
%     .refdip: location of the dipole for computing coherence to.
%

%
% IN:
%  data in /PROJ/subjects/SUBJ/MOD/NICK/
%
% OUT
%  [info.derp 'erpsource_SUBJ_COND'] 'erpsource_s_A': source data for period of interest for each subject
%  [info.derp 'erpsource_SUBJ_COND'] 'erpsource_s_B': source data for baseline for each subject
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

%---------------------------%
%-default cfg
if ~isfield(opt, 'lmcv'); opt.lmcv = []; end
if ~isfield(opt, 'keepfilter'); opt.keepfilter = false; end
%---------------------------%

%---------------------------%
%-dir and files
[vol, lead, sens] = load_headshape(info, subj);
%---------------------------%

erp_peak = get_peak(info, opt.peak, 'erp');

%-------------------------------------%
%-loop over conditions
for k = 1:numel(cfg.erpsource.cond)
  cond     = cfg.erpsource.cond{k};
  condname = regexprep(cond, '*', '');
  
  %---------------------------%
  %-read data
  [data badchan] = load_data(cfg, subj, cond);
  if isempty(data)
    output = sprintf('%sCould not find any file for condition %s\n', ...
      output, cond);
    continue
  end

  clear erpsource_s_A erpsource_s_B
  outputfile = sprintf('erpsource_%04d_%s', subj, condname);
  %---------------------------%
  
  %---------------------------%
  %-remove bad channels from leadfield
  datachan = ft_channelselection([{'all'}; cellfun(@(x) ['-' x], badchan, 'uni', false)], data.label);
  [leadchan] = prepare_leadchan(lead, datachan);
  %---------------------------%
  
  for p = 1:numel(erp_peak)
    
    fprintf('\n   ->->-> Running % 2d erp_peak (%s) <-<-<-\n', p, erp_peak(p).name);
    
    %---------------------------%
    %-parameters
    %-----------------%
    %-covariance window
    cfgerp = cfg.erpsource.erp;
    cfgerp.covariance = 'yes';
    cfgerp.feedback = 'none';
    cfgerp.channel = datachan;
    %-----------------%
    
    %-----------------%
    %-source analysis
    cfgsou = [];
    
    cfgsou.method = 'lcmv';
    cfgsou.lcmv = cfg.erpsource.lcmv;
    cfgsou.lcmv.feedback = 'none';
    
    cfgsou.vol = vol;
    cfgsou.grid = leadchan;
    cfgsou.elec = sens;
    cfgsou.feedback = 'none';
    cfgsou.lcmv.keepmom = 'no';
    if cfg.erpsource.keepfilter
      cfgsou.lcmv.keepfilter   = 'yes';
      cfgsou.lcmv.realfilter   = 'yes';
    end
    %-----------------%
    %---------------------------%
    
    %---------------------------%
    %-baseline
    if isfield(cfg.erpsource, 'bline') && ~isempty(cfg.erpsource.bline)
      
      %-----------------%
      %-covariance window
      cfgerp.covariancewindow = cfg.erpsource.bline  + erp_peak(p).wndw * [-.5 .5];
      avgPre = ft_timelockanalysis(cfgerp, data);
      %-----------------%
      
      %-----------------%
      %-source analysis
      source = ft_sourceanalysis(cfgsou, avgPre);
      source.cfg = [];
      %-----------------%
      
      %-----------------%
      %-realign source
      erpsource_s_B(p,:) = realign_source(cfg, subj, source);
      %-----------------%
      
    else
      erpsource_s_B = [];
      
    end
    %---------------------------%
    
    %---------------------------%
    %-main effect
    %-----------------%
    %-covariance
    cfgerp.covariancewindow = erp_peak(p).time + erp_peak(p).wndw * [-.5 .5];
    avgPost = ft_timelockanalysis(cfgerp, data);
    %-----------------%
    
    %-----------------%
    %-source
    source = ft_sourceanalysis(cfgsou, avgPost);
    
    chan = source.cfg.channel;
    source.cfg = [];
    source.cfg.channel = chan;
    %-----------------%
    
    %-----------------%
    %-realign source
    erpsource_s_A(p,:) = realign_source(cfg, subj, source);
    %-----------------%
    %---------------------------%
    
  end
  
  %-----------------%
  %-save source
  save([info.derp outputfile], 'erpsource_s_A', 'erpsource_s_B', '-v7.3')
  %-----------------%
  
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

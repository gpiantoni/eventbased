function erpstat_subj(info, opt, subj)
%ERPSTAT_SUBJ: use LCMV common filters for each condition
%
% CFG
%  .data: path of /data1/projects/PROJ/subjects/
%  .rec: REC in /data1/projects/PROJ/recordings/REC/
%  .nick: NICK in /data1/projects/PROJ/subjects/0001/MOD/NICK/
%  .mod: modality, MOD in /data1/projects/PROJ/subjects/0001/MOD/NICK/
%  .endname: includes preprocessing steps (e.g. '_seldata_gclean_redef')
%
%  .log: name of the file and directory to save log
%  .derp: directory for ERP data
%
%  .vol.type: 'template' or subject-specific ('dipoli' or 'openmeeg' or 'bemcp')
%  (if cfg.vol.type == 'template')
%      .vol.template: file with template containing vol, lead, sens
%
%  .sourcespace: 'surface' 'volume' 'volume_warp'
%  (if cfg.sourcespace == 'surface')
%  .SUBJECTS_DIR: where the Freesurfer data is stored, like the environmental variable
%
%  .erpstat.refcond: string of the condition used for LCMV-filter
%  .erpstat.cond: cell with conditions to calculate (e.g. {'*cond1' '*cond2'})
%
%  .erpsource.areas: how to speficy peaks to analyze, 'manual' or 'erp_peak' (peaks from granderp)
%  (if .erpsource.areas == 'manual')
%  .erpsource.erp_peak: struct with multiple elements 
%              .name: string ('name_of_the_time_window')
%              .time: scalar (center of the time window in s)
%              .wndw: scalar (length of the time window in s)
%  (if .erpsource.areas == 'erp_peak')
%  .erpsource.refcomp: cell of string(s) of the comparison whose peaks
%                     will be localized (one of the cells of cfg.gerp.comp)
%
%  .erpsource.erp: a structure with cfg to pass to ft_timelockanalysis
%  .erpsource.bline: one number in s, the center of the covariance window
%                    of the baseline (the window length depends on erp_peak) 
%                    If empty, no baseline.
%
%  .erpsource.lcmv: options that will be passed to beamformer. Examples:
%     .lambda: regularization parameter of beamformer ('25%')
%     .powmethod: power method of beamformer ('trace' or 'lambda1')
%
% IN:
%  data in /PROJ/subjects/SUBJ/MOD/NICK/
%  [info.derp 'erpsource_SUBJ_COND'] 'erpsource_s_A': source data for period of interest for each subject
%  [info.derp 'erpsource_SUBJ_COND'] 'erpsource_s_B': source data for baseline for each subject
%
% OUT
%  [info.derp 'erpstat_SUBJ_COND'] 'erpstat_s_A': source data for period of interest for each subject, after common filters
%  [info.derp 'erpstat_SUBJ_COND'] 'erpstat_s_B': source data for baseline for each subject, after common filters
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

%---------------------------%
%-dir and files
[vol, lead, sens] = load_headshape(info, opt, subj);

%-----------------%
%-load source
souname = regexprep(cfg.erpstat.refcond, '*', '');
sourcefile = sprintf('erpsource_%04d_%s', subj, souname);
load([info.derp sourcefile], 'erpsource_s_A', 'erpsource_s_B')
souchan = sourceA{1}.cfg.channel;
%-----------------%

erp_peak = getpeak(cfg, 'erp');
%---------------------------%

%-------------------------------------%
%-loop over conditions
for k = 1:numel(cfg.erpstat.cond)
  cond     = cfg.erpstat.cond{k};
  condname = regexprep(cond, '*', '');
  
  %---------------------------%
  %-read data
  [data badchan] = load_data(cfg, subj, cond);
  if isempty(data)
    output = sprintf('%sCould not find any file for test %s\n', ...
      output, cond);
    continue
  end
  
  outputfile = sprintf('erpstat_%04d_%s', subj, condname);
  %---------------------------%
  
  %---------------------------%
  %-remove bad channels from leadfield
  datachan = ft_channelselection([{'all'}; cellfun(@(x) ['-' x], badchan, 'uni', false)], data.label);
  
  %-------%
  %-check if the channels match
  % It can happen that they don't match if they come from slighly different
  % datasets. However, it's no problem if setdiff(datachan, souchan) is not
  % empty, which means that they are good channels in the data that were
  % not used in the filter.
  % It's a "problem" if setdiff(souchan, datachan) is not empty, because it
  % means that channels that were used in the computation of the filter are
  % interpolated channels in the data here. However, it's just
  % interpolation, not really a huge problem.
  chandiff = setdiff(souchan, datachan);
  if ~isempty(chandiff)
    output = sprintf('%sUsing interpolated channels for source-reconstruction: %s\n', ...
      output, sprintf('%s, ', chandiff{:}));
  end
  datachan = souchan;
  %-------%
  
  [leadchan] = prepare_leadchan(lead, datachan);
  %---------------------------%
  
  for p = 1:numel(erp_peak)
    
    fprintf('\n   ->->-> Running % 2d erpstat (%s) <-<-<-\n', p, erp_peak(p).name);
    
    %---------------------------%
    %-baseline
    %-----------------%
    %-covariance window
    cfgerp = cfg.erpsource.erp;
    cfgerp.covariance = 'yes';
    cfgerp.covariancewindow = cfg.erpsource.bline  + erp_peak(p).wndw * [-.5 .5];
    cfgerp.feedback = 'none';
    cfgerp.channel = datachan;
    
    avgPre = ft_timelockanalysis(cfgerp, data);
    %-----------------%
    
    %-----------------%
    %-source analysis
    cfgsou              = [];
    
    cfgsou.method   = 'lcmv';
    cfgsou.lcmv = cfg.erpsource.lcmv;
    cfgsou.lcmv.feedback = 'none';
    
    cfgsou.vol = vol;
    cfgsou.grid = leadchan;
    cfgsou.elec = sens;
    cfgsou.feedback = 'none';
    cfgsou.lcmv.keepmom = 'no';
    
    cfgsou.grid.filter  = erpsource_s_B{p}.avg.filter;
    
    source = ft_sourceanalysis(cfgsou, avgPre);
    source.cfg = [];
    %-----------------%
    
    %-----------------%
    %-realign source
    erpstat_s_B(p,:) = realign_source(cfg, subj, source);
    %-----------------%
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
    cfgsou.grid.filter  = erpsource_s_A{p}.avg.filter;
    source = ft_sourceanalysis(cfgsou, avgPost);
    
    chan = source.cfg.channel;
    source.cfg = [];
    source.cfg.channel = chan;
    %-----------------%
    
    %-----------------%
    %-realign source
    erpstat_s_A(p,:) = realign_source(cfg, subj, source);
    %-----------------%
    %---------------------------%

  end
  
  %-----------------%
  %-save source
  save([info.derp outputfile], 'erpstat_s_A', 'erpstat_s_B', '-v7.3')
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

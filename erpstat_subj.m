function erpstat_subj(cfg, subj)
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
%  .vol.type: 'template' or subject-specific ('dipoli' or 'openmeeg')
%    if 'template'
%      .vol.template: file with template containing vol, lead, sens
%    if ~ 'template'
%      .bnd2lead.mni.warp: logical (optional. Instead of transforming the
%      brain into MNI coordinates, you can wrap the grid onto it)
%
%  .erpstat.refcond: string of the condition used for LCMV-filter
%  .erpstat.cond: cell with conditions to calculate (e.g. {'*cond1' '*cond2'})
%
%  .erpsource.peaks: how to speficy peaks to analyze, 'manual' or 'erp_peak' (peaks from granderp)
%    if 'manual'
%      .erpsource.erp_peak(1).name: string ('name_of_the_time_window')
%      .erpsource.erp_peak(1).time: scalar (center of the time window in s)
%      .erpsource.erp_peak(1).wndw: scalar (length of the time window in s)
%    if 'erp_peak'
%      .erpsource.refcomp: string of the condition whose peaks will be localized
%
%  .erpsource.bline: one number in s, the center of the covariance window of the baseline (the window length depends on erp_peak)
%
%  .erpsource.lcmv: options that will be passed to beamformer. Examples:
%     .lambda: regularization parameter of beamformer ('25%')
%     .powmethod: power method of beamformer ('trace' or 'lambda1')
%
% IN:
%  data in /PROJ/subjects/SUBJ/MOD/NICK/
%  [cfg.derp 'erpsource_SUBJ_COND'] 'erpsource_s_A': source data for period of interest for each subject
%  [cfg.derp 'erpsource_SUBJ_COND'] 'erpsource_s_B': source data for baseline for each subject
%
% OUT
%  [cfg.derp 'erpstat_SUBJ_COND'] 'erpstat_s_A': source data for period of interest for each subject, after common filters
%  [cfg.derp 'erpstat_SUBJ_COND'] 'erpstat_s_B': source data for baseline for each subject, after common filters
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
%-dir and files
[vol, lead, sens] = load_headshape(cfg, subj);

%-----------------%
%-load source
souname = regexprep(cfg.erpstat.refcond, '*', '');
sourcefile = sprintf('erpsource_%04d_%s', subj, souname);
load([cfg.derp sourcefile], 'erpsource_s_A', 'erpsource_s_B')
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
    cfg2 = cfg.erpsource.erp;
    cfg2.covariance = 'yes';
    cfg2.covariancewindow = cfg.erpsource.bline  + erp_peak(p).wndw * [-.5 .5];
    cfg2.feedback = 'none';
    cfg2.channel = datachan;
    
    avgPre = ft_timelockanalysis(cfg2, data);
    %-----------------%
    
    %-----------------%
    %-source analysis
    cfg3              = [];
    
    cfg3.method   = 'lcmv';
    cfg3.lcmv = cfg.erpsource.lcmv;
    cfg3.lcmv.feedback = 'none';
    
    cfg3.vol = vol;
    cfg3.grid = leadchan;
    cfg3.elec = sens;
    cfg3.feedback = 'none';
    cfg3.lcmv.keepmom = 'no';
    
    cfg3.grid.filter  = erpsource_s_B{p}.avg.filter;
    
    erpstat_s_B{p} = ft_sourceanalysis(cfg3, avgPre);
    erpstat_s_B{p}.cfg = [];
    %-----------------%
    
    %-----------------%
    %-load MNI grid
    if ~strcmp(cfg.vol.type, 'template') ...
        && isfield(cfg, 'bnd2lead') && isfield(cfg.bnd2lead, 'mni') ...
        && isfield(cfg.bnd2lead.mni, 'warp') && cfg.bnd2lead.mni.warp
      
      load(sprintf('%s/template/sourcemodel/standard_grid3d%dmm.mat', ...
        fileparts(which('ft_defaults')), cfg.bnd2lead.mni.resolution), 'grid');
      grid = ft_convert_units(grid, 'mm');
      
      erpstat_s_B{p}.pos = grid.pos;
    end
    %-----------------%
    %---------------------------%
    
    %---------------------------%
    %-main effect
    %-----------------%
    %-covariance
    cfg2.covariancewindow = erp_peak(p).time + erp_peak(p).wndw * [-.5 .5];
    avgPost = ft_timelockanalysis(cfg2, data);
    %-----------------%
    
    %-----------------%
    %-source
    cfg3.grid.filter  = erpsource_s_A{p}.avg.filter;
    erpstat_s_A{p} = ft_sourceanalysis(cfg3, avgPost);
    
    chan = erpstat_s_A{p}.cfg.channel;
    erpstat_s_A{p}.cfg = [];
    erpstat_s_A{p}.cfg.channel = chan;
    %-----------------%
    
    %-----------------%
    %-load MNI grid
    if ~strcmp(cfg.vol.type, 'template') ...
        && isfield(cfg, 'bnd2lead') && isfield(cfg.bnd2lead, 'mni') ...
        && isfield(cfg.bnd2lead.mni, 'warp') && cfg.bnd2lead.mni.warp
      
      erpstat_s_A{p}.pos = grid.pos;
    end
    %-----------------%
    %---------------------------%

  end
  
  %-----------------%
  %-save source
  save([cfg.derp outputfile], 'erpstat_s_A', 'erpstat_s_B', '-v7.3')
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
fid = fopen([cfg.log '.txt'], 'a');
fwrite(fid, output);
fclose(fid);
%-----------------%
%---------------------------%

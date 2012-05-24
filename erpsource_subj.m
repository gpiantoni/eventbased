function erpsource_subj(cfg, subj)
%ERPSOURCE_SUBJ: identify sources from erp peaks using LCMV
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
%  .erpsource.cond: cell with conditions (e.g. {'*cond1' '*cond2'})
%
%  .vol.type: 'template' or subject-specific ('dipoli' or 'openmeeg')
%    if 'template'
%      .vol.template: file with template containing vol, lead, sens
%    if ~ 'template'
%      .bnd2lead.mni.warp: logical (optional. Instead of transforming the
%      brain into MNI coordinates, you can wrap the grid onto it)
%
%  .erpsource.areas: how to speficy peaks to analyze, 'manual' or 'erp_peak' (peaks from granderp)
%    if 'manual'
%      .erpsource.erp_peak(1).name: string ('name_of_the_time_window')
%      .erpsource.erp_peak(1).time: scalar (center of the time window in s)
%      .erpsource.erp_peak(1).wndw: scalar (length of the time window in s)
%    if 'erp_peak'
%      .erp.refcond: string of the comparison whose peaks will be localized
%
%  .erpsource.erp: a structure with cfg to pass to ft_timelockanalysis
%  .erpsource.bline: one number in s, the center of the covariance window of the baseline (the window length depends on erp_peak)
%
%  .erpsource.lcmv: options that will be passed to beamformer. Examples:
%     .lambda: regularization parameter of beamformer ('25%')
%     .powmethod: power method of beamformer ('trace' or 'lambda1')
%     .refdip: location of the dipole for computing coherence to.
%
%  .erpsource.keepfilter: logical, to keep filters or not (keep them only
%                          if you plan to use erpstat or conn analyses)
%
% IN:
%  data in /PROJ/subjects/SUBJ/MOD/NICK/
%
% OUT
%  [cfg.derp 'erpsource_SUBJ_COND'] 'erpsource_subj_A': source data for period of interest for each subject
%  [cfg.derp 'erpsource_SUBJ_COND'] 'erpsource_subj_B': source data for baseline for each subject
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

%---------------------------%
%-dir and files
[vol, lead, sens] = load_headshape(cfg, subj);
%---------------------------%

%-------------------------------------%
%-loop over conditions
for k = 1:numel(cfg.erpsource.cond)
  cond     = cfg.erpsource.cond{k};
  condname = regexprep(cond, '*', '');
  
  %---------------------------%
  %-use predefined or erp_peaks for areas of interest
  if strcmp(cfg.erpsource.areas, 'manual')
    erp_peak = cfg.erpsource.erp_peak;
  elseif strcmp(cfg.erpsource.areas, 'erp_peak')
    peakname = regexprep(cfg.erp.refcond, '*', '');
    load([cfg.derp 'erp_peak_' peakname], 'erp_peak')
  end
  %---------------------------%
  
  %---------------------------%
  %-read data
  [data badchan] = load_data(cfg, subj, cond);
  if isempty(data)
    output = sprintf('%sCould not find any file for condition %s\n', ...
      output, cond);
    continue
  end

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
    if cfg.erpsource.keepfilter
      cfg3.lcmv.keepfilter   = 'yes';
      cfg3.lcmv.realfilter   = 'yes';
    end
    
    erpsource_subj_B{p} = ft_sourceanalysis(cfg3, avgPre);
    erpsource_subj_B{p}.cfg = [];
    %-----------------%
    
    %-----------------%
    %-load MNI grid
    if ~strcmp(cfg.vol.type, 'template') ...
        && isfield(cfg, 'bnd2lead') && isfield(cfg.bnd2lead, 'mni') ...
        && isfield(cfg.bnd2lead.mni, 'warp') && cfg.bnd2lead.mni.warp
      
      load(sprintf('%s/template/sourcemodel/standard_grid3d%dmm.mat', ...
        fileparts(which('ft_defaults')), cfg.bnd2lead.mni.resolution), 'grid');
      grid = ft_convert_units(grid, 'mm');
      
      erpsource_subj_B{p}.pos = grid.pos;
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
    erpsource_subj_A{p} = ft_sourceanalysis(cfg3, avgPost);

    chan = erpsource_subj_A{p}.cfg.channel;
    erpsource_subj_A{p}.cfg = [];
    erpsource_subj_A{p}.cfg.channel = chan;
    %-----------------%
    
    %-----------------%
    %-load MNI grid
    if ~strcmp(cfg.vol.type, 'template') ...
        && isfield(cfg, 'bnd2lead') && isfield(cfg.bnd2lead, 'mni') ...
        && isfield(cfg.bnd2lead.mni, 'warp') && cfg.bnd2lead.mni.warp
      
      erpsource_subj_A{p}.pos = grid.pos;
    end
    %-----------------%
    %---------------------------%
    
  end
  
  %-----------------%
  %-save source
  save([cfg.derp outputfile], 'erpsource_subj_A', 'erpsource_subj_B', '-v7.3')
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
function erpsource_subj(cfg, subj)
%ERPSOURCE_SUBJ: identify sources from erp peaks using LCMV
%
% CFG
%  .data: name of projects/PROJNAME/subjects/
%  .mod: name of the modality used in recordings and projects
%  .cond: name to be used in projects/PROJNAME/subjects/0001/MOD/CONDNAME/
%  .endname: includes previous steps '_seldata_gclean_preproc_redef'
%  .log: name of the file and directory with analysis log
%  .test: a cell with the condition defined by redef. This function will loop over cfg.test
%  .derp: directory to save ERP data
%
%  .erpeffect: effect of interest for source reconstruction, can be a vector (this field is shared with erp_grand.m, maybe it's not a good idea)
%  .erp: a structure with cfg to pass to ft_timelockanalysis
%
%  .vol.type: 'template' or subject-specific ('dipoli' or 'openmeeg')
%    if template, specify template .vol.template (should contain vol, lead, sens)
%    if not template, specify
%      .vol.mod: 'smri'
%      .vol.cond: 't1'
%      .proj: because the project name is part of the MRI name
%
%  .erpsource.areas: how to speficy peaks to analyze, 'manual' or 'erppeak' (peaks from granderp)
%    if 'manual'
%      .erpsource.erppeak(1).name = 'name_of_the_time_ window';
%      .erpsource.erppeak(1).time = 0.10; % center of the time window
%      .erpsource.erppeak(1).wndw = 0.05; % length of the time window
%    if 'erppeak', it reads the significant peaks calculated by erp_grand
%                  erppeak is specific to each condition
%
%  .erpsource.erp: a structure with cfg to pass to ft_timelockanalysis (it's better if identical to cfg.erp)
%  .erpsource.bline: one number in s, the center of the covariance window of the baseline (the window length depends on erppeak)
%
%  .erpsource.lambda: regularization parameter of beamformer ('10%')
%  .erpsource.powmethod: power method of beamformer ('trace' or 'lambda1')
%
% OUT
%  [cfg.derp 'erpsource_001_TEST']: source data for period of interest and baseline for each subject
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

%---------------------------%
%-dir and files
[vol, lead, sens] = load_headshape(cfg, subj);
%---------------------------%

%-------------------------------------%
%-loop over conditions
for k = 1:numel(cfg.erpsource.cond) % DOC: CFG.ERPSOURCE.COND
  cond     = cfg.erpsource.cond{k};
  condname = regexprep(cond, '*', '');
  
  %---------------------------%
  %-use predefined or erppeaks for areas of interest
  if strcmp(cfg.erpsource.areas, 'manual')
    erppeak = cfg.erpsource.erppeak;
  elseif strcmp(cfg.erpsource.areas, 'erppeak')
    peakname = regexprep(cfg.erp.refcond, '*', ''); % DOC: CFG.ERP.REFCOND
    load([cfg.derp cfg.cond peakname '_erppeak'], 'erppeak')
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

  outputfile = sprintf('erpsource_%02.f_%s', subj, condname);
  %---------------------------%
  
  %---------------------------%
  %-remove bad channels from leadfield
  datachan = ft_channelselection([{'all'}; cellfun(@(x) ['-' x], badchan, 'uni', false)], data.label);
  [leadchan] = prepare_leadchan(lead, datachan);
  %---------------------------%
  
  for f = 1:numel(erppeak)
    
    %---------------------------%
    %-baseline
    %-----------------%
    %-covariance window
    cfg2 = cfg.erpsource.erp;
    cfg2.covariance = 'yes';
    cfg2.covariancewindow = cfg.erpsource.bline  + erppeak(f).wndw * [-.5 .5];
    cfg2.feedback = 'none';
    cfg2.channel = datachan;
    
    avgPre = ft_timelockanalysis(cfg2, data);
    %-----------------%
    
    %-----------------%
    %-source analysis
    cfg3              = [];
    
    cfg3.method   = 'lcmv';
    cfg3.lcmv.lambda = cfg.erpsource.lambda;
    cfg3.lcmv.powmethod = cfg.erpsource.powmethod;
    cfg3.lcmv.feedback = 'none';
    
    cfg3.vol      = vol;
    cfg3.grid     = leadchan;
    cfg3.elec     = sens;
    cfg3.feedback = 'none';
    
    souPre{f} = ft_sourceanalysis(cfg3, avgPre);
    souPre{f}.cfg = [];
    %-----------------%
    
    %-----------------%
    %-load MNI grid
    if ~strcmp(cfg.vol.type, 'template') ...
        && isfield(cfg, 'bnd2lead') && isfield(cfg.bnd2lead, 'mni') ...
        && isfield(cfg.bnd2lead.mni, 'warp') && cfg.bnd2lead.mni.warp
      
      load(sprintf('%s/template/sourcemodel/standard_grid3d%dmm.mat', ...
        fileparts(which('ft_defaults')), cfg.bnd2lead.mni.resolution), 'grid');
      grid = ft_convert_units(grid, 'mm');
      
      souPre{f}.pos = grid.pos;
    end
    %-----------------%
    %---------------------------%
    
    %---------------------------%
    %-main effect
    %-----------------%
    %-covariance
    cfg2.covariancewindow = erppeak(f).time + erppeak(f).wndw * [-.5 .5];
    avgPost = ft_timelockanalysis(cfg2, data);
    %-----------------%
    
    %-----------------%
    %-source
    cfg3.lcmv.keepfilter = 'yes'; % it's used by gosdconn
    cfg3.lcmv.realfilter = 'yes';
    cfg3.lcmv.keepmom    = 'no';
    source{f} = ft_sourceanalysis(cfg3, avgPost);
    
    source{f}.avg.nai = log(source{f}.avg.pow ./ souPre{f}.avg.pow);
    chan = source{f}.cfg.channel;
    source{f}.cfg = [];
    source{f}.cfg.channel = chan;
    %-----------------%
    
    %-----------------%
    %-load MNI grid
    if ~strcmp(cfg.vol.type, 'template') ...
        && isfield(cfg, 'bnd2lead') && isfield(cfg.bnd2lead, 'mni') ...
        && isfield(cfg.bnd2lead.mni, 'warp') && cfg.bnd2lead.mni.warp
      
      source{f}.pos = grid.pos;
    end
    %-----------------%
    %---------------------------%
    
  end
  
  %-----------------%
  %-save source
  save([cfg.derp outputfile], 'source', 'souPre', '-v7.3')
  %-----------------%
  
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
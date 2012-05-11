function powsource_subj(cfg, subj)
%POWSOURCE_SUBJ: identify sources from pow peaks using DICS
%
% CFG
%  .data: path of /data1/projects/PROJNAME/subjects/
%  .rec: RECNAME in /data1/projects/PROJNAME/recordings/RECNAME/
%  .nick: NICKNAME in /data1/projects/PROJNAME/subjects/0001/MOD/NICKNAME/
%  .mod: modality, MOD in /data1/projects/PROJNAME/subjects/0001/MOD/NICKNAME/
%  .endname: includes preprocessing steps (e.g. '_seldata_gclean_preproc_redef')
%
%  .log: name of the file and directory to save log
%  .dpow: directory for POW data
%  .powsource.cond: cell with conditions (e.g. {'*cond1' '*cond2'})'
%
%  .vol.type: 'template' or subject-specific ('dipoli' or 'openmeeg')
%    if 'template'
%      .vol.template: file with template containing vol, lead, sens
%    if ~ 'template'
%      .bnd2lead.mni.warp: logical (optional. Instead of transforming the
%      brain into MNI coordinates, you can wrap the grid onto it)
%
%  .powsource.areas: how to speficy peaks to analyze, 'manual' or 'erppeak' (peaks from granderp)
%    if 'manual'
%      .powsource.erppeak(1).name: string ('name_of_the_time_window')
%      .powsource.erppeak(1).time: scalar (center of the time window in s)
%      .powsource.erppeak(1).wndw: scalar (length of the time window in s)
%      .powsource.powpeak(1).freq = 10; % center of the frequency
%      .powsource.powpeak(1).band = 4; % width of the frequency band
%    if 'powpeak'
%      .pow.refcond: string of the condition whose peaks will be localized
%    if 'powcorrpeak'
%      .powcorr.refcond: string of the condition whose peaks will be localized
%
%  .powsource.bline: one number in s, the center of the covariance window of the baseline (the window length depends on powpeak)
%
%  .powsource.lambda: regularization parameter of beamformer ('10%')
%  .powsource.powmethod: power method of beamformer ('trace' or 'lambda1')
%
% IN:
%  data in /PROJNAME/subjects/SUBJCODE/MOD/NICKNAME/
%
% OUT
%  [cfg.dpow 'powsource_SUBJCODE_COND']: source data for period of interest and baseline for each subject
%
% Part of EVENTBASED single-subject
% see also ERP_SUBJ, ERP_GRAND, ERPSOURCE_SUBJ, ERPSOURCE_GRAND,
% POW_SUBJ, POW_GRAND, POWSOURCE_SUBJ, POWSOURCE_GRAND,
% POWCORR_SUBJ, POWCORR_GRAND, POWSTAT_SUBJ, POWSTAT_GRAND, 
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
for k = 1:numel(cfg.powsource.cond) % DOC: CFG.POWSOURCE.COND
  cond     = cfg.powsource.cond{k};
  condname = regexprep(cond, '*', '');

  %---------------------------%
  %-use predefined or power-peaks for areas of interest
  if strcmp(cfg.powsource.areas, 'manual')
    powpeak = cfg.powsource.powpeak;
    
  elseif strcmp(cfg.powsource.areas, 'powpeak')
    peakname = regexprep(cfg.pow.refcond, '*', ''); % DOC: CFG.POW.REFCOND
    load([cfg.dpow cfg.cond peakname '_powpeak'], 'powpeak')
    
  elseif strcmp(cfg.powsource.areas, 'powcorrpeak')
    peakname = regexprep(cfg.powcorr.refcond, '*', ''); % DOC: CFG.POWCORR.REFCOND
    load([cfg.dpow cfg.cond peakname '_powcorrpeak'], 'powcorrpeak')
    powpeak = powcorrpeak;
    
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
  
  outputfile = sprintf('powsource_%02.f_%s', subj, condname);
  %---------------------------%

  %---------------------------%
  %-remove bad channels from leadfield
  datachan = ft_channelselection([{'all'}; cellfun(@(x) ['-' x], badchan, 'uni', false)], data.label);
  [leadchan] = prepare_leadchan(lead, datachan);
  %---------------------------%
  
  for f = 1:numel(powpeak)
    
    fprintf('\n   ->->-> Running % 2.f powpeak (%s) <-<-<-\n', f, powpeak(f).name);
    
    %---------------------------%
    %-more precise freq analysis reconstruction
    freqparam = prepare_freqpeak(cfg, powpeak(f), data.time{1}(1));
    %---------------------------%
    
    %---------------------------%
    %-freq analysis
    cfg1 = [];
    cfg1.method = 'mtmconvol';
    cfg1.output = 'fourier';
    
    cfg1.toi = [cfg.powsource.bline freqparam.time];
    cfg1.t_ftimwin = freqparam.wndw;

    cfg1.foi = freqparam.freq;
    if freqparam.dpss
      cfg1.taper = 'dpss';
      cfg1.tapsmofrq = freqparam.band;
    else
      cfg1.taper = 'hanning';
    end
    
    cfg1.feedback = 'none';
    cfg1.channel = datachan;
    freq = ft_freqanalysis(cfg1, data);
    %---------------------------%
    
    %---------------------------%
    %-baseline
    %-----------------%
    %-source analysis
    cfg1 = [];
    cfg1.latency = cfg.powsource.bline;
    cfg1.frequency = freqparam.freq;
    
    cfg1.method = 'dics';
    cfg1.dics.feedback = 'none';
    cfg1.dics.lambda = cfg.powsource.lambda;
    cfg1.dics.powmethod = cfg.powsource.powmethod;
    
    cfg1.vol          = vol;
    cfg1.grid         = leadchan;
    cfg1.elec         = sens;
    
    souPre{f}       = ft_sourceanalysis(cfg1, freq);
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
      source{f}.pos = grid.pos;
    end
    %-----------------%
    %---------------------------%

    %---------------------------%
    %-main analysis
    %-----------------%
    %-keep filters only for source analysis
    cfg1.(cfg1.method).keepfilter   = 'yes';
    cfg1.(cfg1.method).realfilter   = 'yes';
    
    cfg1.latency    = freqparam.time;
    source{f}       = ft_sourceanalysis(cfg1, freq);
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
  save([cfg.dpow outputfile], 'source', 'souPre', '-v7.3')
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

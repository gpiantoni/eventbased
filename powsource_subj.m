function powsource_subj(cfg, subj)
%POWSOURCE_SUBJ: identify sources from pow peaks using DICS
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
%  .poweffect: index of interest to create powpeak, can be a row vector (this field is shared with pow_grand.m, maybe it's not a good idea)
%
%  .vol.type: 'template' or subject-specific ('dipoli' or 'openmeeg')
%    if template, specify template .vol.template (should contain vol, lead, sens)
%    if not template, specify 
%      .vol.mod: 'smri'
%      .vol.cond: 't1'
%      .proj: because the project name is part of the MRI name
%
%  .powsource.areas: how to speficy peaks to analyze, 'manual' or 'powpeak' (peaks from grandpow)
%    if 'manual'
%      .powsource.powpeak(1).name = 'name_of_the_time_ window';
%      .powsource.powpeak(1).time = 0.10; % center of the time window
%      .powsource.powpeak(1).wndw = 0.05; % length of the time window
%      .powsource.powpeak(1).freq = 10; % center of the frequency
%      .powsource.powpeak(1).band = 4; % width of the frequency band
%    if 'powpeak', it reads the significant peaks calculated by pow_grand
%                  powpeak is specific to each condition
%
%  .powsource.bline: one number in s, the center of the covariance window of the baseline (the window length depends on powpeak)
%
%  .powsource.lambda: regularization parameter of beamformer ('10%')
%  .powsource.powmethod: power method of beamformer ('trace' or 'lambda1')
%
% OUT
%  [cfg.dpow 'powsource_001_TEST']: source data for period of interest and baseline for each subject
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
for p = cfg.poweffect

  %---------------------------%
  %-use predefined or power-peaks for areas of interest
  if strcmp(cfg.powsource.areas, 'manual')
    powpeak = cfg.powsource.powpeak;
    
  elseif strcmp(cfg.powsource.areas, 'powpeak')
    condname = regexprep(cfg.test{p}, '*', '');
    load([cfg.dpow cfg.cond condname '_powpeak'], 'powpeak')
    
  elseif strcmp(cfg.powsource.areas, 'powcorrpeak')
    load([cfg.dpow cfg.cond '_powcorrpeak'], 'powcorrpeak')
    powpeak = powcorrpeak;
    
  end
  %---------------------------%
  
  %---------------------------%
  %-read data
  [data badchan] = load_data(cfg, subj, p);
  if isempty(data)
    output = sprintf('%sCould not find any file for test %s\n', ...
      output, cfg.test{p});
    continue
  end
  
  condname = regexprep(cfg.test{p}, '*', '');
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

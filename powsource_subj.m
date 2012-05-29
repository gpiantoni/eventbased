function powsource_subj(cfg, subj)
%POWSOURCE_SUBJ: identify sources from pow peaks using DICS
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
%  .powsource.cond: cell with conditions (e.g. {'*cond1' '*cond2'})'
%
%  .vol.type: 'template' or subject-specific ('dipoli' or 'openmeeg')
%    if 'template'
%      .vol.template: file with template containing vol, lead, sens
%    if ~ 'template'
%      .bnd2lead.mni.warp: logical (optional. Instead of transforming the
%      brain into MNI coordinates, you can wrap the grid onto it)
%
%  .powsource.peaks: how to speficy peaks to analyze, 'manual' or 'pow_peak'
%          (peaks from grandpow) or 'powcorr_peak' (peaks from grandpowcorr)
%    if 'manual'
%      .powsource.pow_peak(1).name: string ('name_of_the_time_window')
%      .powsource.pow_peak(1).time: scalar (center of the time window in s)
%      .powsource.pow_peak(1).wndw: scalar (length of the time window in s)
%      .powsource.pow_peak(1).freq = 10; % center of the frequency
%      .powsource.pow_peak(1).band = 4; % width of the frequency band
%    if 'pow_peak' or 'powcorr_peak'
%      .powsource.refcomp: cell of string(s) of the comparison whose peaks
%                     will be localized (one of the cells of cfg.gpow.comp or cfg.gpowcorr.comp))
%
%  .powsource.bline: one number in s, the center of the covariance window of the baseline (the window length depends on pow_peak)
%
%  .powsource.dics: options that will be passed to beamformer. Examples:
%     .lambda: regularization parameter of beamformer ('25%')
%     .powmethod: power method of beamformer ('trace' or 'lambda1')
%     .refdip: location of the dipole for computing coherence to.
%
%  .powsource.keepfilter: logical, to keep filters or not (keep them only
%                          if you plan to use powstat or conn analyses)
%
% IN:
%  data in /PROJ/subjects/SUBJ/MOD/NICK/
%
% OUT
%  [cfg.dpow 'powsource_SUBJ_COND'] 'powsource_s_A': source data for period of interest for each subject
%  [cfg.dpow 'powsource_SUBJ_COND'] 'powsource_s_B': source data for baseline for each subject
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

pow_peak = getpeak(cfg, 'pow');

%-------------------------------------%
%-loop over conditions
for k = 1:numel(cfg.powsource.cond)
  cond     = cfg.powsource.cond{k};
  condname = regexprep(cond, '*', '');

  %---------------------------%
  %-read data
  [data badchan] = load_data(cfg, subj, cond);
  if isempty(data)
    output = sprintf('%sCould not find any file for condition %s\n', ...
      output, cond);
    continue
  end
  
  outputfile = sprintf('powsource_%04d_%s', subj, condname);
  %---------------------------%

  %---------------------------%
  %-remove bad channels from leadfield
  datachan = ft_channelselection([{'all'}; cellfun(@(x) ['-' x], badchan, 'uni', false)], data.label);
  [leadchan] = prepare_leadchan(lead, datachan);
  %---------------------------%
  
  for p = 1:numel(pow_peak)
    
    fprintf('\n   ->->-> Running % 2d pow_peak (%s) <-<-<-\n', p, pow_peak(p).name);
    
    %---------------------------%
    %-more precise freq analysis reconstruction
    freqparam = prepare_freqpeak(cfg, pow_peak(p), data.time{1}(1));
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
    cfg1.dics = cfg.powsource.dics;
    
    cfg1.vol = vol;
    cfg1.grid = leadchan;
    cfg1.elec = sens;
    
    if cfg.powsource.keepfilter
      cfg1.dics.keepfilter   = 'yes';
      cfg1.dics.realfilter   = 'yes';
      if isfield(cfg1.dics, 'refdip')
        cfg1.dics = rmfield(cfg1.dics, 'refdip');
      end
    end
    
    powsource_s_B{p} = ft_sourceanalysis(cfg1, freq);
    powsource_s_B{p}.cfg = [];
    %-----------------%
    
    %-----------------%
    %-load MNI grid
    if ~strcmp(cfg.vol.type, 'template') ...
        && isfield(cfg, 'bnd2lead') && isfield(cfg.bnd2lead, 'mni') ...
        && isfield(cfg.bnd2lead.mni, 'warp') && cfg.bnd2lead.mni.warp
      
      load(sprintf('%s/template/sourcemodel/standard_grid3d%dmm.mat', ...
        fileparts(which('ft_defaults')), cfg.bnd2lead.mni.resolution), 'grid'); 
      
      grid = ft_convert_units(grid, 'mm');
      powsource_s_B{p}.pos = grid.pos;
    end
    %-----------------%
    %---------------------------%

    %---------------------------%
    %-main analysis
    %-----------------%
    cfg1.latency = freqparam.time;
    powsource_s_A{p} = ft_sourceanalysis(cfg1, freq);
    chan = powsource_s_A{p}.cfg.channel;
    powsource_s_A{p}.cfg = [];
    powsource_s_A{p}.cfg.channel = chan;
    %-----------------%
    
    %-----------------%
    %-load MNI grid
    if ~strcmp(cfg.vol.type, 'template') ...
        && isfield(cfg, 'bnd2lead') && isfield(cfg.bnd2lead, 'mni') ...
        && isfield(cfg.bnd2lead.mni, 'warp') && cfg.bnd2lead.mni.warp
      powsource_s_A{p}.pos = grid.pos;
    end
    %-----------------%
    %---------------------------%
    
  end

  %-----------------%
  %-save source
  save([cfg.dpow outputfile], 'powsource_s_A', 'powsource_s_B', '-v7.3')
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

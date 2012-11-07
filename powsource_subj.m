function powsource_subj(info, opt, subj)
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
%  .vol.type: 'template' or subject-specific ('dipoli' or 'openmeeg' or 'bemcp')
%  (if cfg.vol.type == 'template')
%      .vol.template: file with template containing vol, lead, sens
%
%  .sourcespace: 'surface' 'volume' 'volume_warp'
%  (if cfg.sourcespace == 'surface')
%  .SUBJECTS_DIR: where the Freesurfer data is stored (like the environmental variable)
%
%  .powsource.areas: how to specify peaks to analyze, 'manual' or 'pow_peak'
%          (peaks from grandpow) or 'powcorr_peak' (peaks from grandpowcorr)
%    if 'manual'
%      .powsource.pow_peak(1).name: string ('name_of_the_time_window')
%      .powsource.pow_peak(1).time: scalar (center of the time window in s)
%      .powsource.pow_peak(1).wndw: scalar (length of the time window in s)
%      .powsource.pow_peak(1).freq = 10; % center of the frequency
%      .powsource.pow_peak(1).band = 4; % width of the frequency band
%    if 'pow_peak' or 'powcorr_peak'
%      .powsource.refcomp: cell of string(s) of the comparison whose peaks
%                    will be localized (one of the cells of cfg.gpow.comp or cfg.gpowcorr.comp))
%
%  .powsource.bline: one number in s, the center of the covariance window
%                    of the baseline (the window length depends on pow_peak)
%                    If empty, no baseline.
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
%  [info.dpow 'powsource_SUBJ_COND'] 'powsource_s_A': source data for period of interest for each subject
%  [info.dpow 'powsource_SUBJ_COND'] 'powsource_s_B': source data for baseline for each subject
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
%-default cfg
if ~isfield(cfg.powsource, 'dics'); cfg.powsource.dics = []; end
if ~isfield(cfg.powsource, 'keepfilter')
  cfg.powsource.keepfilter = false;
end
%---------------------------%

%---------------------------%
%-dir and files
[vol, lead, sens] = load_headshape(info, opt, subj);
%---------------------------%

pow_peak = get_peak(cfg, 'pow');

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
  
  clear powsource_s_A powsource_s_B
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
    [freqparam outtmp] = prepare_freqpeak(cfg, pow_peak(p), data.time{1}(1));
    output = [output outtmp];
    %---------------------------%
    
    %---------------------------%
    %-parameters
    %-----------------%
    %-covariance window
    cfgpow = [];
    cfgpow.method = 'mtmconvol';
    cfgpow.output = 'fourier';
    cfgpow.t_ftimwin = freqparam.wndw;
    
    cfgpow.foi = freqparam.freq;
    if freqparam.dpss
      cfgpow.taper = 'dpss';
      cfgpow.tapsmofrq = freqparam.band;
    else
      cfgpow.taper = 'hanning';
    end
    
    cfgpow.feedback = 'none';
    cfgpow.channel = datachan;
    
    cfgpow.toi = [cfg.powsource.bline freqparam.time];
    %-----------------%
    
    %-----------------%
    %-source analysis
    cfgsou = [];
    
    cfgsou.frequency = freqparam.freq;
    
    cfgsou.method = 'dics';
    cfgsou.dics.feedback = 'none';
    cfgsou.dics = cfg.powsource.dics;
    
    cfgsou.vol = vol;
    cfgsou.grid = leadchan;
    cfgsou.elec = sens;
    
    if cfg.powsource.keepfilter
      cfgsou.dics.keepfilter = 'yes';
      if isfield(cfgpow.dics, 'refdip')
        cfgsou.dics = rmfield(cfgsou.dics, 'refdip');
      end
    end
    %-----------------%
    %---------------------------%
    
    %---------------------------%
    %-freq analysis (for baseline and main together)
    freq = ft_freqanalysis(cfgpow, data);
    %---------------------------%
    
    %---------------------------%
    %-baseline
    if ~isempty(cfg.powsource.bline)
      
      %-----------------%
      cfgsou.latency = cfg.powsource.bline;
      source = ft_sourceanalysis(cfgsou, freq);
      source.cfg = [];
      %-----------------%
      
      %-----------------%
      %-realign source
      powsource_s_B(p,:) = realign_source(cfg, subj, source);
      %-----------------%
      
    else
      powsource_s_B = [];
      
    end
    %---------------------------%

    %---------------------------%
    %-main analysis
    %-----------------%
    cfgsou.latency = freqparam.time;
    source = ft_sourceanalysis(cfgsou, freq);
    chan = source.cfg.channel;
    source.cfg = [];
    source.cfg.channel = chan;
    %-----------------%
    
    %-----------------%
    %-realign source
    powsource_s_A(p,:) = realign_source(cfg, subj, source);
    %-----------------%
    %---------------------------%
    
  end

  %-----------------%
  %-save source
  save([info.dpow outputfile], 'powsource_s_A', 'powsource_s_B', '-v7.3')
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

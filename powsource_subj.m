function powsource_subj(info, opt, subj)
%POWSOURCE_SUBJ: identify sources of POW using DICS for each subject
%
% INFO
%  .log: name of the file and directory to save log
%  .dpow: directory for POW data
%
% OPT
%  .cond*: cell with conditions (e.g. {'*cond1' '*cond2'})
%  .bline*: one number in s, the center of the covariance window of the
%          baseline (the window length depends on pow_peak). If empty, no
%          baseline.
%  .dics: options that will be passed to beamformer. Examples:
%     .lambda: regularization parameter of beamformer ('25%')
%     .powmethod: power method of beamformer ('trace' or 'lambda1')
%     .refdip: location of the dipole for computing coherence to.
%  .keepfilter: keep filters or not, keep them only if you plan to use
%               powstat or conn analyses (logical) 
%
% IN:
%  data in /PROJ/subjects/SUBJ/MOD/NICK/
%
% OUT
%  [info.dpow 'powsource_SUBJ_COND'] 'powsource_s_A': source data for period of interest for each subject
%  [info.dpow 'powsource_SUBJ_COND'] 'powsource_s_B': source data for baseline for each subject
%
% * indicates obligatory parameter
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
%-default opt
if ~isfield(opt, 'dics'); opt.dics = []; end
if ~isfield(opt, 'keepfilter'); opt.keepfilter = false; end
%---------------------------%

%---------------------------%
%-dir and files
[vol, lead, sens] = load_headshape(info, subj);
%---------------------------%

pow_peak = get_peak(info, opt.peak, 'pow');

%-------------------------------------%
%-loop over conditions
for k = 1:numel(opt.cond)
  cond = opt.cond{k};
  condname = regexprep(cond, '*', '');

  %---------------------------%
  %-read data
  [data badchan] = load_data(info, subj, cond);
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
    [freqparam outtmp] = prepare_freqpeak(opt, pow_peak(p), data.time{1}(1));
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
    
    cfgpow.toi = [opt.bline freqparam.time];
    %-----------------%
    
    %-----------------%
    %-source analysis
    cfgsou = [];
    
    cfgsou.frequency = freqparam.freq;
    
    cfgsou.method = 'dics';
    cfgsou.dics.feedback = 'none';
    cfgsou.dics = opt.dics;
    
    cfgsou.vol = vol;
    cfgsou.grid = leadchan;
    cfgsou.elec = sens;
    
    if opt.keepfilter
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
    if ~isempty(opt.bline)
      
      %-----------------%
      cfgsou.latency = opt.bline;
      source = ft_sourceanalysis(cfgsou, freq);
      source.cfg = [];
      %-----------------%
      
      %-----------------%
      %-realign source
      powsource_s_B(p,:) = realign_source(info, subj, source);
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
    powsource_s_A(p,:) = realign_source(info, subj, source);
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

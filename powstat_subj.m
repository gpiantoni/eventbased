function powstat_subj(cfg, subj)
%POWSTAT_SUBJ: use DICS common filters for each condition
% Notice that this function reuses the same filters for all the conditions
% but the beamforming depends on the CSD, which is different across
% conditions. You'd need to rewrite the function, where you load all the
% data at once and reseparate conditions after beamforming and rawtrial
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
%  .powsource.refcond: string of the condition used for DICS-filter
%  .powstat.cond: cell with conditions (e.g. {'*cond1' '*cond2'})
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
%      .powsource.pow_peak(1).freq: scalar (center of the frequency)
%      .powsource.pow_peak(1).band: scalar (total width of the frequency band)
%    if 'pow_peak' or 'powcorr_peak'
%      .powsource.refcomp: cell of string(s) of the comparison whose peaks
%                     will be localized (one of the cells of cfg.gpow.comp or cfg.gpowcorr.comp)
%
%  .powsource.bline: one number in s, the center of the covariance window of the baseline (the window length depends on pow_peak)
%
%  .powsource.dics: options that will be passed to beamformer. Examples:
%     .lambda: regularization parameter of beamformer ('25%')
%     .powmethod: power method of beamformer ('trace' or 'lambda1')
%
% IN:
%  data in /PROJ/subjects/SUBJ/MOD/NICK/
%  [cfg.dpow 'powsource_SUBJ_COND'] 'powsource_s_A': source data for period of interest for each subject
%  [cfg.dpow 'powsource_SUBJ_COND'] 'powsource_s_B': source data for baseline for each subject
%
% OUT
%  [cfg.dpow 'powstat_SUBJ_COND'] 'powstat_s_A': source data for period of interest for each subject, after common filters
%  [cfg.dpow 'powstat_SUBJ_COND'] 'powstat_s_B': source data for baseline for each subject, after common filters
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

%-----------------%
%-load source
souname = regexprep(cfg.powsource.refcond, '*', '');
sourcefile = sprintf('powsource_%04d_%s', subj, souname);
load([cfg.dpow sourcefile], 'powsource_s_A', 'powsource_s_B')
souchan = powsource_s_A{1}.cfg.channel;
%-----------------%

pow_peak = getpeak(cfg, 'pow');
%---------------------------%

%-------------------------------------%
%-loop over statistics conditions
for t = 1:numel(cfg.powstat.comp)
  
  %---------------------------%
  %-statistics for effects of interest
  if numel(cfg.powstat.comp{t}) == 1
    
    %-----------------%
    %-compare against baseline
    cond = cfg.powstat.comp{t}{1};
    comp = regexprep(cond, '*', '');
    output = sprintf('%s\n   COMPARISON %s\n', output, cond);
    
    %-------%
    %-load data
    [outtmp data] = load_data(cfg, subj, cond);
    output = [output outtmp];
    if isempty(data); continue; end
    %-------%
    
    design = ones(numel(data.trial),1);
    %-----------------%
    
  else
    
    %-----------------%
    %-compare two conditions (cond1 - cond2)
    cond1 = cfg.powstat.comp{t}{1};
    cond2 = cfg.powstat.comp{t}{2};
    comp = [regexprep(cond1, '*', '') '_' regexprep(cond2, '*', '')];
    output = sprintf('%s\n   COMPARISON %s vs %s\n', output, cond1, cond2);
    
    %-------%
    %-load data
    [outtmp data1] = load_data(cfg, subj, cond1);
    output = [output outtmp];
    if isempty(data1); continue; end
    %-------%
    
    %-------%
    %-load data
    [outtmp data2] = load_data(cfg, subj, cond2);
    output = [output outtmp];
    if isempty(data2); continue; end
    %-------%
    
    data = ft_appenddata([], data1, data2);
    
    design = [ones(numel(data1.trial),1); 2 * ones(numel(data2.trial),1)];
    %-----------------%
    
  end
  
  outputfile = sprintf('powstat_%04d_%s', subj, comp)
  clear data1 data2
  %---------------------------%
  
  %-------------------------------------%
  %-loop over conditions
  for k = 1:numel(cfg.powstat.cond)
    cond     = cfg.powstat.cond{k};
    condname = regexprep(cond, '*', '');
    
    %---------------------------%
    %-read data
    [data badchan] = load_data(cfg, subj, cond);
    if isempty(data)
      output = sprintf('%sCould not find any file for test %s\n', ...
        output, cond);
      continue
    end
    
    ;
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
  
  for p = 1:numel(pow_peak)
    
    fprintf('\n   ->->-> Running % 2d powstat (%s) <-<-<-\n', p, pow_peak(p).name);
    
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
    
    cfg1.grid.filter  = powsource_s_B{p}.avg.filter;
    
    powstat_s_B{p} = ft_sourceanalysis(cfg1, freq);
    powstat_s_B{p}.cfg = [];
    powstat_s_B{p}.cfg.ntapers = sum(freq.cumtapcnt);
    %-----------------%
    
    %-----------------%
    %-load MNI grid
    if ~strcmp(cfg.vol.type, 'template') ...
        && isfield(cfg, 'bnd2lead') && isfield(cfg.bnd2lead, 'mni') ...
        && isfield(cfg.bnd2lead.mni, 'warp') && cfg.bnd2lead.mni.warp
      
      load(sprintf('%s/template/sourcemodel/standard_grid3d%dmm.mat', ...
        fileparts(which('ft_defaults')), cfg.bnd2lead.mni.resolution), 'grid');
      
      grid = ft_convert_units(grid, 'mm');
      powstat_s_B{p}.pos = grid.pos;
    end
    %-----------------%
    %---------------------------%
    
    %---------------------------%
    %-main analysis
    %-----------------%
    cfg1.latency = freqparam.time;
    cfg1.grid.filter = powsource_s_A{p}.avg.filter;
    powstat_s_A{p} = ft_sourceanalysis(cfg1, freq);
    chan = powstat_s_A{p}.cfg.channel;
    powstat_s_A{p}.cfg = [];
    powstat_s_A{p}.cfg.channel = chan;
    powstat_s_A{p}.cfg.ntapers = sum(freq.cumtapcnt);
    %-----------------%
    
    %-----------------%
    %-load MNI grid
    if ~strcmp(cfg.vol.type, 'template') ...
        && isfield(cfg, 'bnd2lead') && isfield(cfg.bnd2lead, 'mni') ...
        && isfield(cfg.bnd2lead.mni, 'warp') && cfg.bnd2lead.mni.warp
      powstat_s_A{p}.pos = grid.pos;
    end
    %-----------------%
    %---------------------------%
    
  end
  
  %-----------------%
  %-save source
  save([cfg.dpow outputfile], 'powstat_s_A', 'powstat_s_B', '-v7.3')
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
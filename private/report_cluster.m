function [clpeak stat output] = reportcluster(cfg, gdat, gdat2, paired)
%REPORTCLUSTER report which clusters are significant
% It compares against zero or two conditions. The significant
% time-(freq)-sensor points are clustered and then statistics is run on
% these clusters. Significant clusters will be reported and can be used for
% further analysis, for example, for source reconstruction.
%
% Not to be called directly, this function is called by ERP_GRAND,
% POW_GRAND and POWCORR_GRAND.
%
% Use as:
%   [clpeak output] = reportcluster(cfg, gdat, gdat2, unpaired)
%
% CFG
%  .sens.file: file with EEG sensors. It can be sfp or mat.
%  .sens.dist: distance of two sensors to be considered neighbors
%
%  You can then specify the time window and frequency band to do statistics on:
%  .gerp.stat.time: for ERP only, two scalars for the beginning and end of the time window
%
%  .gpow.stat.time: for POW only, two scalars for the beginning and end of the time window
%  .gpow.stat.freq: for POW only, two scalars for the beginning and end of the frequency band
% 
%  .gpowcorr.stat.time: for POWCORR only, two scalars for the beginning and end of the time window
%  .gpowcorr.stat.freq: for POWCORR only, two scalars for the beginning and end of the frequency band
% 
%  .cluster.thr: threshold 
%
%  .pow:
%  
% GDAT, GDAT2
%  One or two datasets obtained from ft_timelockgrandaverage or
%  ft_freqgrandaverage.
%  If one dataset, the dataset is compared against zero
%  If two datasets, the two datasets are compared against each other in a
%  PAIRED or UNPAIRED t-test.
%  More complicated designs are not possible.
%
% PAIRED
%  Flag for PAIRED (true) or UNPAIRED (false) design. Default is true
%
% CLPEAK
%  which is the cluster peak, with fields:
%    .name: the name, based on the peak time for ERP data or on the peak
%           frequency and peak time for POW data.
%    .time: center of the time window (not the peak)
%    .wndw: length of the time window
%    .freq: center of the frequency band (not the peak)
%    .band: width in frequency of the frequency band (take half of that for dpss)
% The sizes of the frequency band and time window are defined by FWHM
%
% Note that the name uses the peak frequency and peak time window, which
% should not be very different from .freq and .time 
% This function is not very robust against one big significant cluster
% which connects two blobs (beta and gamma for example). This can seen when
% the peak freq (in the name) is very different from .freq
%
% Part of EVENTBASED/PRIVATE

%-----------------%
%-check CFG
if ~isfield(cfg, 'cluster') || ~isfield(cfg.cluster, 'thr')
  cfg.cluster.thr = 0.05;
end
if nargin < 4
  paired = true;
end
%-----------------%

%-----------------%
%-check data
if ~isfield(gdat, 'powspctrm')
  iserp = true;
  param = 'individual';
else
  iserp = false;
  param = 'powspctrm';
end

output = '';
nsubj1 = size(gdat.(param),1);

if paired
  nsubj2 = nsubj1;
else
  nsubj2 = size(gdat2.(param),1);
end

if nargin == 2
  gdat2 = gdat;
  gdat2.(param) = zeros(size(gdat.(param)));
end
%-----------------%

%-----------------%
%-sensors
if ~isempty(cfg.sens.file)
  
  %-------%
  %-create neighbors from file
  sens = ft_read_sens(cfg.sens.file);
  sens.label = upper(sens.label);
  
  tmpcfg = [];
  tmpcfg.elec = sens;
  tmpcfg.method = 'distance';
  tmpcfg.neighbourdist = cfg.sens.dist;
  neigh = ft_prepare_neighbours(tmpcfg);
  %-------%
  
else
  
  %-------%
  %-create fake neighbors (no real neighbors)
  neigh = [];
  for i = 1:numel(gdat.label)
    neigh(i).label = gdat.label{i};
    neigh(i).neighblabel(1) = gdat.label(i);
  end
  %-------%
  
end
%-----------------%

%-----------------%
%-calc clusters
tmpcfg = [];
tmpcfg.method      = 'montecarlo';
if paired
  tmpcfg.statistic   = 'depsamplesT';
else
  tmpcfg.statistic   = 'indepsamplesT';
end
tmpcfg.alpha       = 0.05;
tmpcfg.correctm    = 'cluster';
tmpcfg.numrandomization = 100; 

if paired
  tmpcfg.design = [ones(1,nsubj1) ones(1,nsubj2)*2; 1:nsubj1 1:nsubj2];
  tmpcfg.ivar   = 1;
  tmpcfg.uvar   = 2;
else
  tmpcfg.design = [ones(1,nsubj1) ones(1,nsubj2)*2];
  tmpcfg.ivar   = 1;
end

if iserp && isfield(cfg.gerp, 'stat') && isfield(cfg.gerp.stat, 'time') && ~isempty(cfg.gerp.stat.time)
  tmpcfg.latency = cfg.gerp.stat.time;
elseif ~iserp && isfield(cfg.gpow, 'stat') && isfield(cfg.gpow.stat, 'time') && ~isempty(cfg.gpow.stat.time)
  tmpcfg.latency = cfg.gpow.stat.time;
else
  tmpcfg.latency = gdat.time([1 end]);
end

if ~iserp && isfield(cfg.gpow, 'stat') && isfield(cfg.gpow.stat, 'freq') && ~isempty(cfg.gpow.stat.freq)
  tmpcfg.frequency = cfg.gpow.stat.freq;
end

tmpcfg.neighbours = neigh;
tmpcfg.feedback = 'etf';

if iserp
  % tmpcfg.minnbchan = 5; % to avoid huge clusters
  stat = ft_timelockstatistics(tmpcfg, gdat, gdat2); 
else
  stat = ft_freqstatistics(tmpcfg, gdat, gdat2);
end
%-----------------%

%-----------------%
%-if there are no clusters at all
if ~isfield(stat, 'posclusters') || ~isfield(stat, 'negclusters') 
  clpeak = [];
  output = sprintf('no clusters at all in the data\n');
  return
end  
%-----------------%

%-------------------------------------%
%-report cluster
clpeak = [];

%-----------------%
%-positive cluster
if isempty(stat.posclusters)
  poscl = [];
else
  poscl = find([stat.posclusters.prob] < cfg.cluster.thr);
end

for i = 1:numel(poscl)
  
  %-------%
  %-find important elements of cluster
  clmat = stat.posclusterslabelmat == i;
  clst  = stat.stat .* clmat;
  
  if iserp
    cltime = squeeze(nansum(clst, 1));
    
  else
    cltime = squeeze(nansum(nansum(clst, 1),2));

    clfreq = squeeze(nansum(nansum(clst, 1),3));
    [~, maxfreq] = max(clfreq);
    [fbeg fend] = fwhm(clfreq);
    fbeg = stat.freq(fbeg); % index -> absolute value
    fend = stat.freq(fend); % index -> absolute value
    
  end
  
  [~, maxtime] = max(cltime);
  [tbeg tend] = fwhm(cltime);
  tbeg = stat.time(tbeg); % index -> absolute value
  tend = stat.time(tend); % index -> absolute value
  %-------%
  
  %-------%
  %-write description
  if iserp
    outtmp = sprintf('pos cluster % 2.f: P = %4.3f, size =% 6.f, time (% 3.2f) [% 3.2f % 3.2f]\n', ...
      i, stat.posclusters(i).prob, numel(find(clmat)), ...
      stat.time(maxtime), tbeg, tend);
    
  else
    outtmp = sprintf('pos cluster % 2.f: P = %4.3f, size =% 6.f, freq (% 5.1f) [% 5.1f % 5.1f], time (% 3.2f) [% 3.2f % 3.2f]\n', ...
      i, stat.posclusters(i).prob, numel(find(clmat)), ...
      stat.freq(maxfreq), fbeg, fend, ...
      stat.time(maxtime), tbeg, tend);
    
  end
  output = [output outtmp];
  %-------%
  
  %-------%
  %-prepare peak
  clpeak(i).pval = stat.posclusters(i).prob;
  clpeak(i).time = mean([tbeg tend]);
  clpeak(i).wndw = tend - tbeg;
  
  if iserp
    clpeak(i).name = sprintf('pos_at%04.f', stat.time(maxtime)*1e3);
    
  else
    clpeak(i).freq = mean([fbeg fend]);
    clpeak(i).band = fend - fbeg;
    clpeak(i).name = sprintf('pos_freq%02.fat%04.f', stat.freq(maxfreq), stat.time(maxtime)*1e3);
    
  end
  %-------%
end

if isempty(i)
  ipos = 0;
else
  ipos = i;
end
%-----------------%

%-----------------%
%-negative cluster
if isempty(stat.negclusters)
  negcl = [];
else
  negcl = find([stat.negclusters.prob] < cfg.cluster.thr);
end

for i = 1:numel(negcl)
  
  %-------%
  %-find important elements of cluster
  clmat = stat.negclusterslabelmat == i;
  clst  = -1 * stat.stat .* clmat;
  
  if iserp
    cltime = squeeze(nansum(clst, 1));
    
  else
    cltime = squeeze(nansum(nansum(clst, 1),2));

    clfreq = squeeze(nansum(nansum(clst, 1),3));
    [~, maxfreq] = max(clfreq);
    [fbeg fend] = fwhm(clfreq);
    fbeg = stat.freq(fbeg); % index -> absolute value
    fend = stat.freq(fend); % index -> absolute value
    
  end
  
  [~, maxtime] = max(cltime);
  [tbeg tend] = fwhm(cltime);
  tbeg = stat.time(tbeg); % index -> absolute value
  tend = stat.time(tend); % index -> absolute value
  %-------%
  
  %-------%
  %-write description
  if iserp
    outtmp = sprintf('neg cluster % 2.f: P = %4.3f, size =% 6.f, time (% 3.2f) [% 3.2f % 3.2f]\n', ...
      i, stat.negclusters(i).prob, numel(find(clmat)), ...
      stat.time(maxtime), tbeg, tend);
    
  else
    outtmp = sprintf('neg cluster % 2.f: P = %4.3f, size =% 6.f, freq (% 5.1f) [% 5.1f % 5.1f], time (% 3.2f) [% 3.2f % 3.2f]\n', ...
      i, stat.negclusters(i).prob, numel(find(clmat)), ...
      stat.freq(maxfreq), fbeg, fend, ...
      stat.time(maxtime), tbeg, tend);
    
  end
  output = [output outtmp];
  %-------%
  
  %-------%
  %-prepare peak
  clpeak(i+ipos).pval = stat.negclusters(i).prob;
  clpeak(i+ipos).time = mean([tbeg tend]);
  clpeak(i+ipos).wndw = tend - tbeg;
  
  if iserp
    clpeak(i+ipos).name = sprintf('neg_at%04.f', stat.time(maxtime)*1e3);
    
  else
    clpeak(i+ipos).freq = mean([fbeg fend]);
    clpeak(i+ipos).band = fend - fbeg;
    clpeak(i+ipos).name = sprintf('neg_freq%02.fat%04.f', stat.freq(maxfreq), stat.time(maxtime)*1e3);
    
  end
  %-------%
end
%-----------------%

%-----------------%
%-add time window used for the calculation of power into clpeak
% This time window is extremely tricky to work with: powpeak is
% selected based on the significant time points of the TFR. However,
% each of those points represents an FFT around that time window, with
% length defined by cfg.pow.t_ftimwin. powpeak(f).wndw reports only the
% significant timepoints (it can 0, if only one time point is
% significant), so we add cfg.pow.t_ftimwin here (half in the beginning
% and half in the end, but it's centered around powpeak(f).time anyway)
if ~iserp
  for f = 1:numel(clpeak)
    ifreq = nearest(cfg.pow.foi, clpeak(f).freq); % length of cfg.pow.t_ftimwin depends on the frequency
    if strcmp(cfg.pow.method, 'wavelet') || strcmp(cfg.pow.method, 'tfr');
      clpeak(f).wndw = clpeak(f).wndw + cfg.pow.width/cfg.pow.foi(ifreq); 
    elseif strcmp(cfg.pow.method, 'mtmconvol');
      clpeak(f).wndw = clpeak(f).wndw + cfg.pow.t_ftimwin(ifreq); 
    end
  end
end
%-----------------%
%-------------------------------------%

%-------------------------------------%
%-SUBFUNCTION: fwhm
%-------------------------------------%
function [tbeg tend] = fwhm(clsum)
%FWHM: find the full witdh at half maximum of time and frequency
hm = max(clsum)/2;

tbeg = find(clsum >= hm, 1, 'first'); % from left 
tend = find(clsum(end:-1:1) >= hm, 1, 'first'); % from right
tend = numel(clsum) - tend + 1;

% figure
% plot(clsum)
% hold on;
% plot(tbeg, hm, 'r+')
% plot(tend, hm, 'r+')
%-------------------------------------%
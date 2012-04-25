function [clpeak output] = reportcluster(gdat, cfg)
%REPORTCLUSTER get cluster which are different from zero, even if not significant
% The clusters to determine the main results of the analysis, for example
% to concentrate the source reconstruction
% It returns clpeak
%   .freq = center of the frequency band
%   .band = size in frequency of the frequency band (take half of that for dpss)
%   .time = center of the time window
%   .wndw = length of the time window
%   .name = the name, based on the peak frequency and peak time
% Note that the name refers to the peak frequency and peak time window,
% they should not be very different from .freq and .time
% This function is not very robust against one big significant cluster
% which connects two blobs (beta and gamma for example). This can seen when
% the peak freq (in the name) is very different from .freq
% However, peak detection is very hard. It now relies on cluster analysis,
% it might be made more flexible.
%
% The sizes of the frequency band and time window are defined by FWHM

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
nsubj = size(gdat.(param),1);

gzero = gdat;
gzero.(param) = zeros(size(gdat.(param)));
%-----------------%

%-----------------%
%-sensors
if ~isempty(cfg.sens.file)
  
  %-------%
  %-create neighbors from file
  sens = ft_read_sens(cfg.sens.file);
  sens.label = upper(sens.label);
  
  cfg1 = [];
  cfg1.elec = sens;
  cfg1.method = 'distance';
  cfg1.neighbourdist = cfg.sens.dist;
  neigh = ft_prepare_neighbours(cfg1);
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
cfg3 = [];
cfg3.method      = 'montecarlo';
cfg3.statistic   = 'depsamplesT';
cfg3.alpha       = 0.05;
cfg3.correctm    = 'cluster';
cfg3.numrandomization = 5000; 
cfg3.design = [ones(1,nsubj) ones(1,nsubj).*2; 1:nsubj 1:nsubj];
cfg3.ivar   = 1;
cfg3.uvar   = 2;

if iserp && isfield(cfg.gerp, 'time') && ~isempty(cfg.gerp.time)
  cfg3.latency = cfg.gerp.time;
elseif ~iserp && isfield(cfg.gpow, 'time') && ~isempty(cfg.gpow.time)
  cfg3.latency = cfg.gpow.time;
else
  cfg3.latency = gdat.time([1 end]);
end

if ~iserp && isfield(cfg.gpow, 'freq') && ~isempty(cfg.gpow.freq)
  cfg3.frequency = cfg.gpow.freq;
end

cfg3.neighbours = neigh;
cfg3.feedback = 'etf';

if iserp
  % cfg3.minnbchan = 5; % to avoid huge clusters
  stat = ft_timelockstatistics(cfg3, gdat, gzero); 
else
  stat = ft_freqstatistics(cfg3, gdat, gzero);
end
%-----------------%

%-----------------%
%-if there are no clusters at all
if ~isfield(stat, 'posclusters') 
  clpeak = [];
  output = sprintf('no clusters at all in the data\n');
  return
end  
%-----------------%

%-------------------------------------%
%-report cluster
clpeak = [];
clusterthr = .9;

%-----------------%
%-positive cluster
if isempty(stat.posclusters)
  poscl = [];
else
  poscl = find([stat.posclusters.prob] < clusterthr);
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
  negcl = find([stat.negclusters.prob] < clusterthr);
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
    if isfield(cfg.pow.method, 'wavelet') || isfield(cfg.pow.method, 'tfr');
      clpeak(f).wndw = clpeak(f).wndw + cfg.pow.width/cfg.pow.foi(ifreq); 
    elseif isfield(cfg.pow.method, 'mtmconvol');
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
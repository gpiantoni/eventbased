function [stat output] = report_source(cfg, gdat1, gdat2)
%REPORT_SOURCE get clusters in the comparison between two conditions or against baseline
%
% It only reports either positive or negative clusters. It does not make
% sense to have both (how can one TFR element be associated with activation
% and disactivation from baseline?)  
% 
% Use as:
%   [stat output] = report_source(cfg, gdat1, gdat2)
%
% CFG
%  .clusterstatistics: 'maxsize' or 'max'
%  .clusteralpha: level to select sensors (default 0.05)
%                   it can be a string in format '5%' to take top 5 voxels
%                   and put them in a cluster. 
%  .numrandomization: n of randomizations (default 1e5)
%  .clusterthr: threshold to report clusters in output (default 0.5)
%
% GDAT1/GDAT2
%  grandaverage datasets for source
% The results are GDAT1 - GDAT2
% If positive, more source in GDAT1 than GDAT2
% If negative, more source in GDAT2 than GDAT1
%
% Part of EVENTBASED/PRIVATE

%-------------------------------------%
%---------------------------%
%-cfg default values
if ~isfield(cfg, 'clusterstatistics'); cfg.clusterstatistics = 'maxsize'; end
if ~isfield(cfg, 'numrandomization'); cfg.numrandomization = 1e5; end
if ~isfield(cfg, 'clusteralpha'); cfg.clusteralpha = 0.05; end
if ~isfield(cfg, 'clusterthr'); cfg.clusterthr = .5; end
if ~isfield(cfg, 'minnbchan'); cfg.minnbchan = 0; end
if ~isfield(cfg, 'atlas'); cfg.atlas = 1; end
%---------------------------%

%---------------------------%
%-check data
nsubj = numel(gdat1.trial);
output = sprintf('on field %s, with %d subjects\n', cfg.parameter, nsubj);

if nargin == 2
  gdat2 = gdat1;
  for i = 1:numel(gdat2.trial)
    gdat2.trial(i).(cfg.parameter) = zeros(size(gdat2.trial(i).(cfg.parameter)));
  end
end
%-----------------%
%---------------------------%
%-------------------------------------%

%-------------------------------------%
%-calc clusters
tmpcfg = [];
tmpcfg.method      = 'montecarlo';
if strcmp(cfg.parameter, 'coh')
  tmpcfg.statistic   = 'diff';
  tmpcfg.clusterthreshold = 'nonparametric_common';
else
  tmpcfg.statistic   = 'depsamplesT';
end
tmpcfg.correctm    = 'cluster';
tmpcfg.clusterstatistic = cfg.clusterstatistics; % 'maxsize' or 'max' ('max' might be better for focal sources)
tmpcfg.numrandomization = cfg.numrandomization;
tmpcfg.design = [ones(1,nsubj) ones(1,nsubj).*2; 1:nsubj 1:nsubj];
tmpcfg.ivar = 1;
tmpcfg.uvar = 2;
tmpcfg.feedback = 'etf';

tmpcfg.parameter = cfg.parameter;
if isfield(gdat1, 'dim')
  tmpcfg.dim = gdat1.dim;
  issurf = false;
  
else
  tmpcfg.dim = size(gdat1.avg.pow);
  tmpcfg.channeighbstructmat = cfg.channeighbstructmat;
  issurf = true;
  
end
tmpcfg.minnbchan = cfg.minnbchan;
tmpcfg.alpha = 0.05;

%---------------------------%
%-get value for clusteralpha
if isnumeric(cfg.clusteralpha)
  
  tmpcfg.clusteralpha = cfg.clusteralpha;
  
elseif ischar(cfg.clusteralpha) && cfg.clusteralpha(end)=='%'
  
  %-----------------%
  %-use adaptive threshold
  ratio = sscanf(cfg.clusteralpha, '%f%%');
  ratio = 1 - ratio/100;
  
  %-------%
  %-get t-stat
  tstat = @(x) mean(x,2) ./ (std(x,[],2)) *sqrt(size(x,2));
  x1 = cat(2, gdat2.trial.pow);
  x2 = cat(2, gdat1.trial.pow);
  
  t = tstat(x2 - x1);
  t_thr = quantile(abs(t(~isnan(t))), ratio); % t-stat at the cfg.clusteralpha quantile
  %-------%
  
  %-------%
  %-threshold as p-value
  df = size(x1,2)-1;
  tmpcfg.clusteralpha = 2 * (1 - tcdf(t_thr, df));
  %-------%
  
  %-------%
  %-output
  outtmp = sprintf('using adaptive threshold at %s, which corresponds to t-stat % 6.4f, P-value % 10.2e\n', ...
    cfg.clusteralpha, t_thr, tmpcfg.clusteralpha);
  outtmp = regexprep(outtmp, '%', '%%'); % otherwise fprint gets confused for normal % sign
  output = [output outtmp];
  %-------%
  
  %-----------------%
  
end
%---------------------------%

stat = ft_sourcestatistics(tmpcfg, gdat1, gdat2);
%-------------------------------------%

%-------------------------------------%
%-find clusters
%-----------------%
%-if there are no clusters at all
if ~isfield(stat, 'posclusters')
  soupeak = [0 0 0];
  output = sprintf('no clusters at all in the data\n');
  stat.image = zeros(stat.dim);
  return
end
%-----------------%

%-----------------%
%-report cluster
%-------%
%-find out if more positive or negative
% it does not make sense to have both increase and decrease in power for
% the same frequency
if isempty(stat.posclusters)
  posneg = false;
elseif isempty(stat.negclusters)
  posneg = true;
else
  posneg = stat.posclusters(1).prob <= stat.negclusters(1).prob;
end

if posneg
  posnegtxt = 'pos';
  clusters = stat.posclusters;
  clusterslabelmat = stat.posclusterslabelmat;
else
  posnegtxt = 'neg';
  clusters = stat.negclusters;
  clusterslabelmat = stat.negclusterslabelmat;
end
%-------%

%-------%
%-report significant cluster
signcl = find([clusters.prob] < cfg.clusterthr);

for i = 1:numel(signcl)
  
  clmat = find(clusterslabelmat == i);
  
  if issurf
    outtmp = sprintf('    %s cluster% 3.f: P = %4.3f, size =% 5.f\n', ...
      posnegtxt, i, clusters(i).prob, numel(clmat));
    
  else
    pos = stat.pos(clmat,:);
    atlas = mni2ba(pos, cfg.atlas);
    atlasname = fieldnames(atlas);
    atlas = {atlas.(atlasname{1})};
    atlas = atlas(~strcmp(atlas, ''));
    
    if isempty(atlas)
      s_atlas = '';
    else
      [atlas, ~, i_atl] = unique(atlas);
      s_atlas = atlas{mode(i_atl)};
    end
    
    outtmp = sprintf('    %s cluster% 3.f: P = %4.3f, size =% 5.f, [% 5.1f % 5.1f % 5.1f], %s\n', ...
      posnegtxt, i, clusters(i).prob, numel(clmat), mean(pos(:,1)), mean(pos(:,2)), mean(pos(:,3)), s_atlas);
    
  end
  
  output = [output outtmp];
  
end
%-------%
%-----------------%
%-------------------------------------%

%-------------------------------------%
%-simplify the output
stat = ft_checkdata(stat, 'datatype', 'source'); % from volume into source
stat.posneg = posneg;
%-------------------------------------%

function [soupeak stat output] = reportsource(gdat, gpre)
%REPORTSOURCE get cluster which are different from baseline, even if not significant
% The clusters to determine the main results of the analysis, for example
% to concentrate the source reconstruction
%
% For source, it only reports either positive or negative clusters. It does
% not make sense to have both (how can one TFR element be associated with
% activation and disactivation from baseline?)

%06 12/02/14 export stat as well
%05 12/02/12 use absolute z-score threshold instead of cluster
%04 12/01/31 plot image with most significant voxels
%03 12/01/19 if there are no clusters at all, use [0 0 0]
%02 12/01/12 don't compare against zero, but against baseline
%01 12/01/11 created from reportcluster

addpath /data1/toolbox/helpers/ % mni2ba

maxclust = 150;

%-----------------%
%-check data
output = '';
nsubj = numel(gdat.trial);

%-------%
%-pow or nai
if isfield(gdat.avg, 'nai')
  param = 'nai';
  
elseif isfield(gdat.avg, 'pow')
  param = 'pow';
  
else
  error('field in .avg. not recognized')
end
%-------%

% gzero = gdat;
% 
% for i = 1:nsubj
%   gzero.trial(i).(param) = zeros(size(gzero.trial(i).(param)));
%   gzero.trial(i).(param)(isnan(gzero.trial(i).(param))) = NaN;
% end
%-----------------%

%-----------------%
%-calc clusters
cfg3 = [];
cfg3.method      = 'montecarlo';
cfg3.statistic   = 'depsamplesT';
cfg3.alpha       = 0.05;
cfg3.correctm    = 'cluster';
cfg3.clusterstatistic = 'maxsize'; % 'maxsize' or 'max' ('max' might be better for focal sources)
cfg3.numrandomization = 1e4;
cfg3.design = [ones(1,nsubj) ones(1,nsubj).*2; 1:nsubj 1:nsubj];
cfg3.ivar   = 1;
cfg3.uvar   = 2;
cfg3.feedback = 'none';

cfg3.parameter = param;
cfg3.dim = gdat.dim;
stat = ft_sourcestatistics(cfg3, gdat, gpre);
%-----------------%

%-----------------%
%-if there are no clusters at all
if ~isfield(stat, 'posclusters') 
  soupeak = [0 0 0];
  output = sprintf('no clusters at all in the data\n');
  return
end  
%-----------------%

%-----------------%
%-avoid huge clusters
alphaval = [0.025 0.01 0.005 0.001 0.0005 0.0001 0.00005 0.00001 0.000005 0.000001];
for alpha = alphaval
  stat_orig  = stat;
  
  if numel(find(stat.posclusterslabelmat == 1)) > maxclust ...
      || numel(find(stat.negclusterslabelmat == 1)) > maxclust
    
    cfg3.clusteralpha       = alpha;
    stat = ft_sourcestatistics(cfg3, gdat, gpre);
  else
    break
  end
  
  %-------%
  %-no clusters at all, we went too far
  if ~isfield(stat, 'posclusters')
    stat = stat_orig;
    break
  end
  %-------%
  
end
%-----------------%
  
%-----------------%
%-report cluster
clusterthr = .8;

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
%-show only first source for connectivity analysis
soupeak = stat.pos(clusterslabelmat==1, :);
%-------%

%-------%
%-report significant cluster
signcl = find([clusters.prob] < clusterthr);

for i = 1:numel(signcl)
  
  clmat = find(clusterslabelmat == i);
  pos = stat.pos(clmat,:);

  areas = mni2ba(pos);
  areastext = unique({areas.IBASPM116});
   
  outtmp = sprintf('    (a: %1.e) %s cluster% 3.f: P = %4.3f, size =% 4.f, [% 5.1f % 5.1f % 5.1f], %s\n', ...
    stat.cfg.clusteralpha, posnegtxt, i, clusters(i).prob, numel(clmat), mean(pos(:,1)), mean(pos(:,2)), mean(pos(:,3)), sprintf(' %s', areastext{:}));
  output = [output outtmp];
  
end
%-------%
%-----------------%

%-----------------%
%-plot main cluster
%-------%
%-prepare figure
backgrnd = isnan(clusterslabelmat); % separate NaN to be used as background
% clmat = clusterslabelmat == 1; % only main cluster
clmat = abs(stat.stat) > 3; % only main cluster

%-prepare axis 
xpos = unique(stat.pos(:,1));
ypos = unique(stat.pos(:,2));
zpos = unique(stat.pos(:,3));
%-------%

%-------%
%-x-axis
subplot(2,2,1)
[~, imax] = max(sum(sum(clmat,2),3));
toplot = nansum(cat(1, -1 * backgrnd(imax,:,:),  clmat(imax, :, :) .* abs(stat.stat(imax,:,:))), 1);

imagesc(ypos, zpos, squeeze(toplot)', [-.5 7])
axis xy equal
colormap hot

title(sprintf('x =% 3.f', xpos(imax)))
%-------%

%-------%
%-y-axis
subplot(2,2,2)
[~, imax] = max(sum(sum(clmat,1),3));
toplot = nansum(cat(2, -1 * backgrnd(:,imax,:),  clmat(:,imax,:) .* abs(stat.stat(:,imax,:))), 2);

imagesc(xpos, zpos, squeeze(toplot)', [-.5 7])
axis xy equal
colormap hot

title(sprintf('y =% 3.f', ypos(imax)))
%-------%

%-------%
%-z-axis
subplot(2,2,3)
[~, imax] = max(sum(sum(clmat,1),2));
toplot = nansum(cat(3, -1 * backgrnd(:,:,imax),  clmat(:,:,imax) .* abs(stat.stat(:,:,imax))), 3);

imagesc(xpos, ypos, squeeze(toplot)', [-.5 7])
axis xy equal
colormap hot

title(sprintf('z =% 3.f', zpos(imax)))
%-------%
%-----------------%
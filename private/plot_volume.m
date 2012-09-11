function plot_volume(stat, template, toplot)
%PLOT_VOLUME plot significant voxels
%
% TODO: clean-up this  function, plot image as mask

tmpcfg = [];
tmpcfg.parameter = toplot;
souinterp = ft_sourceinterpolate(tmpcfg, stat{1}, template);

tmpcfg = [];
tmpcfg.method = 'slice';
tmpcfg.funparameter = toplot;
tmpcfg.funcolorlim = [-2.5 2.5];
ft_sourceplot(tmpcfg, souinterp)

return




%-------------------------------------%
%-----------------%
%-plot main cluster
%-prepare figure
backgrnd = isnan(clusterslabelmat); % separate NaN to be used as background

%-prepare axis
xpos = unique(stat.pos(:,1));
ypos = unique(stat.pos(:,2));
zpos = unique(stat.pos(:,3));
%-----------------%


%-----------------%
%-plot
%-------%
%-x-axis
subplot(2,2,1)
[~, imax] = max(sum(sum(clmat==2,2),3));
toplot = nansum(cat(1, -1 * backgrnd(imax,:,:),  clmat(imax, :, :)), 1);

imagesc(ypos, zpos, squeeze(toplot)', [-.5 2])
axis xy equal
colormap hot

title(sprintf('x =% 3.f', xpos(imax)))
%-------%

%-------%
%-y-axis
subplot(2,2,2)
[~, imax] = max(sum(sum(clmat==2,1),3));
toplot = nansum(cat(2, -1 * backgrnd(:,imax,:),  clmat(:,imax,:)), 2);

imagesc(xpos, zpos, squeeze(toplot)', [-.5 2])
axis xy equal
colormap hot

title(sprintf('y =% 3.f', ypos(imax)))
%-------%

%-------%
%-z-axis
subplot(2,2,3)
[~, imax] = max(sum(sum(clmat==2,1),2));
toplot = nansum(cat(3, -1 * backgrnd(:,:,imax),  clmat(:,:,imax)), 3);

imagesc(xpos, ypos, squeeze(toplot)', [-.5 2])
axis xy equal
colormap hot

title(sprintf('z =% 3.f', zpos(imax)))
%-------%
%-----------------%
%-------------------------------------%

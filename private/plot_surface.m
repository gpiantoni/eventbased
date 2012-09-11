function plot_surface(stat, surfplot, toplot)
%PLOT_SURFACE plot surface
%
% TODO: allow for masking

cnt = 0;
for i = 1:2 % viewpoint
  for h = 1:numel(stat) % hemisphere
    cnt = cnt + 1;
    subplot(2,2, cnt)
    
    ft_plot_mesh(surfplot{h}, 'edgecolor', 'none', 'vertexcolor', double(stat{h}.(toplot)))
    set(gca, 'clim', [-3 3]);
    
    if ismember(cnt, [1 4])
        view(90,0);
    else
        view(-90,0);
    end
    camlight
    
  end
end



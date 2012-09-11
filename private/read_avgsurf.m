function [surfstat, surfplot] = read_avgsurf(sdir, ratio)
%READ_AVGSURF read average surface and downsample it
%
% Use as:
%   [surfstat, surfplot] = read_avgsurf(sdir, ratio)
%
% SDIR is the directory with the surfaces of fsaverage
%
% RATIO is the amount of downsampling
%
% SURFSTAT is the mesh for the sphere
% 
% SURFPLOT is the mesh for the plotting surface

hemi = {'lh' 'rh'};
surf4stat = 'sphere.reg';
surf4plot = 'pial';

%-------------------------------------%
%-loop over hemisphere
for h = 1:numel(hemi)

  %-----------------%
%-reduce one stat mesh
  surfbase = ft_read_headshape([sdir hemi{h} '.' surf4stat]);
  surfstat{h} = reducemesh(surfbase, ratio);
  [~, i_full, i_reduced] = intersect(surfbase.pnt, surfstat{h}.pnt, 'rows');
  surfstat{h}.inside = true(size(surfstat{h}.pnt,1),1);
  %-----------------%

    %-----------------%
  %-create neighbor structure
  tri = surfstat{h}.tri;
  
  ndip = size(surfstat{h}.pnt,1);
  neigh = single(zeros(ndip, ndip));

  % mark neighbours according to triangulation
  for i=1:size(tri, 1)
    neigh(tri(i, 1), tri(i, 2)) = 1;
    neigh(tri(i, 1), tri(i, 3)) = 1;
    neigh(tri(i, 2), tri(i, 1)) = 1;
    neigh(tri(i, 3), tri(i, 1)) = 1;
    neigh(tri(i, 2), tri(i, 3)) = 1;
    neigh(tri(i, 3), tri(i, 2)) = 1;
  end
  surfstat{h}.neigh = neigh;
  %-----------------%
  
  %-----------------%
  %-use the same reduction
  surfbase = ft_read_headshape([sdir hemi{h} '.' surf4plot]);
  surfplot{h}.pnt(i_reduced,:) = surfbase.pnt(i_full,:);
  surfplot{h}.tri = surfstat{h}.tri;
  %-----------------%

end
%-------------------------------------%

%-------------------------------------%
%-function
function [mesh] = reducemesh(mesh, ratio)

p.vertices = mesh.pnt;
p.faces = mesh.tri;
p = reducepatch(p, ratio);

mesh.pnt = p.vertices;
mesh.tri = p.faces;
%-------------------------------------%
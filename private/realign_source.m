function sourceout = realign_source(cfg, subj, source)
%REALIGN_SOURCE realign to MNI space or to freesurfer sphere
% 
% To match the location of the sources, there are three approaches (defined
% by cfg.sourcespace)
%   'volume': realign MRI to MNI with SPM/FSL, then compute forward model.
%              The dipoles are already in MNI space
%   'volume_warp': don't realign MRI to MNI, but try to compute the
%                  location of the dipoles in MNI equivalent space (it's
%                  unstable in my test)
%   'surface': use freesurfer surface and do 2-D statistics. Interpolate
%              onto high-resolution mesh in freesurfer and assign the
%              sphere coordinates. Afterwards you can interpolate and
%              average the sphere coordinates.
%
% Part of EVENTBASED/PRIVATE

switch cfg.sourcespace
  
  case 'volume_warp'
    %-------------------------------------%
    %-case 1: dipoles were warped onto MNI, so unwrap them here
    
    load(sprintf('%s/template/sourcemodel/standard_grid3d%dmm.mat', ...
      fileparts(which('ft_defaults')), cfg.bnd2lead.mni.resolution), 'grid');
    
    grid = ft_convert_units(grid, 'mm');
    sourceout = source;
    sourceout.pos = grid.pos;
    
    return
    %-------------------------------------%
    
  case 'surface'
    %-------------------------------------%
    %-case 2: project onto a sphere in freesurfer
    
    %-----------------%
    %-dir
    %-------%
    %-freesurface surfaces
    if ~isfield(cfg, 'surftype'); cfg.surftype = 'smoothwm'; end
    sdir = sprintf('%s%04d/%s', cfg.SUBJECTS_DIR, subj, 'surf/');
    %-------%
    
    %-------%
    %-low res surface (not sure if it's necessary)
    mdir = sprintf('%s%04.f/%s/%s/', cfg.data, subj, cfg.vol.mod, cfg.vol.cond); % mridata dir
    mfile = sprintf('%s_%04.f_%s_%s', cfg.rec, subj, cfg.vol.mod, cfg.vol.cond); % mridata
    gridfile = [mdir mfile '_grid'];
    %-------%
    %-----------------%
    
    %---------------------------%
    %-loop over hemisphere
    hemi = {'lh' 'rh'};
    for i = 1:numel(hemi)
      
      %-----------------%
      %-load mesh
      sphere = ft_read_headshape([sdir hemi{i} '.' 'sphere.reg']);
      load(gridfile, 'lowres', 'interpmat')
      %-----------------%
      
      %-----------------%
      %-source out (one for each hemisphere)
      [~, i_sou, i_surf] = intersect(source.pos, lowres{i}.pnt, 'rows');
      sourceout{1,i} = rmfield(source, {'pos' 'avg' 'inside'});
      sourceout{1,i}.pos = sphere.pnt;
      sourceout{1,i}.inside = 1:size(sourceout{1,i}.pos,1);
      %-----------------%
      
      %-----------------%
      %-interpolate mat
      pow = [];
      pow(i_surf,1) = source.avg.pow(i_sou);
      sourceout{1,i}.avg.pow = interpmat{i} * pow;
      %-----------------%
      
    end
    %---------------------------%
    %-------------------------------------%
    
end

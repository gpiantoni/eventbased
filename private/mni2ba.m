function [BA] = mni2ba(mni, wh_atl, Atlas)
%
%
% http://fmri.wfubmc.edu/cms/software#PickAtlas
%
% Use this reference:
% Maldjian, JA, Laurienti, PJ, Burdette, JB, Kraft RA. An Automated Method
% for Neuroanatomic and Cytoarchitectonic Atlas-based Interrogation of fMRI
% Data Sets. NeuroImage 2003. 19:1233-1239.

% 11/09/01 does not depend on warp_apply
% 10/08/30 allows for legacy use (-1)
% 10/08/26 BA based on WFU_PickAtlas
% 10/04/15 created

if nargin == 1 || wh_atl ~= -1
  % FROM WFU_PickAtlas manual 3.7
  % The segmented atlases are saved as unsigned byte or integer data ANALYZE
  % format (Mayo Clinic, Rochester USA) volumes in the MNI_atlas_templates
  % subdirectory with their corresponding lookup tables. The atlases are in
  % MNI space with dimensions of 91x109x91 sampled at 2 mm intervals,
  % corresponding to the SPM MNI templates. The atlases are in neurologic
  % convention (right of image = right of subject). Coordinates were
  % converted from Talairach space using a nonlinear transformation
  % originally described by Matthew Brett (www.mrc-cbu.cam.ac.uk/Imaging/mnispace.html)
  
  
  if nargin <= 2
    warning('off', 'MATLAB:unknownElementsNowStruc') % non informative warning (Atlas was a class, it uses as a field)
    
    load Atlas_from_WFU_PickAtlas
    AtlasName = {'brodmann', 'lobes', 'hemispheres', 'labels', 'type', 'aal', ...
      'IBASPM71', 'IBASPM116'};
    
    if nargin == 1
      wh_atl = 1:8; % all the atlas if not specified
    end
    
    for k1 = 1:size(mni,1)
      for k2 = 1:numel(wh_atl)
        BA(k1).(AtlasName{ wh_atl(k2)}) = mni2ba( mni(k1,:), wh_atl(k2), Atlas);
      end
    end
    
  else
    
    idx_mni = inv(Atlas(1).Iheader.mat) * [mni ones(1,size(mni,1))]';
    idx_mni = round(idx_mni(1:3,:))';
    idx_atl = find(Atlas( wh_atl ).Region.SubregionValues == Atlas( wh_atl ).Atlas(idx_mni(1), idx_mni(2), idx_mni(3)));
    % very strange behavior from cell, with empty index. It returns no
    % argument instead of []
    if ~isempty( idx_atl)
      BA = Atlas( wh_atl ).Region.SubregionNames{ idx_atl};
    else
      BA = '';
    end
    
    warning('on', 'MATLAB:unknownElementsNowStruc')
  end
  
else % use legacy
  % OLD VERSION: based on MRICRON
  % brod_hdr  = spm_vol('brodmann.nii'); % brodmann.nii comes from mricron templates
  % brod_area = spm_read_vols(brod_hdr);
  % save('brod_area', 'brod_area')
  % this is not very accurate, and area 48, which should be pretty small, is
  % actually pretty large.
  % Furthermore WFU_PickAtlas contains more information
  
  if nargin == 2
    load('brod_area')
    for k = 1:size(mni,1)
      BA(k) = MNI2BA( mni(k,:), -1, brod_area);
    end
  else
    
    mni  = round(mni);
    mni  = mni +[91 126 72]; % these 3 values come from brod_hdr.mat
    BA   = Atlas(mni(1), mni(2), mni(3));
  end
  
end
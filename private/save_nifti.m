function save_nifti(cfg, powsource)
%SAVE_NIFTI: save nifti with activation
% It's only a sketch. The point of this function is to create a mask for
% DTI. It should handle both volume and surface info, in the same way as
% probtrackx can handle both.

%--------%
%-prepare nifti image
if isfield(cfg.powsource, 'nifti') && ~isempty(cfg.powsource.nifti)
  
  dtimri = ft_read_mri(cfg.mriref);
  
  cfg1 = [];
  cfg1.parameter = 'image';
  souinterp = ft_sourceinterpolate(cfg1, powsource{p}, dtimri);
  
  mriname = [cfg.powsource.nifti '_' condname '_' powsource_peak(p).name];
  cfg1 = [];
  cfg1.parameter = 'image';
  cfg1.filename = mriname;
  ft_sourcewrite(cfg1, souinterp);
  gzip([mriname '.nii'])
  delete([mriname '.nii'])
end
%--------%
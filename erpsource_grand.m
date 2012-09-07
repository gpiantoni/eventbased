function erpsource_grand(cfg)
%ERPSOURCE_GRAND group-analysis of ERP source data
%
% CFG
%-Average
%  .nick: NICK to save files specific to each NICK
%  .log: name of the file and directory with analysis log
%  .subjall: index of the number of subjects
%
%  .derp: directory with ERP data
%  .erpsource.cond: cell to make averages
%
%-Statistics
%  The comparison is always against baseline. If you want to compare
%  conditions, use powstats_subj.
%
%  .erpsource.areas: how to speficy peaks to analyze, 'manual' or 'erp_peak' (peaks from granderp)
%    if 'manual'
%      .erpsource.erp_peak(1).name: string ('name_of_the_time_window')
%      .erpsource.erp_peak(1).time: scalar (center of the time window in s)
%      .erpsource.erp_peak(1).wndw: scalar (length of the time window in s)
%    if 'erp_peak'
%      .erp.refcond: string of the comparison whose peaks will be localized
%
% Options for reportsource:
%  .erpsource.clusterstatistics: 'maxsize' or 'max'
%  .erpsource.clusteralpha: level to select sensors (default 0.05)
%                           it can be a string in format '5%' to take top 5 voxels and put them in a cluster.
%  .erpsource.maxvox: max number of significant voxels to be used in soupeak
%  .erpsource.clusterthr: threshold to report clusters in output
%
%  .rslt: directory images are saved into
%
% Options if you want to create significance mask
%  .erpsource.nifti: directory and initial part of the name where you want to save the masks
%  .mriref: template for mask ('/usr/share/data/fsl-mni152-templates/MNI152_T1_1mm_brain.nii.gz')
%
% IN
%  [cfg.derp 'erpsource_SUBJ_COND'] 'erpsource_subj_A': source data for period of interest for each subject
%  [cfg.derp 'erpsource_SUBJ_COND'] 'erpsource_subj_B': source data for baseline for each subject
%
% OUT
%  [cfg.derp 'erpsource_COND'] 'erpsource': source analysis for all subject
%  [cfg.derp 'erpsource_peak_COND'] 'erpsource_peak': significant source peaks in the ERP
%
% FIGURES
%  gerp_peak_COND_ERPPEAK: 3d plot of the source for one peak
%
% Part of EVENTBASED group-analysis
% see also ERP_SUBJ, ERP_GRAND, 
% ERPSOURCE_SUBJ, ERPSOURCE_GRAND, ERPSTAT_SUBJ, ERPSTAT_GRAND,
% POW_SUBJ, POW_GRAND, POWCORR_SUBJ, POWCORR_GRAND,
% POWSOURCE_SUBJ, POWSOURCE_GRAND, POWSTAT_SUBJ, POWSTAT_GRAND,
% SOURCE_SUBJ, CONN_SUBJ, CONN_GRAND, CONN_STAT

%---------------------------%
%-start log
output = sprintf('%s began at %s on %s\n', ...
  mfilename,  datestr(now, 'HH:MM:SS'), datestr(now, 'dd-mmm-yy'));
tic_t = tic;
%---------------------------%

erp_peak = getpeak(cfg, 'erp');

%---------------------------%
%-prepare two hemisphere if surface
if strcmp(cfg.sourcespace, 'surface')
  hemi = {'lh' 'rh'};
  
  %-----------------%
  %-average sphere info
  sdir = sprintf('%s%s/%s', cfg.SUBJECTS_DIR, 'fsaverage', 'surf/');
  for h = 1:numel(hemi)
    avgsphere{h} = ft_read_headshape([sdir hemi{h} '.' 'sphere.reg']);
    avgsphere{h}.inside = true(size(avgsphere{h}.pnt,1),1);
    
    surfplot{h} = ft_read_headshape([sdir hemi{h} '.' 'pial']);
  end
  %-----------------%
  
else
  hemi = {''}; % no hemisphere
  
end
%---------------------------%

%---------------------------------------------------------%
%-statistics for main effects
%---------------------------%
%-loop over conditions
for k = 1:numel(cfg.erpsource.cond)
  cond     = cfg.erpsource.cond{k};
  condname = regexprep(cond, '*', '');
  
  %-----------------%
  %-file for each cond
  [outtmp data] = load_subj(cfg, 'erpsource', cond);
  output = [output outtmp];
  if isempty(data); continue; end
  %-----------------%
  
  %-----------------%
  %-pow or coh
  if isfield(data{1}.avg, 'coh') % coh wins
    cfg.erpsource.parameter = 'coh';
  else
    cfg.erpsource.parameter = 'pow';
  end
  %-----------------%
  
  %-----------------%
  %-loop over peaks
  erpsource_peak = [];
  erpsource = [];
  
  for p = 1:numel(erp_peak)
    output = sprintf('%s\n%s:\n', output, erp_peak(p).name);
    
    %-----------------%
    %-grand average
    cfg1 = [];
    cfg1.keepindividual = 'yes';
    cfg1.parameter = cfg.erpsource.parameter;
    gerpsouPre = ft_sourcegrandaverage(cfg1, data{:,1,p});
    gerpsource = ft_sourcegrandaverage(cfg1, data{:,2,p});
    %-----------------%
    
    %--------%
    %-do stats and figure
    h = figure;
    [soupos erpsource{p} outtmp] = reportsource(cfg.erpsource, gerpsource, gerpsouPre);
    erpsource_peak(p).pos = soupos;
    erpsource_peak(p).center = mean(soupos,1);
    erpsource_peak(p).name = erp_peak(p).name;
    output = [output outtmp];
    %--------%
    
    %--------%
    pngname = sprintf('gerp_peak_%s_%s', condname, erp_peak(p).name);
    saveas(gcf, [cfg.log filesep pngname '.png'])
    close(gcf); drawnow
    
    [~, logfile] = fileparts(cfg.log);
    system(['ln ' cfg.log filesep pngname '.png ' cfg.rslt pngname '_' logfile '.png']);
    %--------%
    
    %--------%
    %-prepare nifti image
    if isfield(cfg.erpsource, 'nifti') && ~isempty(cfg.erpsource.nifti)
      
      dtimri = ft_read_mri(cfg.mriref);
      
      cfg1 = [];
      cfg1.parameter = 'image';
      souinterp = ft_sourceinterpolate(cfg1, erpsource{p}, dtimri);
      
      mriname = [cfg.erpsource.nifti '_' condname '_' erpsource_peak(p).name];
      cfg1 = [];
      cfg1.parameter = 'image';
      cfg1.filename = mriname;
      ft_sourcewrite(cfg1, souinterp);
      gzip([mriname '.nii'])
      delete([mriname '.nii'])
    end
    %--------%
    
  end
  %-----------------%
  
  %-----------------%
  %-save
  save([cfg.derp 'erpsource_peak_' condname], 'erpsource_peak')
  
  for p = 1:numel(erpsource)
    erpsource{p}.cfg = []; % this is huge
  end
  save([cfg.derp 'erpsource_' condname], 'erpsource', '-v7.3')
  %-----------------%
  
end
%---------------------------%
%---------------------------------------------------------%

%---------------------------%
%-end log
toc_t = toc(tic_t);
outtmp = sprintf('%s ended at %s on %s after %s\n\n', ...
  mfilename, datestr(now, 'HH:MM:SS'), datestr(now, 'dd-mmm-yy'), ...
  datestr( datenum(0, 0, 0, 0, 0, toc_t), 'HH:MM:SS'));
output = [output outtmp];

%-----------------%
fprintf(output)
fid = fopen([cfg.log '.txt'], 'a');
fwrite(fid, output);
fclose(fid);
%-----------------%
%---------------------------%

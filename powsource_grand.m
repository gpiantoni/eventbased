function powsource_grand(cfg)
%POWSOURCE_GRAND group-analysis of POW source data
%
% CFG
%-Average
%  .nick: NICK to save files specific to each NICK
%  .log: name of the file and directory with analysis log
%  .subjall: index of the number of subjects
%
%  .dpow: directory with POW data
%  .powsource.cond: cell to make averages
%
%-Statistics
%  The comparison is always against baseline. If you want to compare
%  conditions, use powstats_subj.
% 
%  .powsource.areas: how to speficy peaks to analyze, 'manual' or 'powpeak'
%          (peaks from grandpow) or 'powcorrpeak' (peaks from grandpowcorr)
%    if 'manual'
%      .powsource.powpeak(1).name: string ('name_of_the_time_window')
%      .powsource.powpeak(1).time: scalar (center of the time window in s)
%      .powsource.powpeak(1).wndw: scalar (length of the time window in s)
%      .powsource.powpeak(1).freq = 10; % center of the frequency
%      .powsource.powpeak(1).band = 4; % width of the frequency band
%    if 'powpeak'
%      .pow.refcond: string of the comparison whose peaks will be localized
%    if 'powcorrpeak'
%      .powcorr.refcond: string of the comparison whose peaks will be localized
%
% Options for reportsource:
%  .powsource.clusterstatistics: 'maxsize' or 'max'
%  .powsource.clusteralpha: level to select sensors (default 0.05)
%                           it can be a string in format '5%' to take top 5 voxels and put them in a cluster.
%  .powsource.maxvox: max number of significant voxels to be used in soupeak
%  .powsource.param: 'pow' or 'coh' ('coh' only works if you specified cfg.powsource.dics.refdip)
%
%  .rslt: directory images are saved into
%
% Options if you want to create significance mask
%  .powsource.nifti: directory and initial part of the name where you want to save the masks
%  .dti.ref: template for mask ('/usr/share/data/fsl-mni152-templates/MNI152_T1_1mm_brain.nii.gz')
%
% IN
%  [cfg.dpow 'powsource_SUBJ_COND'] 'powsource_subj_A': source data for period of interest for each subject
%  [cfg.dpow 'powsource_SUBJ_COND'] 'powsource_subj_B': source data for baseline for each subject
%
% OUT
%  [cfg.dpow 'powsource_COND'] 'powsource': source analysis for all subject
%  [cfg.dpow 'powsource_peak_COND'] 'powsource_peak': significant source peaks in POW
%
% FIGURES
%  gpowpeak_COND_POWPEAK: 3d plot of the source for one peak
%
% Part of EVENTBASED group-analysis
% see also ERP_SUBJ, ERP_GRAND, 
% ERPSOURCE_SUBJ, ERPSOURCE_GRAND, ERPSTAT_SUBJ, ERPSTAT_GRAND,
% POW_SUBJ, POW_GRAND, POWCORR_SUBJ, POWCORR_GRAND,
% POWSOURCE_SUBJ, POWSOURCE_GRAND, POWSTAT_SUBJ, POWSTAT_GRAND,
% CONN_SUBJ, CONN_GRAND, CONN_STAT

%---------------------------%
%-start log
output = sprintf('%s began at %s on %s\n', ...
  mfilename,  datestr(now, 'HH:MM:SS'), datestr(now, 'dd-mmm-yy'));
tic_t = tic;
%---------------------------%

%---------------------------%
%-use predefined or power-peaks for areas of interest
if strcmp(cfg.powsource.areas, 'manual')
  powpeak = cfg.powsource.powpeak;
  
elseif strcmp(cfg.powsource.areas, 'powpeak')
  peakname = regexprep(cfg.pow.refcond, '*', '');
  load([cfg.dpow 'powpeak_' peakname], 'powpeak')
  
elseif strcmp(cfg.powsource.areas, 'powcorrpeak')
  peakname = regexprep(cfg.powcorr.refcond, '*', '');
  load([cfg.dpow 'powcorrpeak_' peakname], 'powcorrpeak')
  powpeak = powcorrpeak;
  
end
%---------------------------%

%---------------------------------------------------------%
%-statistics for main effects
%---------------------------%
%-loop over conditions
for k = 1:numel(cfg.powsource.cond)
  cond     = cfg.powsource.cond{k};
  condname = regexprep(cond, '*', '');
  
  %-----------------%
  %-file for each cond
  [outtmp data] = load_subj(cfg, 'powsource', cond);
  output = [output outtmp];
  if isempty(data); continue; end
  %-----------------%
  
  %-----------------%
  %-loop over peaks
  powsource_peak = [];
  powsource = [];
  for p = 1:numel(powpeak)
    output = sprintf('%s\n%s:\n', output, powpeak(p).name);
    
    %--------%
    %-grand average
    cfg1 = [];
    cfg1.keepindividual = 'yes';
    if ~isfield(cfg.powsource, 'param')
      cfg1.parameter = 'pow';
    else
      cfg1.parameter = cfg.powsource.param;
    end
    gpowsouPre = ft_sourcegrandaverage(cfg1, data{:,1,p});
    gpowsource = ft_sourcegrandaverage(cfg1, data{:,2,p});
    %--------%
    
    %--------%
    %-do stats and figure
    h = figure;
    [soupos powsource{p} outtmp] = reportsource(cfg.powsource, gpowsource, gpowsouPre);
    powsource_peak(p).pos = soupos;
    powsource_peak(p).center = mean(soupos,1);
    powsource_peak(p).name = powpeak(p).name;
    output = [output outtmp];
    
    %--------%
    pngname = sprintf('gpowpeak_%s_%s', condname, powpeak(p).name);
    saveas(gcf, [cfg.log filesep pngname '.png'])
    close(gcf); drawnow
    
    [~, logfile] = fileparts(cfg.log);
    system(['ln ' cfg.log filesep pngname '.png ' cfg.rslt pngname '_' logfile '.png']);
    %--------%
    
    %--------%
    %-prepare nifti image
    if isfield(cfg.powsource, 'nifti') && ~isempty(cfg.powsource.nifti)
      
      dtimri = ft_read_mri(cfg.dti.ref);
      
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
    
  end
  %-----------------%
  
  %-----------------%
  %-save
  save([cfg.dpow 'powsource_peak_' condname], 'powsource_peak')
  
  for p = 1:numel(powsource)
    powsource{p}.cfg = []; % this is huge
  end
  save([cfg.dpow 'powsource_' condname], 'powsource', '-v7.3')
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
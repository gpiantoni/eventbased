function powsource_grand(cfg)
%POWSOURCE_GRAND group-analysis of POW source data
%
% CFG
%-Average
%  .nick: NICK to save files specific to each NICK
%  .log: name of the file and directory with analysis log
%  .subjall: index of the number of subjects
%
%  .dpow: directory with ERP data
%  .powsource.cond: cell to make averages
%
%-Statistics
%  The comparison is always against baseline. If you want to compare
%  conditions, use powstats_subj.
% 
%  .powsource.areas: how to speficy peaks to analyze, 'manual' or 'powpeak'
%          (peaks from grandpow) or 'powcorrpeak' (peaks from grandpowcorr)
%    if 'manual'
%      .powsource.erppeak(1).name: string ('name_of_the_time_window')
%      .powsource.erppeak(1).time: scalar (center of the time window in s)
%      .powsource.erppeak(1).wndw: scalar (length of the time window in s)
%      .powsource.powpeak(1).freq = 10; % center of the frequency
%      .powsource.powpeak(1).band = 4; % width of the frequency band
%    if 'powpeak'
%      .pow.refcond: string of the condition whose peaks will be localized
%    if 'powcorrpeak'
%      .powcorr.refcond: string of the condition whose peaks will be localized
%
% Options for reportsource:
%  .powsource.clusterstatistics: 'maxsize' or 'max'
%  .powsource.clusteralpha: level to select sensors (default 0.05)
%                           it can be a string in format '5%' to take top 5 voxels and put them in a cluster.
%  .powsource.maxvox: max number of significant voxels to be used in soupeak
%
%  .rslt: directory images are saved into
%
% Options if you want to create significance mask
%  .powsource.nifti: directory and initial part of the name where you want to save the masks
%  .dti.ref: template for mask ('/usr/share/data/fsl-mni152-templates/MNI152_T1_1mm_brain.nii.gz')
%
% IN 
%  [cfg.dpow 'powsource_SUBJ_COND']: source data for period of interest and baseline for each subject
%
% OUT
%  [cfg.dpow 'NICK_grandpowsource']: source analysis for all subject
%  [cfg.dpow 'NICK_soupeak']: significant source peaks in POW
%
% Part of EVENTBASED group-analysis
% see also ERP_SUBJ, ERP_GRAND, ERPSOURCE_SUBJ, ERPSOURCE_GRAND,
% POW_SUBJ, POW_GRAND, POWSOURCE_SUBJ, POWSOURCE_GRAND,
% POWCORR_SUBJ, POWCORR_GRAND, POWSTAT_SUBJ, POWSTAT_GRAND,
% CONN_SUBJ, CONN_GRAND, CONN_STAT

%---------------------------%
%-start log
output = sprintf('%s began at %s on %s\n', ...
  mfilename,  datestr(now, 'HH:MM:SS'), datestr(now, 'dd-mmm-yy'));
tic_t = tic;
%---------------------------%

%---------------------------%
%-loop over conditions
for k = 1:numel(cfg.powsource.cond)
  cond     = cfg.powsource.cond{k};
  condname = regexprep(cond, '*', '');
  
  %-----------------%
  %-file for each cond
  subjfile = @(s) sprintf('%spowsource_%04d_%s.mat', cfg.dpow, s, condname);
  allname = cellfun(subjfile, num2cell(cfg.subjall), 'uni', 0);
  
  allfiles = true(1, numel(allname));
  for i = 1:numel(allname)
    if ~exist(allname{i}, 'file')
      output = [output sprintf('%s does not exist\n', allname{i})];
      allfiles(i) = false;
    end
  end
  allname = allname(allfiles);
  %-----------------%
  
  %-----------------%
  %-read data
  for s = 1:numel(allname)
    load(allname{s});
    spre(s,:) = souPre;
    sall(s,:) = source;
    clear source souPre
  end
  %-----------------%
  
  %-----------------%
  %-powsource over subj: loop over areas
  for a = 1:size(sall,2) % this is powpeak, but implicit
    cfg1 = [];
    cfg1.keepindividual = 'yes';
    gpowsouPre{k,a} = ft_sourcegrandaverage(cfg1, spre{:,a});
    gpowsource{k,a} = ft_sourcegrandaverage(cfg1, sall{:,a});
  end
  %-----------------%
  
  clear sall spre
end
%---------------------------%

%---------------------------------------------------------%
%-statistics for main effects
%-----------------%
%-use predefined or power-peaks for areas of interest
if strcmp(cfg.powsource.areas, 'manual')
  powpeak = cfg.powsource.powpeak;
  
elseif strcmp(cfg.powsource.areas, 'powpeak')
  peakname = regexprep(cfg.pow.refcond, '*', ''); % DOC: CFG.POW.REFCOND
  load([cfg.dpow cfg.cond peakname '_powpeak'], 'powpeak')
  
elseif strcmp(cfg.powsource.areas, 'powcorrpeak')
  peakname = regexprep(cfg.powcorr.refcond, '*', ''); % DOC: CFG.POWCORR.REFCOND
  load([cfg.dpow cfg.cond peakname '_powcorrpeak'], 'powcorrpeak')
  powpeak = powcorrpeak;
  
end
%-----------------%

%-------------------------------------%
%-loop over statistics conditions
for i_cond = 1:numel(cfg.powsource.cond)
  cond     = cfg.powsource.cond{i_cond};
  condname = regexprep(cond, '*', '');
  
  %-----------------%
  %-loop over peaks
  soupeak = [];
  for p = 1:numel(powpeak)
    output = sprintf('%s\n%s:\n', output, powpeak(p).name);
    
    %--------%
    %-do stats and figure
    h = figure;
    [soupos powstat{p} outtmp] = reportsource(cfg.powsource, gpowsource{i_cond, p}, gpowsouPre{i_cond, p});
    soupeak(p).pos = soupos;
    soupeak(p).center = mean(soupos,1);
    soupeak(p).name = powpeak(p).name;
    output = [output outtmp];
    %--------%
    
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
      souinterp = ft_sourceinterpolate(cfg1, powstat{p}, dtimri);
      
      cfg1 = [];
      cfg1.parameter = 'image';
      cfg1.filename = [cfg.powsource.nifti '_' condname soupeak(p).name];
      ft_sourcewrite(cfg1, souinterp);
      gzip([cfg.powsource.nifti soupeak(p).name '.nii'])
      delete([cfg.powsource.nifti soupeak(p).name '.nii'])
    end
    %--------%
    
  end
  %-----------------%
  
  %-----------------%
  %-save
  save([cfg.dpow cfg.nick '_' condname '_soupeak'], 'soupeak')
  
  for p = 1:numel(powstat)
    powstat{p}.cfg = []; % this is huge
  end
  save([cfg.dpow cfg.nick '_grandpowsource'], 'powstat', '-v7.3')
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
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
%  .erpsource.areas: how to speficy peaks to analyze, 'manual' or 'erppeak' (peaks from granderp)
%    if 'manual'
%      .erpsource.erppeak(1).name: string ('name_of_the_time_window')
%      .erpsource.erppeak(1).time: scalar (center of the time window in s)
%      .erpsource.erppeak(1).wndw: scalar (length of the time window in s)
%    if 'erppeak'
%      .erp.refcond: string of the condition whose peaks will be localized
%
% Options for reportsource:
%  .erpsource.clusterstatistics: 'maxsize' or 'max'
%  .erpsource.clusteralpha: level to select sensors (default 0.05)
%                           it can be a string in format '5%' to take top 5 voxels and put them in a cluster.
%  .erpsource.maxvox: max number of significant voxels to be used in soupeak
%
%  .rslt: directory images are saved into
%
% Options if you want to create significance mask
%  .erpsource.nifti: directory and initial part of the name where you want to save the masks
%  .dti.ref: template for mask ('/usr/share/data/fsl-mni152-templates/MNI152_T1_1mm_brain.nii.gz')
%
% IN 
%  [cfg.derp 'erpsource_SUBJ_COND']: source data for period of interest and baseline for each subject
%
% OUT
%  [cfg.derp 'NICK_granderpsource']: source analysis for all subject
%  [cfg.derp 'NICK_COND_soupeak']: significant source peaks in the ERP
%
% FIGURES
%  gerppeak_COND_ERPPEAK: 3d plot of the source for one peak
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
for k = 1:numel(cfg.erpsource.cond)
  cond     = cfg.erpsource.cond{k};
  condname = regexprep(cond, '*', '');
  
  %-----------------%
  %-file for each cond
  subjfile = @(s) sprintf('%serpsource_%04d_%s.mat', cfg.derp, s, condname);
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
  for a = 1:size(sall,2) % this is erppeak, but implicit
    cfg1 = [];
    cfg1.keepindividual = 'yes';
    cfg1.parameter = 'pow'; % instead of nai
    gerpsouPre{k,a} = ft_sourcegrandaverage(cfg1, spre{:,a});
    gerpsource{k,a} = ft_sourcegrandaverage(cfg1, sall{:,a});
  end
  %-----------------%
  
  clear sall spre
end
%---------------------------%

%---------------------------------------------------------%
%-statistics for main effects
%-----------------%
%-use predefined or erppeaks for areas of interest
if strcmp(cfg.erpsource.areas, 'manual')
  erppeak = cfg.erpsource.erppeak;
elseif strcmp(cfg.erpsource.areas, 'erppeak')
  peakname = regexprep(cfg.erp.refcond, '*', '');
  load([cfg.derp cfg.nick peakname '_erppeak'], 'erppeak')
end
%-----------------%

%-------------------------------------%
%-loop over statistics conditions
for i_cond = 1:numel(cfg.erpsource.cond)
  cond     = cfg.erpsource.cond{i_cond};
  condname = regexprep(cond, '*', '');
  
  %-----------------%
  %-loop over peaks
  soupeak = [];
  for p = 1:numel(erppeak)
    output = sprintf('%s\n%s:\n', output, erppeak(p).name);
    
    %--------%
    %-do stats and figure
    h = figure;
    [soupos erpstat{p} outtmp] = reportsource(cfg.erpsource, gerpsource{i_cond, p}, gerpsouPre{i_cond, p});
    soupeak(p).pos = soupos;
    soupeak(p).center = mean(soupos,1);
    soupeak(p).name = erppeak(p).name;
    output = [output outtmp];
    %--------%
    
    %--------%
    pngname = sprintf('gerppeak_%s_%s', condname, erppeak(p).name);
    saveas(gcf, [cfg.log filesep pngname '.png'])
    close(gcf); drawnow
    
    [~, logfile] = fileparts(cfg.log);
    system(['ln ' cfg.log filesep pngname '.png ' cfg.rslt pngname '_' logfile '.png']);
    %--------%
    
    %--------%
    %-prepare nifti image
    if isfield(cfg.erpsource, 'nifti') && ~isempty(cfg.erpsource.nifti)
      
      dtimri = ft_read_mri(cfg.dti.ref);
      
      cfg1 = [];
      cfg1.parameter = 'image';
      souinterp = ft_sourceinterpolate(cfg1, erpstat{p}, dtimri);
      
      mriname = [cfg.erpsource.nifti '_' condname '_' soupeak(p).name];
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
  save([cfg.derp cfg.nick '_' condname '_soupeak'], 'soupeak')
  
  for p = 1:numel(erpstat)
    erpstat{p}.cfg = []; % this is huge
  end
  save([cfg.derp cfg.nick '_granderpsource'], 'erpstat', '-v7.3')
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
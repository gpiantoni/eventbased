function erpstat_grand(info, opt)
%POWSTAT_GRAND group-analysis of POW source data
%
% CFG
%-Average
%  .nick: NICK to save files specific to each NICK
%  .log: name of the file and directory with analysis log
%  .subjall: index of the number of subjects
%
%  .derp: directory with POW data
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
%  .powstat.clusterstatistics: 'maxsize' or 'max'
%  .powstat.clusteralpha: level to select sensors (default 0.05)
%                           it can be a string in format '5%' to take top 5 voxels and put them in a cluster.
%  .powstat.maxvox: max number of significant voxels to be used in soupeak
%
%  .rslt: directory images are saved into
%
% Options if you want to create significance mask
%  .powstat.nifti: directory and initial part of the name where you want to save the masks
%  .dti.ref: template for mask ('/usr/share/data/fsl-mni152-templates/MNI152_T1_1mm_brain.nii.gz')
%
% IN
%  [info.dpow 'powstat_SUBJ_COND']: source data for period of interest and baseline for each subject
%
% OUT
%  [info.dpow 'powstat_COND']: source analysis for all subject
%  [info.dpow 'powsoupeak_COND']: significant source peaks in POW
%
% FIGURES
%  gpowpeak_COND_POWPEAK: 3d plot of the source for one peak
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

%---------------------------%
%-use predefined or power-peaks for areas of interest
if strcmp(cfg.powsource.areas, 'manual')
  powpeak = cfg.powsource.powpeak;
  
elseif strcmp(cfg.powsource.areas, 'powpeak')
  peakname = regexprep(cfg.pow.refcond, '*', '');
  load([info.dpow 'powpeak_' peakname], 'powpeak')
  
elseif strcmp(cfg.powsource.areas, 'powcorrpeak')
  peakname = regexprep(cfg.powcorr.refcond, '*', '');
  load([info.dpow 'powcorrpeak_' peakname], 'powcorrpeak')
  powpeak = powcorrpeak;
  
end
%---------------------------%

%-------------------------------------%
%-loop over statistics conditions
for t = 1:numel(cfg.powstat.comp) % DOC: cfg.powstat.comp
  
  %---------------------------%
  %-statistics for effects of interest
  if numel(cfg.powstat.comp{t}) == 1
    
    %-----------------%
    %-compare against baseline
    cond = cfg.powstat.comp{t}{1};
    comp = regexprep(cond, '*', '');
    output = sprintf('%s\n   COMPARISON %s\n', output, cond);
    
    %-------%
    %-pow over subj
    [outtmp data] = load_subj(info, 'powstat', cond);
    output = [output outtmp];
    if isempty(data); continue; end
    %-------%
    
    sousubj1 = squeeze(data(:,2,:)); % period of interest
    sousubj2 = squeeze(data(:,1,:)); % baseline
    %-----------------%
    
  else
    
    %-----------------%
    %-compare two conditions (cond1 - cond2)
    cond1 = cfg.gpow.comp{t}{1};
    cond2 = cfg.gpow.comp{t}{2};
    comp = [regexprep(cond1, '*', '') '_' regexprep(cond2, '*', '')];
    output = sprintf('%s\n   COMPARISON %s vs %s\n', output, cond1, cond2);
    
    %-------%
    %-pow over subj
    [outtmp data1] = load_subj(info, 'powstat', cond1);
    output = [output outtmp];
    if isempty(data); continue; end
    %-------%
    
    %-------%
    %-pow over subj
    [outtmp data2] = load_subj(info, 'powstat', cond2);
    output = [output outtmp];
    if isempty(data); continue; end
    %-------%
    
    if cfg.powstat.nai % DOC: cfg.powstat.nai
      for i1 = 1:size(data1,1)
        for i2 = 1:size(data1,2)
          data1{:,2,:}.avg.pow = data1{:,2,:}.avg.pow ./ data1{:,1,:}.avg.pow;
          data2{:,2,:}.avg.pow = data2{:,2,:}.avg.pow ./ data2{:,1,:}.avg.pow;
        end
      end
    end
    
    sousubj1 = squeeze(data1(:,2,:)); % cond1
    sousubj2 = squeeze(data2(:,2,:)); % cond2
    %-----------------%
    
  end
  clear data*
  %---------------------------%
  
  %---------------------------%
  %-loop over peaks
  soupeak = [];
  for p = 1:numel(powpeak)
    output = sprintf('%s\n%s:\n', output, powpeak(p).name);

    %-----------------%
    %-grand average
    cfg1 = [];
    cfg1.keepindividual = 'yes';
    cfg1.parameter = 'pow';
    soustat1 = ft_sourcegrandaverage(cfg1, sousubj1{:,p});
    soustat2 = ft_sourcegrandaverage(cfg1, sousubj2{:,p});
    %-----------------%    
    
    %-----------------%
    %-do stats and figure
    h = figure;
    [soupos powstat{p} outtmp] = reportsource(cfg.powstat, soustat1, soustat2); % DOC: cfg.powstat
    soupeak(p).pos = soupos;
    soupeak(p).center = mean(soupos,1);
    soupeak(p).name = powpeak(p).name;
    output = [output outtmp];
    
    %--------%
    pngname = sprintf('gpowstat_%s_%s', comp, powpeak(p).name);
    saveas(gcf, [info.log filesep pngname '.png'])
    close(gcf); drawnow
    
    [~, logfile] = fileparts(info.log);
    system(['ln ' info.log filesep pngname '.png ' info.rslt pngname '_' logfile '.png']);
    %--------%
    %-----------------%
    
  end
  %---------------------------%
  
  %-----------------%
  %-save
  save([info.dpow 'powstatpeak_' condname], 'soupeak')
  
  for p = 1:numel(powstat)
    powstat{p}.cfg = []; % this is huge
  end
  save([info.dpow 'powstat_' condname], 'powstat', '-v7.3')
  %-----------------%
  
end
%-------------------------------------%

%---------------------------%
%-end log
toc_t = toc(tic_t);
outtmp = sprintf('%s ended at %s on %s after %s\n\n', ...
  mfilename, datestr(now, 'HH:MM:SS'), datestr(now, 'dd-mmm-yy'), ...
  datestr( datenum(0, 0, 0, 0, 0, toc_t), 'HH:MM:SS'));
output = [output outtmp];

%-----------------%
fprintf(output)
fid = fopen([info.log '.txt'], 'a');
fwrite(fid, output);
fclose(fid);
%-----------------%
%---------------------------%

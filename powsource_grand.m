function powsource_grand(info, opt)
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
%  .sourcespace: 'surface' 'volume' 'volume_warp'
%  (if cfg.sourcespace == 'surface')
%  .SUBJECTS_DIR: where the Freesurfer data is stored (like the environmental variable)
%  .surfdownsample: ratio of downsampling of the surface
%
%-Statistics
%  The comparison is always against baseline. If you want to compare
%  conditions, use powstats_subj.
%
%  .powsource.peaks: how to speficy peaks to analyze, 'manual' or 'pow_peak'
%          (peaks from grandpow) or 'powcorrpeak' (peaks from grandpowcorr)
%    if 'manual'
%      .powsource.pow_peak(1).name: string ('name_of_the_time_window')
%      .powsource.pow_peak(1).time: scalar (center of the time window in s)
%      .powsource.pow_peak(1).wndw: scalar (length of the time window in s)
%      .powsource.pow_peak(1).freq = 10; % center of the frequency
%      .powsource.pow_peak(1).band = 4; % width of the frequency band
%    if 'pow_peak'
%      .pow.refcomp: string of the comparison whose peaks will be localized
%    if 'powcorrpeak'
%      .powcorr.refcomp: string of the comparison whose peaks will be localized
%
% Options for reportsource:
%  .powsource.clusterstatistics: 'maxsize' or 'max'
%  .powsource.clusteralpha: level to select sensors (default 0.05)
%                           it can be a string in format '5%' to take top 5 voxels and put them in a cluster.
%  .powsource.maxvox: max number of significant voxels to be used in soupeak
%  .powsource.clusterthr: threshold to report clusters in output
%  .powsource.atlas: index of the atlas to report labels
% 
%-Plot
%  .rslt: directory images are saved into
%
%  (if cfg.sourcespace == 'volume')
%  .template.volume: absolute path to the reference mri 
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
%  gpow_peak_COND_POWPEAK: 3d plot of the source for one peak
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

pow_peak = get_peak(cfg, 'pow');

%---------------------------%
%-prepare two hemisphere if surface
if strcmp(cfg.sourcespace, 'surface')
  hemi = {'lh' 'rh'};
  if ~isfield(cfg, 'surfdownsample'); cfg.surfdownsample = 0.01; end
  
  %-----------------%
  %-mesh for stats and plotting
  sdir = sprintf('%s%s/%s', cfg.SUBJECTS_DIR, 'fsaverage', 'surf/');
  [avgsphere, template] = read_avgsurf(sdir, cfg.surfdownsample);
  %-----------------%
  
else
  hemi = {''}; % no hemisphere
  
  %-----------------%
  %-load template
  template = ft_read_mri(cfg.template.volume);
  %-----------------%
  
end
%---------------------------%

%---------------------------------------------------------%
%-statistics for main effects
%-----------------------------------------------%
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
  %-pow or coh
  if isfield(data{1}.avg, 'coh') % coh wins
    cfg.powsource.parameter = 'coh';
  else
    cfg.powsource.parameter = 'pow';
  end
  %-----------------%
  
  %-------------------------------------%
  %-loop over peaks
  powsource_peak = [];
  powsource = [];
  for p = 1:numel(pow_peak)
    output = sprintf('%s\n%s:\n', output, pow_peak(p).name);
    
    %---------------------------%
    %-loop over hemisphere
    for h = 1:numel(hemi)
      
      output = [output hemi{h} ' '];
      
      %-----------------%
      %-interpolate to average sphere
      if strcmp(cfg.sourcespace, 'surface')
        tmpcfg = [];
        %tmpcfg.method = TODO: it's possible to specify other interpolation methods
        tmpcfg.parameter = ['avg.' cfg.powsource.parameter];
        
        for i1 = 1:size(data,1)
          for i2 = 1:size(data,2)
            data{i1,i2, p, h} = ft_sourceinterpolate(tmpcfg, data{i1,i2,p,h}, avgsphere{h});
          end
        end
      end
      %-----------------%
      
      %-----------------%
      %-grand average
      tmpcfg = [];
      tmpcfg.keepindividual = 'yes';
      tmpcfg.parameter = cfg.powsource.parameter;
      gpowsouPre = ft_sourcegrandaverage(tmpcfg, data{:,1,p,h});
      gpowsource = ft_sourcegrandaverage(tmpcfg, data{:,2,p,h});
      %-----------------%
      
      %-----------------%
      %-do stats
      cfg.powsource.channeighbstructmat = avgsphere{h}.neigh;
      [soupos powsource{p,h} outtmp] = report_source(cfg.powsource, gpowsource, gpowsouPre);
      powsource_peak(p,h).pos = soupos;
      powsource_peak(p,h).center = mean(soupos,1);
      powsource_peak(p,h).name = pow_peak(p).name;
      output = [output outtmp];
      %-----------------%
      
    end
    
    %-----------------%
    %-plot source
    h = figure;
    if strcmp(cfg.sourcespace, 'surface')
      plot_surface(powsource(p,:), template, 'stat')
      
    else
      plot_volume(powsource(p, :), template, 'stat')
      
    end
    
    %--------%
    pngname = sprintf('gpow_peak_%s_%s', condname, pow_peak(p).name);
    saveas(h, [info.log filesep pngname '.png'])
    close(h); drawnow
    
    [~, logfile] = fileparts(info.log);
    system(['ln ' info.log filesep pngname '.png ' cfg.rslt pngname '_' logfile '.png']);
    %--------%
    %-----------------%
    %---------------------------%
    
  end
  %-------------------------------------%
  
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
fid = fopen([info.log '.txt'], 'a');
fwrite(fid, output);
fclose(fid);
%-----------------%
%---------------------------%
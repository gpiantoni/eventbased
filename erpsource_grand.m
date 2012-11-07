function erpsource_grand(info, opt)
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
%  .sourcespace: 'surface' 'volume' 'volume_warp'
%  (if cfg.sourcespace == 'surface')
%  .SUBJECTS_DIR: where the Freesurfer data is stored (like the environmental variable)
%  .surfdownsample: ratio of downsampling of the surface
%
%-Statistics
%  The comparison is always against baseline. If you want to compare
%  conditions, use erpstats_subj.
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
%  .erpsource.atlas: index of the atlas to report labels
%
%  .rslt: directory images are saved into
%
% IN
%  [info.derp 'erpsource_SUBJ_COND'] 'erpsource_subj_A': source data for period of interest for each subject
%  [info.derp 'erpsource_SUBJ_COND'] 'erpsource_subj_B': source data for baseline for each subject
%
% OUT
%  [info.derp 'erpsource_COND'] 'erpsource': source analysis for all subject
%  [info.derp 'erpsource_peak_COND'] 'erpsource_peak': significant source peaks in the ERP
%
% FIGURES
%  gerp_peak_COND_ERPPEAK: 3d plot of the source for one peak
%
% Part of EVENTBASED group-analysis
% see also ERP_SUBJ, ERP_GRAND, 
% ERPSOURCE_SUBJ, ERPSOURCE_GRAND, ERPSTAT_SUBJ, ERPSTAT_GRAND,
% POW_SUBJ, POW_GRAND, POW_GRP, POWCORR_SUBJ, POWCORR_GRAND,
% POWSOURCE_SUBJ, POWSOURCE_GRAND, POWSTAT_SUBJ, POWSTAT_GRAND,
% SOURCE_SUBJ, CONN_SUBJ, CONN_GRAND, CONN_STAT

%---------------------------%
%-start log
output = sprintf('%s began at %s on %s\n', ...
  mfilename,  datestr(now, 'HH:MM:SS'), datestr(now, 'dd-mmm-yy'));
tic_t = tic;
%---------------------------%

erp_peak = get_peak(cfg, 'erp');

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
for k = 1:numel(cfg.erpsource.cond)
  cond     = cfg.erpsource.cond{k};
  condname = regexprep(cond, '*', '');
  
  %-----------------%
  %-file for each cond
  [outtmp data] = load_subj(info, 'erpsource', cond);
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
  
  %-------------------------------------%
  %-loop over peaks
  erpsource_peak = [];
  erpsource = [];
  for p = 1:numel(erp_peak)
    output = sprintf('%s\n%s:\n', output, erp_peak(p).name);
    
    %---------------------------%
    %-loop over hemisphere
    for h = 1:numel(hemi)
      
      %-----------------%
      %-interpolate to average sphere
      if strcmp(cfg.sourcespace, 'surface')
        tmpcfg = [];
        %tmpcfg.method = TODO: it's possible to specify other interpolation methods
        tmpcfg.parameter = ['avg.' cfg.erpsource.parameter];
        
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
      tmpcfg.parameter = cfg.erpsource.parameter;
      gerpsouPre = ft_sourcegrandaverage(tmpcfg, data{:,1,p,h});
      gerpsource = ft_sourcegrandaverage(tmpcfg, data{:,2,p,h});
      %-----------------%
      
      %-----------------%
      %-do stats
      cfg.erpsource.channeighbstructmat = avgsphere{h}.neigh;
      [soupos erpsource{p,h} outtmp] = report_source(cfg.erpsource, gerpsource, gerpsouPre);
      erpsource_peak(p,h).pos = soupos;
      erpsource_peak(p,h).center = mean(soupos,1);
      erpsource_peak(p,h).name = erp_peak(p).name;
      output = [output outtmp];
      %-----------------%
      
    end
    
    %-----------------%
    %-plot source
    h = figure;
    if strcmp(cfg.sourcespace, 'surface')
      plot_surface(erpsource(p,:), template, 'stat')
      
    else
      plot_volume(erpsource(p, :), template, 'stat')
      
    end
    
    %--------%
    pngname = sprintf('gerp_peak_%s_%s', condname, erp_peak(p).name);
    saveas(h, [info.log filesep pngname '.png'])
    close(h); drawnow
    
    [~, logfile] = fileparts(info.log);
    system(['ln ' info.log filesep pngname '.png ' info.rslt pngname '_' logfile '.png']);
    %--------%
    %-----------------%
    %---------------------------%
    
  end
  %-------------------------------------%
  
  %-----------------%
  %-save
  save([info.derp 'erpsource_peak_' condname], 'erpsource_peak')
  
  for p = 1:numel(erpsource)
    erpsource{p}.cfg = []; % this is huge
  end
  save([info.derp 'erpsource_' condname], 'erpsource', '-v7.3')
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

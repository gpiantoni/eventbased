function erp_grand(cfg)
%ERP_GRAND grand time lock analysis.
% 1) read single subject-data and create gerp in cfg.derp
% 2) do statistics for condition indicated by cfg.gerp.comp
% 3) plot the topoplot over time and singleplot for some electrodes
%
% CFG
%-Average
%  .nick: NICK to save files specific to each NICK
%  .log: name of the file and directory with analysis log
%  .subjall: index of the number of subjects
%
%  .derp: directory with ERP data
%  .erp.cond: conditions to make averages
%
%-Statistics
%  .gerp.comp: comparisons to test
%        (cell within cell, e.g. {{'cond1' 'cond2'} {'cond1'} {'cond2'}})
%        but you cannot have more than 2 conditions (it's always a t-test).
%   If empty, not statistics and no plots
%   If stats,
%     .gerp.stat.time: time limit for statistics (two scalars)
%     .cluster.thr: threshold to consider clusters are erppeaks
%
%-Plot
%  .gerp.bline: two scalars indicating the time window for baseline in s (only for plotting, TODO: check if necessary for normal analysis as well)
%  .gerp.chan(1).name: 'name_of_channels'
%  .gerp.chan(1).chan: cell with labels of channels of interest
%
%  .sens.layout: file with layout. It should be a mat containing 'layout'
%                If empty, it does not plot topo.
%
%  .rslt: directory images are saved into
%
% IN
%  [cfg.derp 'erp_SUBJ_COND'] 'erp_subj': timelock analysis for single-subject
%
% OUT
%  [cfg.derp 'erp_COND'] 'erp': the ERP for all subjects
%  [cfg.derp 'erp_peak_COMP'] 'erp_peak': significant peaks in the ERP for the comparison
%
% FIGURES (saved in cfg.log and, if not empty, cfg.rslt)
%  gerp_erp_COMP_CHAN: singleplot ERP, all conditions, for one channel group
%  gerp_topo_COMP: topoplot ERP for each comparison, over time
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
%-loop over conditions
for k = 1:numel(cfg.erp.cond)
  cond     = cfg.erp.cond{k};
  condname = regexprep(cond, '*', '');
  
  %-----------------%
  %-erp over subj
  [outtmp data] = load_subj(cfg, 'erp', cond);
  output = [output outtmp];
  if isempty(data); continue; end
  
  erp = ft_timelockgrandaverage([], data{:});
  erp.cfg = [];
  %-----------------%
  
  %-----------------%
  %-save
  save([cfg.derp 'erp_' condname], 'erp')
  %-----------------%
  
end
clear erp
%---------------------------%

%-----------------------------------------------%
%-stats and plots
if ~isempty(cfg.sens.layout) && ...
    ~(isfield(cfg.pow, 'source') && cfg.pow.source)
  haslay = true;
  load(cfg.sens.layout, 'layout');
  
else
  haslay = false;
  
end

if isfield(cfg.gerp, 'comp')
  
  %-------------------------------------%
  %-loop over statistics conditions
  for t = 1:numel(cfg.gerp.comp)
    
    %---------------------------%
    %-statistics for effects of interest
    clear gerp gerpall* erp_peak
    
    if numel(cfg.gerp.comp{t}) == 1
      
      %-----------------%
      %-compare against baseline
      cond = cfg.gerp.comp{t}{1};
      comp = regexprep(cond, '*', '');
      output = sprintf('%s\n   COMPARISON %s\n', output, cond);
      
      %-------%
      %-erp over subj
      [outtmp data] = load_subj(cfg, 'erp', cond);
      output = [output outtmp];
      if isempty(data); continue; end
      
      tmpcfg = [];
      gerp{1} = ft_timelockgrandaverage([], data{:});
      tmpcfg.keepindividual = 'yes';
      gerpall1 = ft_timelockgrandaverage(tmpcfg, data{:});
      %-------%
      
      %-------%
      %-data to plot
      gplot = gerp{1};
      %-------%
      
      [erp_peak stat outtmp] = report_cluster(cfg, gerpall1);
      %-----------------%
      
    else
      
      %-----------------%
      %-compare two conditions
      cond1 = cfg.gerp.comp{t}{1};
      cond2 = cfg.gerp.comp{t}{2};
      comp = [regexprep(cond1, '*', '') '_' regexprep(cond2, '*', '')];
      output = sprintf('%s\n   COMPARISON %s vs %s\n', output, cond1, cond2);
      
      %-------%
      %-erp over subj
      [outtmp, data1, data2] = load_subj(cfg, 'erp', cfg.gerp.comp{t});
      output = [output outtmp];
      if isempty(data1); continue; end
      
      tmpcfg = [];
      gerp{1} = ft_timelockgrandaverage(tmpcfg, data1{:});
      gerp{2} = ft_timelockgrandaverage(tmpcfg, data2{:});
      
      tmpcfg.keepindividual = 'yes';
      gerpall1 = ft_timelockgrandaverage(tmpcfg, data1{:});
      gerpall2 = ft_timelockgrandaverage(tmpcfg, data2{:});
      %-------%
      
      %-------%
      %-data for plot
      gplot = gerp{2};
      gplot.avg = gerp{2}.avg - gerp{1}.avg;
      %-------%
      
      [erp_peak stat outtmp] = reportcluster(cfg, gerpall1, gerpall2);
      %-----------------%
      
    end
    
    save([cfg.derp 'erp_peak_' comp], 'erp_peak', 'stat')
    output = [output outtmp];
    %---------------------------%
    
    %---------------------------%
    %-singleplotER (multiple conditions at once)
    for c = 1:numel(cfg.gerp.chan)
      
      %-----------------%
      %-figure
      h = figure;
      set(h, 'Renderer', 'painters')
      
      %--------%
      %-plot
      tmpcfg = [];
      tmpcfg.channel = cfg.gerp.chan(c).chan;
      tmpcfg.baseline = cfg.gerp.bline;
      tmpcfg.ylim = 'maxabs';
      ft_singleplotER(tmpcfg, gerp{:});
      
      legend('cond1', 'cond2')
      
      title([comp ' ' cfg.gerp.chan(c).name], 'Interpreter', 'none')
      %--------%
      %-----------------%
      
      %-----------------%
      %-save and link
      pngname = sprintf('gerp_erp_%s_%s', comp, cfg.gerp.chan(c).name);
      saveas(gcf, [cfg.log filesep pngname '.png'])
      close(gcf); drawnow
      
      [~, logfile] = fileparts(cfg.log);
      system(['ln ' cfg.log filesep pngname '.png ' cfg.rslt pngname '_' logfile '.png']);
      %-----------------%
      
    end
    %---------------------------%
    
    %---------------------------%
    %-topoplotTFR (loop over tests)
    if haslay
      
      %-----------------%
      %-figure
      h = figure;
      set(h, 'Renderer', 'painters')
      
      %--------%
      %-plot
      tmpcfg = [];
      
      timelim = gplot.time([1 end]);
      tmpcfg.xlim = timelim(1):.1:timelim(2); % one plot every 100 ms
      
      tmpcfg.zlim = [-1 1] * max(gplot.avg(:));
      
      tmpcfg.layout = layout;
      tmpcfg.style = 'straight';
      tmpcfg.marker = 'off';
      tmpcfg.comment = 'xlim';
      tmpcfg.commentpos = 'title';
      
      ft_topoplotER(tmpcfg, gplot);
      %--------%
      
      %--------%
      %-colorbar only for last subplot, to give indication on scaling
      colorbar
      %--------%
      %-----------------%
      
      %-----------------%
      %-save and link
      pngname = sprintf('gerp_topo_%s', comp);
      saveas(gcf, [cfg.log filesep pngname '.png'])
      close(gcf); drawnow
      
      [~, logfile] = fileparts(cfg.log);
      system(['ln ' cfg.log filesep pngname '.png ' cfg.rslt pngname '_' logfile '.png']);
      %-----------------%
      
    end
    %---------------------------%
    
  end
  %-------------------------------------%
  
end
%-----------------------------------------------%

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
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
%  [cfg.derp 'erp_SUBJ_COND']: timelock analysis for single-subject
%
% OUT
%  [cfg.derp 'erp_COND']: timelock analysis for all subjects
%  [cfg.derp 'erppeak_COMP']: significant peaks in the ERP for the comparison
%
% FIGURES (saved in cfg.log and, if not empty, cfg.rslt)
%  gerp_erp_COMP_CHAN: singleplot ERP, all conditions, for one channel group
%  gerp_topo_COMP: topoplot ERP for each comparison, over time
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
for k = 1:numel(cfg.erp.cond)
  cond     = cfg.erp.cond{k};
  condname = regexprep(cond, '*', '');
  
  %-----------------%
  %-erp over subj
  [outtmp data] = load_subj(cfg, 'erp', cond);
  output = [output outtmp];
  if isempty(data); continue; end
  
  cfg1 = [];
  gerp = ft_timelockgrandaverage(cfg1, data{:});
  gerp.cfg = [];
  %-----------------%
  
  %-----------------%
  %-save
  save([cfg.derp 'erp_' condname], 'gerp')
  %-----------------%
  
end
clear gerp
%---------------------------%

%-----------------------------------------------%
%-stats and plots
if ~isempty(cfg.sens.layout)
  load(cfg.sens.layout, 'layout');
end

if isfield(cfg.gerp, 'comp')
  
  %-------------------------------------%
  %-loop over statistics conditions
  for t = 1:numel(cfg.gerp.comp)
    
    %---------------------------%
    %-statistics for effects of interest
    clear gerp gerpall
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
      
      cfg1 = [];
      gerp{1} = ft_timelockgrandaverage(cfg1, data{:});
      cfg1.keepindividual = 'yes';
      gerpall1 = ft_timelockgrandaverage(cfg1, data{:});
      %-------%
      
      %-------%
      %-data to plot
      gplot = gerp{1};
      %-------%
      
      [erppeak outtmp] = reportcluster(cfg, gerpall1);
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
      
      cfg1 = [];
      gerp{1} = ft_timelockgrandaverage(cfg1, data1{:});
      gerp{2} = ft_timelockgrandaverage(cfg1, data2{:});
      
      cfg1.keepindividual = 'yes';
      gerpall1 = ft_timelockgrandaverage(cfg1, data1{:});
      gerpall2 = ft_timelockgrandaverage(cfg1, data2{:});
      %-------%
      
      %-------%
      %-data for plot
      gplot = gerp{2};
      gplot.avg = gerp{2}.avg - gerp{1}.avg;
      %-------%
      
      [erppeak outtmp] = reportcluster(cfg, gerpall1, gerpall2);
      %-----------------%
      
    end
    
    save([cfg.derp 'erppeak_' comp], 'erppeak')
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
      cfg3 = [];
      cfg3.channel = cfg.gerp.chan(c).chan;
      cfg3.baseline = cfg.gerp.bline;
      cfg3.ylim = 'maxabs';
      ft_singleplotER(cfg3, gerp{:});
      
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
    if ~isempty(cfg.sens.layout)
      
      %-----------------%
      %-figure
      h = figure;
      set(h, 'Renderer', 'painters')
      
      %--------%
      %-plot
      cfg4 = [];
      
      timelim = gplot.time([1 end]);
      cfg4.xlim = timelim(1):.1:timelim(2); % one plot every 100 ms
      
      cfg4.zlim = [-1 1] * max(gplot.avg(:));
      
      cfg4.layout = layout;
      cfg4.style = 'straight';
      cfg4.marker = 'off';
      cfg4.comment = 'xlim';
      cfg4.commentpos = 'title';
      
      ft_topoplotER(cfg4, gplot);
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
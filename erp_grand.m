function erp_grand(info, opt)
%ERP_GRAND timelock-analysis over subject
% 1) read single subject-data and create gerp in info.derp
% 2) do statistics for condition indicated by opt.comp
% 3) plot the topoplot over time and singleplot for some electrodes
%
% INFO
%  .log: name of the file and directory with analysis log
%  .derp: directory with ERP data
%  .sens.layout: file with layout. It should be a mat containing 'layout'
%                If empty, it does not plot topo.
%  .rslt: directory images are saved into
%
% CFG.OPT
%  .cond*: cell with conditions (e.g. {'*cond1' '*cond2'})'
%  .comp*: comparisons to test (cell within cell, e.g. {{'cond1' 'cond2'} {'cond1'} {'cond2'}})
%        but you cannot have more than 2 conditions (it's always a t-test). If empty, not statistics and no plots
%
%  .plot.bline: two scalars indicating the time window for baseline in s (only for plotting, TODO: check if necessary for normal analysis as well)
%  .plot.chan(1).name: 'name_of_channels'
%  .plot.chan(1).chan: cell with labels of channels of interest
%
% IN
%  [info.derp 'erp_SUBJ_COND'] 'erp_subj': timelock analysis for single-subject
%
% OUT
%  [info.derp 'erp_COND'] 'erp': the ERP for all subjects
%  [info.derp 'erp_peak_COMP'] 'erp_peak': significant peaks in the ERP for the comparison
%
% FIGURES (saved in info.log and, if not empty, info.rslt)
%  gerp_erp_COMP_CHAN: singleplot ERP, all conditions, for one channel group
%  gerp_topo_COMP: topoplot ERP for each comparison, over time
%
% * indicates obligatory parameter
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
%-by default, no plot
if ~isfield(opt, 'plot'); opt.plot = []; end
if ~isfield(opt.plot, 'chan'); opt.plot.chan = []; end
%---------------------------%

%---------------------------%
%-loop over conditions
for k = 1:numel(opt.cond)
  cond = opt.cond{k};
  condname = regexprep(cond, '*', '');
  
  %-----------------%
  %-erp over subj
  [outtmp data] = load_subj(info, 'erp', cond);
  output = [output outtmp];
  if isempty(data); continue; end
  
  erp = ft_timelockgrandaverage([], data{:});
  erp.cfg = [];
  %-----------------%
  
  %-----------------%
  %-save
  save([info.derp 'erp_' condname], 'erp')
  %-----------------%
  
end
clear erp
%---------------------------%

%-----------------------------------------------%
%-stats and plots
%---------------------------%
%-sensors
if ~isempty(info.sens.layout) %TODO: check that data is not source
  haslay = true;
  load(info.sens.layout, 'layout');
  
else
  haslay = false;
  
end

%-------%
%-sensor information to pass to report_cluster
opt.sens = info.sens;
%-------%
%---------------------------%

if isfield(opt, 'comp')
  
  %-------------------------------------%
  %-loop over statistics conditions
  for t = 1:numel(opt.comp)
    
    %---------------------------%
    %-statistics for effects of interest
    clear gerp gerpall* erp_peak
    
    if numel(opt.comp{t}) == 1
      
      %-----------------%
      %-compare against baseline
      cond = opt.comp{t}{1};
      comp = regexprep(cond, '*', '');
      output = sprintf('%s\n   COMPARISON %s\n', output, cond);
      
      %-------%
      %-erp over subj
      [outtmp data] = load_subj(info, 'erp', cond);
      output = [output outtmp];
      if isempty(data); continue; end
      
      cfg = [];
      gerp{1} = ft_timelockgrandaverage([], data{:});
      cfg.keepindividual = 'yes';
      gerpall1 = ft_timelockgrandaverage(cfg, data{:});
      %-------%
      
      %-------%
      %-data to plot
      gplot = gerp{1};
      %-------%
      
      [erp_peak stat outtmp] = report_cluster(opt, gerpall1);
      %-----------------%
      
    else
      
      %-----------------%
      %-compare two conditions
      cond1 = opt.comp{t}{1};
      cond2 = opt.comp{t}{2};
      comp = [regexprep(cond1, '*', '') '_' regexprep(cond2, '*', '')];
      output = sprintf('%s\n   COMPARISON %s vs %s\n', output, cond1, cond2);
      
      %-------%
      %-erp over subj
      [outtmp, data1, data2] = load_subj(info, 'erp', opt.comp{t});
      output = [output outtmp];
      if isempty(data1); continue; end
      
      cfg = [];
      gerp{1} = ft_timelockgrandaverage(cfg, data1{:});
      gerp{2} = ft_timelockgrandaverage(cfg, data2{:});
      
      cfg.keepindividual = 'yes';
      gerpall1 = ft_timelockgrandaverage(cfg, data1{:});
      gerpall2 = ft_timelockgrandaverage(cfg, data2{:});
      %-------%
      
      %-------%
      %-data for plot
      gplot = gerp{2};
      gplot.avg = gerp{2}.avg - gerp{1}.avg;
      %-------%
      
      [erp_peak stat outtmp] = report_cluster(opt, gerpall1, gerpall2);
      %-----------------%
      
    end
    
    save([info.derp 'erp_peak_' comp], 'erp_peak', 'stat')
    output = [output outtmp];
    %---------------------------%
    
    %---------------------------%
    %-singleplotER (multiple conditions at once)
    for c = 1:numel(opt.plot.chan)
      
      %-----------------%
      %-figure
      h = figure;
      set(h, 'Renderer', 'painters')
      
      %--------%
      %-plot
      cfg = [];
      cfg.channel = opt.plot.chan(c).chan;
      cfg.baseline = opt.plot.bline;
      cfg.ylim = 'maxabs';
      ft_singleplotER(cfg, gerp{:});
      
      legend('cond1', 'cond2')
      
      title([comp ' ' opt.plot.chan(c).name], 'Interpreter', 'none')
      %--------%
      %-----------------%
      
      %-----------------%
      %-save and link
      pngname = sprintf('gerp_erp_%s_%s', comp, opt.plot.chan(c).name);
      saveas(gcf, [info.log filesep pngname '.png'])
      close(gcf); drawnow
      
      [~, logfile] = fileparts(info.log);
      system(['ln ' info.log filesep pngname '.png ' info.rslt pngname '_' logfile '.png']);
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
      cfg = [];
      
      timelim = gplot.time([1 end]);
      cfg.xlim = timelim(1):.1:timelim(2); % one plot every 100 ms
      
      cfg.zlim = [-1 1] * max(gplot.avg(:));
      
      cfg.layout = layout;
      cfg.style = 'straight';
      cfg.marker = 'off';
      cfg.comment = 'xlim';
      cfg.commentpos = 'title';
      
      ft_topoplotER(cfg, gplot);
      %--------%
      
      %--------%
      %-colorbar only for last subplot, to give indication on scaling
      colorbar
      %--------%
      %-----------------%
      
      %-----------------%
      %-save and link
      pngname = sprintf('gerp_topo_%s', comp);
      saveas(gcf, [info.log filesep pngname '.png'])
      close(gcf); drawnow
      
      [~, logfile] = fileparts(info.log);
      system(['ln ' info.log filesep pngname '.png ' info.rslt pngname '_' logfile '.png']);
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
fid = fopen([info.log '.txt'], 'a');
fwrite(fid, output);
fclose(fid);
%-----------------%
%---------------------------%
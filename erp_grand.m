function erp_grand(cfg)
%ERP_GRAND grand time lock analysis.
% 1) read single subject-data and create gerp in cfg.derp
% 2) do statistics for condition indicated by cfg.erpeffect, to create erppeak
% 3) plot the topoplot over time and singleplot for some electrodes
%
% CFG
%  .cond: name to be used in projects/PROJNAME/subjects/0001/MOD/CONDNAME/
%  .test: a cell with the condition defined by redef.
%  .log: name of the file and directory with analysis log
%  .rslt: directory images are saved into
%  .sens.layout: file with layout. It should be a mat containing 'layout'
%                If empty, it does not plot topo.
%
%  .derp: directory to save ERP data
%  .erpeffect: effect of interest to create erppeak. If empty, no stats.
%
%  .gerp.test.time: time limit for statistics (two scalars)
%
%  .gerp.chan(1).name = 'name of channel group';
%  .gerp.chan(1).chan =  cell with labels of channels of interest
%  .gerp.bline = two scalars indicating the time window for baseline in s
%  (only for plotting, TODO: check if necessary for normal analysis as well)
%
% OUT
%  [cfg.derp 'COND_granderp']: timelock analysis for all subjects
%  [cfg.derp 'COND_erppeak']: significant peaks in the ERP
%
% FIGURES
%  gerp_erp_c01: singleplot ERP, all conditions, for one channel group
%  gerp_topo_COND: topoplot ERP for each condition, over time
%
% Part of EVENTBASED group-analysis
% see also ERP_SUBJ, ERP_GRAND, ERPSOURCE_SUBJ, ERPSOURCE_GRAND, 
% POW_SUBJ, POW_GRAND, POWSOURCE_SUBJ, POWSOURCE_GRAND, 
% POWCORR_SUBJ, POWCORR_GRAND,
% CONN_SUBJ, CONN_GRAND, CONN_STAT

%---------------------------%
%-start log
output = sprintf('%s started at %s on %s\n', ...
  mfilename,  datestr(now, 'HH:MM:SS'), datestr(now, 'dd-mmm-yy'));
tic_t = tic;
%---------------------------%

%---------------------------%
%-loop over conditions
gerp = [];
for k = 1:numel(cfg.test)
  
  %-----------------%
  %-file for each cond
  condname = regexprep(cfg.test{k}, '*', '');
  subjfile = @(s) sprintf('%serp_%02.f_%s.mat', cfg.derp, s, condname);
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
  %-erp over subj
  cfg1 = [];
  cfg1.inputfile = allname;
  gerp{k} = ft_timelockgrandaverage(cfg1);
  cfg1.keepindividual = 'yes';
  gerpall{k} = ft_timelockgrandaverage(cfg1);
  %-----------------%
  
end

%-----------------%
%-save
save([cfg.derp cfg.cond '_granderp'], 'gerp')
%-----------------%
%---------------------------%

%-----------------------------------------------%
%-stats and plots
if ~isempty(cfg.sens.layout)
  load(cfg.sens.layout, 'layout');
end

if ~isempty(gerp)
  
  %---------------------------%
  %-statistics for main effects
  for p = cfg.erpeffect
    [erppeak outtmp] = reportcluster(cfg, gerpall{p});
    
    condname = regexprep(cfg.test{p}, '*', '');
    save([cfg.derp cfg.cond condname '_erppeak'], 'erppeak')
    output = [output outtmp];
  end
  %---------------------------%
  
  %---------------------------%
  %-singleplotER (multiple conditions at once)
  for c = 1:numel(cfg.gerp.chan)
    
    %--------%
    %-figure
    h = figure;
    set(h, 'Renderer', 'painters')
    %--------%
    
    %--------%
    %-plot
    cfg3 = [];
    cfg3.channel = cfg.gerp.chan(c).chan;
    cfg3.baseline = cfg.gerp.bline;
    cfg3.ylim = 'maxabs';
    ft_singleplotER(cfg3, gerp{:});
    
    legend(cfg.test)
    
    title(cfg.gerp.chan(c).name)
    %--------%
    
    %--------%
    %-save and link
    pngname = sprintf('gerp_erp_c%02.f', c);
    saveas(gcf, [cfg.log filesep pngname '.png'])
    close(gcf); drawnow
    
    [~, logfile] = fileparts(cfg.log);
    system(['ln ' cfg.log filesep pngname '.png ' cfg.rslt pngname '_' logfile '.png']);
    %--------%
    
  end
  %---------------------------%
  
  %---------------------------%
  %-topoplotTFR (loop over tests)
  if ~isempty(cfg.sens.layout)
    for t = 1:numel(cfg.test)
      
      %--------%
      %-figure
      h = figure;
      set(h, 'Renderer', 'painters')
      %--------%
      
      %--------%
      %-plot
      
      cfg4 = [];
      %-define representative dataset to get timeinfo and max abs color
      if ~isempty(cfg.erpeffect) 
        i_gerp = cfg.erpeffect(1);
      else
        i_gerp = 1;
      end
      
      timelim = gerp{i_gerp}.time([1 end]);
      cfg4.xlim = timelim(1):.1:timelim(2); % one plot every 100 ms
  
      colorlim = max(max(abs(gerp{i_gerp}.avg)));
      cfg4.zlim = [-1 1] * colorlim;

      cfg4.layout = layout;
      cfg4.style = 'straight';
      cfg4.marker = 'off';
      cfg4.comment = 'xlim';
      cfg4.commentpos = 'title';
      
      ft_topoplotER(cfg4, gerp{t});
      %--------%
      
      %--------%
      %-save and link
      condname = regexprep(cfg.test{t}, '*', '');
      pngname = sprintf('gerp_topo_%s', condname);
      saveas(gcf, [cfg.log filesep pngname '.png'])
      close(gcf); drawnow
      
      [~, logfile] = fileparts(cfg.log);
      system(['ln ' cfg.log filesep pngname '.png ' cfg.rslt pngname '_' logfile '.png']);
      %--------%
      
    end
  end
  %---------------------------%
  
end

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
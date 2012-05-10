function erp_grand(cfg)
%ERP_GRAND grand time lock analysis.
% 1) read single subject-data and create gerp in cfg.derp
% 2) do statistics for condition indicated by cfg.gerp.stat.cond
% 3) plot the topoplot over time and singleplot for some electrodes
%
% CFG
%-Average
%  .nick: NICKNAME to save files specific to each NICKNAME
%  .log: name of the file and directory with analysis log
%  .subjall: index of the number of subjects
%
%  .derp: directory with ERP data
%  .erp.cond: conditions to make averages
%
%-Statistics
%  .gerp.stat.cond: cells within cell (e.g. {{'cond1' 'cond2'} {'cond1'} {'cond2'}})
%        but you cannot have more than 2 conditions (it's always a t-test).
%   If empty, not statistics.
%   If stats,
%     .gerp.test.time: time limit for statistics (two scalars)
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
%  [cfg.derp 'erp_SUBJCODE_CONDNAME']: timelock analysis for single-subject
%
% OUT
%  [cfg.derp 'COND_granderp']: timelock analysis for all subjects
%  [cfg.derp 'COND_erppeak']: significant peaks in the ERP
%
% FIGURES (saved in cfg.log and, if not empty, cfg.rslt)
%  gerp_erp_c01: singleplot ERP, all conditions, for one channel group
%  gerp_topo_COND: topoplot ERP for each condition, over time
%
% Part of EVENTBASED group-analysis
% see also ERP_SUBJ, ERP_GRAND, ERPSOURCE_SUBJ, ERPSOURCE_GRAND, 
% POW_SUBJ, POW_GRAND, POWSOURCE_SUBJ, POWSOURCE_GRAND, 
% POWCORR_SUBJ, POWCORR_GRAND, POWSTAT_SUBJ, POWSTAT_GRAND, 
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
for k = 1:numel(cfg.erp.cond) 
  cond     = cfg.erp.cond{k};
  condname = regexprep(cond, '*', '');
  
  %-----------------%
  %-file for each cond
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
save([cfg.derp cfg.nick '_granderp'], 'gerp')
%-----------------%
%---------------------------%

%-----------------------------------------------%
%-stats and plots
if ~isempty(cfg.sens.layout)
  load(cfg.sens.layout, 'layout');
end

if ~isempty(gerp)
  
  %-------------------------------------%
  %-loop over statistics conditions
  for t = 1:numel(cfg.gerp.stat.cond)
    
    %---------------------------%
    %-statistics for effects of interest
    if numel(cfg.gerp.stat.cond{t}) == 1
      
      %-----------------%
      %-compare against baseline
      cond = cfg.gerp.stat.cond{t}{1};
      i_cond = strfind(cfg.erp.cond, cond);
      condname = regexprep(cond, '*', '');
      
      [erppeak outtmp] = reportcluster(cfg, gerpall{i_cond});
      %-----------------%
      
      %-----------------%
      %-data for plot
      gplot = gerp{i_cond};
      gtime{1} = gplot; % to plot freq fluctuations over time
      %-----------------%
      
    else
      
      %-----------------%
      %-compare two conditions
      cond1 = cfg.gerp.stat.cond{t}{1};
      cond2 = cfg.gerp.stat.cond{t}{2};
      i_cond1 = strfind(cfg.erp.cond, cond1);
      i_cond2 = strfind(cfg.erp.cond, cond2);
      condname = [regexprep(cond1, '*', '') '_' regexprep(cond2, '*', '')];
      
      [erppeak outtmp] = reportcluster(cfg, gerpall{i_cond1}, gerpall{i_cond2});
      %-----------------%
      
      %-----------------%
      %-data for plot
      gplot = gerp{i_cond1};
      gplot.avg = log(gerp{i_cond1}.avg ./ gerp{i_cond2}.avg);
      
      gtime{1} = gerp{i_cond1}; % to plot freq fluctuations over time
      gtime{2} = gerp{i_cond2}; % to plot freq fluctuations over time
      %-----------------%
      
    end
    
    save([cfg.derp cfg.nick condname '_erppeak'], 'erppeak')
    output = [output outtmp];
  end
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
    
    title([condname ' ' cfg.gerp.chan(c).name])
    %--------%
    %-----------------%
    
    %-----------------%
    %-save and link
    pngname = sprintf('gerp_erp_%s_c%02.f', condname, c);
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
    
    cfg4.zlim = 'maxabs';
    
    cfg4.layout = layout;
    cfg4.style = 'straight';
    cfg4.marker = 'off';
    cfg4.comment = 'xlim';
    cfg4.commentpos = 'title';
    
    ft_topoplotER(cfg4, gplot);
    %--------%
    %-----------------%
    
    %-----------------%
    %-save and link
    pngname = sprintf('gerp_topo_%s', condname);
    saveas(gcf, [cfg.log filesep pngname '.png'])
    close(gcf); drawnow
    
    [~, logfile] = fileparts(cfg.log);
    system(['ln ' cfg.log filesep pngname '.png ' cfg.rslt pngname '_' logfile '.png']);
    %-----------------%
    
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
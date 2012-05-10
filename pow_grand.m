function pow_grand(cfg)
%POW_GRAND grand power average
% 1) read single subject-data and create gpow in cfg.dpow
% 2) do statistics for condition indicated by cfg.gpow.stat.cond
% 3) plot the topoplot over time, frequency and singleplot for some electrodes
%
% CFG
%-Average
%  .nick: NICKNAME to save files specific to each NICKNAME
%  .log: name of the file and directory with analysis log
%  .subjall: index of the number of subjects
%
%  .dpow: directory with ERP data
%  .pow.cond: conditions to make averages
%
%  Baseline correction at the single-subject level:
%  .pow.bl: if empty, no baseline. Otherwise:
%  .pow.bl.baseline: two scalars with baseline windows
%  .pow.bl.baselinetype: type of baseline ('relchange')
%
%  .gpow.outliers: logical (print tables with number of points above a
%  certain number of standard deviation, experimental code)
%
%-Statistics
%  .gpow.stat.cond: cells within cell (e.g. {{'cond1' 'cond2'} {'cond1'} {'cond2'}})
%        but you cannot have more than 2 conditions (it's always a t-test).
%   If empty, not statistics.
%   If stats,
%     .gpow.test.time: time limit for statistics (two scalars)
%     .gpow.test.freq: freq limit for statistics (two scalars)
%     .cluster.thr: threshold to consider clusters are erppeaks
%
%-Plot
%  .gpow.bline: two scalars indicating the time window for baseline in s (only for plotting, TODO: check if necessary for normal analysis as well)
%  .gpow.chan(1).name: 'name_of_channels'
%  .gpow.chan(1).chan: cell with labels of channels of interest
%  .gpow.freq(1).name: 'name_of_frequency'
%  .gpow.freq(1).freq: two scalars with the frequency limits
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
%  [cfg.dpow 'COND_grandpow']: power analysis for all subjects
%  [cfg.dpow 'COND_powpeak']: significant peaks in the POW
%
% FIGURES (saved in cfg.log and, if not empty, cfg.rslt)
%  gpow_tfr_c01_COND: time-frequency plot POW, for each condition, for one channel group
%  gpow_val_c01_FREQ: singleplot POW, all conditions, for one channel group, one frequency
%  gpow_topo_COND_FREQ: topoplot POW for each condition and frequency, over time
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

%-----------------------------------------------%
%-read the data
%---------------------------%
%-loop over conditions
gpow = [];
for k = 1:numel(cfg.pow.cond)
  cond     = cfg.pow.cond{k};
  condname = regexprep(cond, '*', '');
  
  %-----------------%
  %-file for each cond
  subjfile = @(s) sprintf('%spow_%02.f_%s.mat', cfg.dpow, s, condname);
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
  %-load and apply baseline correction if necessary
  freqall = [];
  for i = 1:numel(allname)
    load(allname{i}, 'freq')
    
    %-------%
    %-baseline correction
    if ~isempty(cfg.pow.bl)
      cfg3 = [];
      cfg3.baseline = cfg.pow.bl.baseline;
      cfg3.baselinetype = cfg.pow.bl.baselinetype;
      freq = ft_freqbaseline(cfg3, freq);
    end
    %-------%
    
    freqall{i} = freq;
    
  end
  %-----------------%
  
  %-----------------%
  %-read data from all subjects
  cfg1 = [];
  cfg1.keepindividual = 'yes';
  gfreq{k} = ft_freqgrandaverage(cfg1, freqall{:});
  %-----------------%
  
  %-----------------%
  %-average across subjects
  cfg2 = [];
  cfg2.variance = 'yes';
  gpow{k} = ft_freqdescriptives(cfg2, gfreq{k});
  gpow{k}.tscore =  gpow{k}.powspctrm ./ gpow{k}.powspctrmsem;
  %-----------------%
  
  %-----------------%
  %-find if there are outliers
  if isfield(cfg.gpow, 'outliers') && cfg.gpow.outliers
    
    thr = [1 2 3 5 10 20]; % threshold above sem
    
    %-------%
    %-calculate mean and standard deviation (not sem, but keep same name)
    % jackknife is more robust towards outliers
    cfg2 = [];
    cfg2.jackknife = 'yes';
    gjack = ft_freqdescriptives(cfg2, gfreq{k});
    gjack.powspctrm = shiftdim(median(gfreq{k}.powspctrm), 1);
    gjack.powspctrmsem = gjack.powspctrmsem * sqrt(numel(allname)); % sem = sd / sqrt(n)
    %-------%
    
    %-------%
    %-header
    output = [output sprintf('Subject-outliers\n')];
    output = [output 'Filename                ' sprintf('>%2.f    ', thr) sprintf('\n')];
    %-------%
    
    %-------%
    %-loop over subjects
    for s = 1:numel(allname)
      
      %-filename
      [~, subjfile] = fileparts(allname{s});
      subjfile = [subjfile repmat(' ', 1, 20-numel(subjfile))];
      
      %-matrix giving the difference in sem from the mean
      jsubj = (shiftdim(gfreq{k}.powspctrm(s,:,:,:), 1) - gjack.powspctrm) ./ gjack.powspctrmsem; % subject-specific distance from mean
      
      %-count how many points are above the threshold
      morethan = @(n)numel(find(jsubj(:) > n));
      abovethr = cellfun(morethan, num2cell(thr));
      
      %-print one row
      outtmp = sprintf('%s%s\n', subjfile, sprintf('%7.f', abovethr));
      output = [output outtmp];
      
    end
    %-------%
    
  end
  %-----------------%
  
end

%-----------------%
%-save
save([cfg.dpow cfg.cond '_grandpow'], 'gpow')
%-----------------%
%---------------------------%
%-----------------------------------------------%

%-----------------------------------------------%
%-stats and plots
if ~isempty(cfg.sens.layout)
  load(cfg.sens.layout, 'layout');
end

if ~isempty(gpow)
  
  %-------------------------------------%
  %-loop over statistics conditions
  for t = 1:numel(cfg.gpow.stat.cond)
    
    %---------------------------%
    %-statistics for effects of interest
    if numel(cfg.gpow.stat.cond{t}) == 1
      
      %-----------------%
      %-compare against baseline
      cond = cfg.gpow.stat.cond{t}{1};
      i_cond = strfind(cfg.pow.cond, cond);
      condname = regexprep(cond, '*', '');
      
      [powpeak outtmp] = reportcluster(cfg, gfreq{i_cond});
      %-----------------%
      
      %-----------------%
      %-data for plot
      gplot = gpow{i_cond};
      gtime{1} = gplot; % to plot freq fluctuations over time
      %-----------------%
      
    else
      
      %-----------------%
      %-compare two conditions
      cond1 = cfg.gpow.stat.cond{t}{1};
      cond2 = cfg.gpow.stat.cond{t}{2};
      i_cond1 = strfind(cfg.pow.cond, cond1);
      i_cond2 = strfind(cfg.pow.cond, cond2);
      condname = [regexprep(cond1, '*', '') '_' regexprep(cond2, '*', '')];
      
      [powpeak outtmp] = reportcluster(cfg, gfreq{i_cond1}, gfreq{i_cond2});
      %-----------------%
      
      %-----------------%
      %-data for plot
      gplot = gpow{i_cond1};
      gplot.powspctrm = log(gpow{i_cond1}.powspctrm ./ gpow{i_cond2}.powspctrm);
      
      gtime{1} = gpow{i_cond1}; % to plot freq fluctuations over time
      gtime{2} = gpow{i_cond2}; % to plot freq fluctuations over time
      %-----------------%
      
    end
    
    save([cfg.dpow cfg.cond condname '_powpeak'], 'powpeak')
    output = [output outtmp];
    %---------------------------%
    
    %---------------------------%
    %-loop over channels
    for c = 1:numel(cfg.gpow.chan)
      
      %-----------------%
      %-figure
      h = figure;
      set(h, 'Renderer', 'painters')
      
      %--------%
      %-options
      cfg3 = [];
      cfg3.channel = cfg.gpow.chan(c).chan;
      cfg3.zlim = 'maxabs';
      cfg3.parameter = 'powspctrm';
      ft_singleplotTFR(cfg3, gplot);
      colorbar
      
      title([condname ' ' cfg.gpow.chan(c).name])
      %--------%
      %-----------------%
      
      %-----------------%
      %-save and link
      pngname = sprintf('gpow_tfr_c%02.f_%s', c, condname);
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
      
      %-loop over freq
      for f = 1:numel(cfg.gpow.freq)
        
        %-----------------%
        %-figure
        h = figure;
        set(h, 'Renderer', 'painters')
        
        %--------%
        %-options
        cfg5 = [];
        cfg5.parameter = 'powspctrm';
        cfg5.layout = layout;
        
        cfg5.ylim = cfg.gpow.freq(f).freq;
        cfg5.zlim = 'maxabs';
        cfg5.style = 'straight';
        cfg5.marker = 'off';
        cfg5.comment = 'xlim';
        cfg5.commentpos = 'title';
        
        %-no topoplot if the data contains NaN
        onedat = squeeze(gpow{t}.powspctrm(1, cfg.gpow.freq(f).freq(1), :)); % take one example, lowest frequency
        cfg5.xlim = gpow{t}.time(~isnan(onedat));
        
        ft_topoplotER(cfg5, gplot);
        %--------%
        %-----------------%
        
        %-----------------%
        %-save and link
        pngname = sprintf('gpow_topo_%s_%s', condname, cfg.gpow.freq(f).name);
        saveas(gcf, [cfg.log filesep pngname '.png'])
        close(gcf); drawnow
        
        [~, logfile] = fileparts(cfg.log);
        system(['ln ' cfg.log filesep pngname '.png ' cfg.rslt pngname '_' logfile '.png']);
        %-----------------%
        
      end
      
    end
    %---------------------------%
    
    %---------------------------%
    %-singleplotER (multiple conditions at once)
    for c = 1:numel(cfg.gpow.chan)
      
      %-loop over freq
      for f = 1:numel(cfg.gpow.freq)
        
        %-----------------%
        %-figure
        h = figure;
        set(h, 'Renderer', 'painters')
        
        %--------%
        %-plot
        cfg4 = [];
        cfg4.channel = cfg.gpow.chan(c).chan;
        cfg4.parameter = 'powspctrm';
        cfg4.zlim = cfg.gpow.freq(f).freq;
        ft_singleplotER(cfg4, gtime{:});
        
        legend('cond1', 'cond2')
        ylabel(cfg4.parameter)
        
        title([condname ' ' cfg.gpow.chan(c).name ' ' cfg.gpow.freq(f).name])
        %--------%
        %-----------------%
        
        %-----------------%
        %-save and link
        pngname = sprintf('gpow_val_%s_%s_%s', condname, cfg.gpow.chan(c).name, cfg.gpow.freq(f).freq);
        saveas(gcf, [cfg.log filesep pngname '.png'])
        close(gcf); drawnow
        
        [~, logfile] = fileparts(cfg.log);
        system(['ln ' cfg.log filesep pngname '.png ' cfg.rslt pngname '_' logfile '.png']);
        %-----------------%
        
      end
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
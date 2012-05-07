function pow_grand(cfg)
%POW_GRAND grand power average
% 1) read single subject-data and create gpow in cfg.dpow
% 2) do statistics for condition indicated by cfg.poweffect, to create powpeak
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
%  Baseline correction at the single-subject level:
%  .pow.bl: if empty, no baseline. Otherwise:
%  .pow.bl.baseline: two scalars with baseline windows
%  .pow.bl.baselinetype: type of baseline ('relchange')
%
%  .gpow.test.time: time limit for statistics (two scalars)
%  .gpow.test.freq: freq limit for statistics (two scalars)
% 
%  .gpow.outliers: logical (print tables with number of points above a
%  certain number of standard deviation, experimental code)
%
%  .dpow: directory to save POW data
%  .poweffect: index of interest to create powpeak, can be a row vector. If empty, no stats.
%
%  .gpow.chan(1).name = 'name of channel group';
%  .gpow.chan(1).chan =  cell with labels of channels of interest
%
% OUT
%  [cfg.dpow 'COND_grandpow']: power analysis for all subjects
%  [cfg.dpow 'COND_powpeak']: significant peaks in the POW
%
% FIGURES
%  gpow_tfr_c01_COND: time-frequency plot POW, for each condition, for one channel group
%  gpow_val_c01_FREQ: singleplot POW, all conditions, for one channel group, one frequency
%  gpow_topo_COND_FREQ: topoplot POW for each condition and frequency, over time
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

%-----------------------------------------------%
%-read the data
%---------------------------%
%-loop over conditions
gpow = [];
for k = 1:numel(cfg.test)
  
  %-----------------%
  %-file for each cond
  condname = regexprep(cfg.test{k}, '*', '');
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
save([cfg.dpow cfg.proj '_grandpow'], 'gpow')
%-----------------%
%---------------------------%
%-----------------------------------------------%

%-----------------------------------------------%
%-stats and plots
if ~isempty(cfg.sens.layout)
  load(cfg.sens.layout, 'layout');
end

if ~isempty(gpow)
  
  %---------------------------%
  %-statistics for main effects
  for p = cfg.poweffect
    [powpeak outtmp] = reportcluster(cfg, gfreq{p});
    
    save([cfg.dpow cfg.cond '_powpeak'], 'powpeak')
    output = [output outtmp];
  end
  %---------------------------%
  
  %---------------------------%
  %-singleplotTFR
  for t = 1:numel(cfg.test)
    
    %-----------------%
    %-loop over channels
    for c = 1:numel(cfg.gpow.chan)
      
      %--------%
      %-figure
      h = figure;
      set(h, 'Renderer', 'painters')
      %--------%
      
      %--------%
      %-options
      cfg3 = [];
      cfg3.channel = cfg.gpow.chan(c).chan;
      cfg3.zlim = 'maxabs';
      cfg3.parameter = 'powspctrm';
      ft_singleplotTFR(cfg3, gpow{t});
      colorbar
      
      title([cfg.test{t} ' ' cfg.gpow.chan(c).name])
      %--------%
      
      %--------%
      %-save and link
      condname = regexprep(cfg.test{t}, '*', '');
      pngname = sprintf('gpow_tfr_c%02.f_%s', c, condname);
      saveas(gcf, [cfg.log filesep pngname '.png'])
      close(gcf); drawnow
      
      [~, logfile] = fileparts(cfg.log);
      system(['ln ' cfg.log filesep pngname '.png ' cfg.rslt pngname '_' logfile '.png']);
      %--------%
      
    end
    %-----------------%
    
  end
  %---------------------------%
  
  %---------------------------%
  %-singleplotER (multiple conditions at once)
  for c = 1:numel(cfg.gpow.chan)
    
    %-----------------%
    %-loop over freq
    for f = 1:numel(cfg.gpow.freq)
      
      %--------%
      %-figure
      h = figure;
      set(h, 'Renderer', 'painters')
      %--------%
      
      %--------%
      %-plot
      cfg4 = [];
      cfg4.channel = cfg.gpow.chan(c).chan;
      cfg4.parameter = 'powspctrm';
      cfg4.zlim = cfg.gpow.freq{f};
      ft_singleplotER(cfg4, gpow{:});
      
      legend(cfg.test)
      ylabel(cfg4.parameter)
      
      freqname = sprintf('f%02.f-%02.f', cfg.gpow.freq{f}(1), cfg.gpow.freq{f}(2));
      title([cfg.gpow.chan(c).name, ' ' freqname])
      %--------%
      
      %--------%
      %-save and link
      pngname = sprintf('gpow_val_c%02.f_%s', c, freqname);
      saveas(gcf, [cfg.log filesep pngname '.png'])
      close(gcf); drawnow
      
      [~, logfile] = fileparts(cfg.log);
      system(['ln ' cfg.log filesep pngname '.png ' cfg.rslt pngname '_' logfile '.png']);
      %--------%
      
    end
    %-----------------%
    
  end
  %---------------------------%
  
  %---------------------------%
  %-topoplotTFR (loop over tests)
  if ~isempty(cfg.sens.layout)
    for t = 1:numel(cfg.test)
      
      %-----------------%
      %-loop over freq
      for f = 1:numel(cfg.gpow.freq)
        
        %--------%
        %-figure
        h = figure;
        set(h, 'Renderer', 'painters')
        %--------%
        
        %--------%
        %-plot
        cfg5 = [];
        cfg5.parameter = 'powspctrm';
        cfg5.layout = layout;
        
        cfg5.ylim = cfg.gpow.freq{f};
        cfg5.zlim = 'maxabs';
        cfg5.style = 'straight';
        cfg5.marker = 'off';
        cfg5.comment = 'xlim';
        cfg5.commentpos = 'title';
        
        %-no topoplot if the data contains NaN
        onedat = squeeze(gpow{t}.powspctrm(1, cfg.gpow.freq{f}(1), :)); % take one example, lowest frequency)
        cfg5.xlim = gpow{t}.time(~isnan(onedat));
        
        ft_topoplotER(cfg5, gpow{t});
        %--------%
        
        %--------%
        %-save and link
        condname = regexprep(cfg.test{t}, '*', '');
        freqname = sprintf('f%02.f-%02.f', cfg.gpow.freq{f}(1), cfg.gpow.freq{f}(2));
        
        pngname = sprintf('gpow_topo_%s_%s', condname, freqname);
        saveas(gcf, [cfg.log filesep pngname '.png'])
        close(gcf); drawnow
        
        [~, logfile] = fileparts(cfg.log);
        system(['ln ' cfg.log filesep pngname '.png ' cfg.rslt pngname '_' logfile '.png']);
        %--------%
        
        %--------%
      end
      %-----------------%
      
    end
  end
  %---------------------------%
  
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
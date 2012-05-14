function powcorr_grand(cfg)
%POWCORR_GRAND grand power average
% 1) read single subject-data and create gpowcorr in cfg.dpow
% 2) do statistics for condition indicated by cfg.gpowcorr.comp
% 3) plot the topoplot over time, frequency and singleplot for some electrodes
%
% CFG
%-Average
%  .nick: NICK to save files specific to each NICK
%  .log: name of the file and directory with analysis log
%  .subjall: index of the number of subjects
%
%  .dpow: directory with ERP data
%  .powcorr.cond: conditions to make averages
%
%-Statistics
%  .gpowcorr.comp: cells within cell (e.g. {{'cond1' 'cond2'} {'cond1'} {'cond2'}})
%        but you cannot have more than 2 conditions (it's always a t-test).
%   If empty, not statistics.
%   If stats,
%     .gpowcorr.test.time: time limit for statistics (two scalars)
%     .gpowcorr.test.freq: freq limit for statistics (two scalars)
%     .cluster.thr: threshold to consider clusters are erppeaks
%
%-Plot
%  .gpowcorr.chan(1).name: 'name_of_channels'
%  .gpowcorr.chan(1).chan: cell with labels of channels of interest
%  .gpowcorr.freq(1).name: 'name_of_frequency'
%  .gpowcorr.freq(1).freq: two scalars with the frequency limits
%
%  .sens.layout: file with layout. It should be a mat containing 'layout'
%                If empty, it does not plot topo.
%
%  .rslt: directory images are saved into
%
% IN
%  [cfg.dpow 'powcorr_SUBJ_COND']: power correlation for single-subject
%
% OUT
%  [cfg.dpow 'powcorr_COND']: power correlation for all subjects
%  [cfg.dpow 'powcorrpeak_COMP']: significant peaks in the POWCORR for the comparison
%
% FIGURES (saved in cfg.log and, if not empty, cfg.rslt)
%  gpowcorr_tfr_COMP_COND: time-frequency plot powcorr, for each comparison, for one channel group
%  gpowcorr_val_CHAN_FREQ: singleplot powcorr, all conditions, for one channel group, one frequency
%  gpowcorr_topo_COMP_FREQ: topoplot powcorr for each condition and frequency, over time
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

%-----------------------------------------------%
%-read the data
%---------------------------%
%-loop over conditions
for k = 1:numel(cfg.powcorr.cond)
  cond     = cfg.powcorr.cond{k};
  condname = regexprep(cond, '*', '');
  
  %-----------------%
  %-powcorr over subj
  [data outtmp] = load_subj(cfg, 'powcorr', cond);
  output = [output outtmp];
  if isempty(data); continue; end
  data = data(~cellfun(@isempty, data));
  %-----------------%
  
  %-----------------%
  %-average across subjects
  cfg1 = [];
  cfg1.keepindividual = 'yes';
  gpowcorrall = ft_freqgrandaverage(cfg1, data{:});
  
  cfg2 = [];
  cfg2.variance = 'yes';
  gpowcorr = ft_freqdescriptives(cfg2, gpowcorrall);
  gpowcorr.tscore =  gpowcorr.powspctrm ./ gpowcorr.powspctrmsem;
  gpowcorr.cfg = []; % remove cfg
  %-----------------%
  
  %-----------------%
  %-save
  save([cfg.dpow 'powcorr_' condname], 'gpowcorr')
  %-----------------%
  
end
clear gpowcorr gpowcorrall
%---------------------------%

%-----------------------------------------------%
%-stats and plots
if ~isempty(cfg.sens.layout)
  load(cfg.sens.layout, 'layout');
end

if isfield(cfg.gpowcorr, 'comp')
  
  %-----------------%
  %-use gpowcorr.test info, not gpow.test info; then pass them to
  %reportcluster
  cfg.gpow.test.time = [];
  cfg.gpow.test.freq = [];
  if isfield(cfg.gpowcorr, 'test') && isfield(cfg.gpowcorr.test, 'time')
    cfg.gpow.test.time = cfg.gpowcorr.test.time;
  end
  
  if isfield(cfg.gpowcorr, 'test') && isfield(cfg.gpowcorr.test, 'freq')
    cfg.gpow.test.freq = cfg.gpowcorr.test.freq;
  end
  %-----------------%
  
  %-------------------------------------%
  %-loop over statistics conditions
  for t = 1:numel(cfg.gpowcorr.comp)
    
    %---------------------------%
    %-statistics for effects of interest
    if numel(cfg.gpowcorr.comp{t}) == 1
      
      %-----------------%
      %-compare against baseline
      cond = cfg.gpowcorr.comp{t}{1};
      comp = regexprep(cond, '*', '');
      
      %-------%
      %-powcorr over subj
      [data] = load_subj(cfg, 'powcorr', cond);
      if isempty(data); continue; end
      data = data(~cellfun(@isempty, data));
      
      cfg1 = [];
      cfg1.keepindividual = 'yes';
      gpowcorrall1 = ft_freqgrandaverage(cfg1, data{:});
      
      cfg2 = [];
      cfg2.variance = 'yes';
      gpowcorr{1} = ft_freqdescriptives(cfg2, gpowcorrall1);
      gpowcorr{1}.tscore =  gpowcorr{1}.powspctrm ./ gpowcorr{1}.powspctrmsem;
      gpowcorr{1}.cfg = []; % remove cfg
      %-------%
      
      %-------%
      %-data to plot
      gplot = gpowcorr{1};
      %-------%
      
      [powcorrpeak outtmp] = reportcluster(cfg, gpowcorrall1);
      %-----------------%
      
    else
      
      %-----------------%
      %-compare two conditions
      cond1 = cfg.gpowcorr.comp{t}{1};
      cond2 = cfg.gpowcorr.comp{t}{2};
      comp = [regexprep(cond1, '*', '') '_' regexprep(cond2, '*', '')];
      
      %-------%
      %-powcorr over subj
      [data] = load_subj(cfg, 'powcorr', cond1);
      if isempty(data); continue; end
      data = data(~cellfun(@isempty, data));
      
      cfg1 = [];
      cfg1.keepindividual = 'yes';
      gpowcorrall1 = ft_freqgrandaverage(cfg1, data{:});
      
      cfg2 = [];
      cfg2.variance = 'yes';
      gpowcorr{1} = ft_freqdescriptives(cfg2, gpowcorrall1);
      gpowcorr{1}.tscore =  gpowcorr{1}.powspctrm ./ gpowcorr{1}.powspctrmsem;
      gpowcorr{1}.cfg = []; % remove cfg
      %-------%
      
      %-------%
      %-powcorr over subj
      [data] = load_subj(cfg, 'powcorr', cond2);
      if isempty(data); continue; end
      data = data(~cellfun(@isempty, data));
      
      cfg1 = [];
      cfg1.keepindividual = 'yes';
      gpowcorrall2 = ft_freqgrandaverage(cfg1, data{:});
      
      cfg2 = [];
      cfg2.variance = 'yes';
      gpowcorr{2} = ft_freqdescriptives(cfg2, gpowcorrall2);
      gpowcorr{2}.tscore =  gpowcorr{2}.powspctrm ./ gpowcorr{2}.powspctrmsem;
      gpowcorr{2}.cfg = []; % remove cfg
      %-------%
      
      %-------%
      %-data for plot
      gplot = gpowcorr{2};
      gplot.tscore = gpowcorr{2}.tscore - gpowcorr{1}.tscore;
      %-------%
      
      [powcorrpeak outtmp] = reportcluster(cfg, gpowcorrall1, gpowcorrall2);
      %-----------------%
      
    end
    
    save([cfg.dpow 'powcorrpeak_' comp], 'powcorrpeak')
    output = [output outtmp];
    %---------------------------%
    
    %---------------------------%
    %-loop over channels
    for c = 1:numel(cfg.gpowcorr.chan)
      
      %-----------------%
      %-figure
      h = figure;
      set(h, 'Renderer', 'painters')
      
      %--------%
      %-options
      cfg3 = [];
      cfg3.channel = cfg.gpowcorr.chan(c).chan;
      cfg3.zlim = [-4 4];
      cfg3.parameter = 'tscore';
      ft_singleplotTFR(cfg3, gplot);
      colorbar
      
      title([comp ' ' cfg.gpowcorr.chan(c).name], 'Interpreter', 'none')
      %--------%
      %-----------------%
      
      %-----------------%
      %-save and link
      pngname = sprintf('gpow_tfr_%s_%s', comp, cfg.gpowcorr.chan(c).name);
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
      for f = 1:numel(cfg.gpowcorr.freq)
        
        %-----------------%
        %-figure
        h = figure;
        set(h, 'Renderer', 'painters')
        
        %--------%
        %-options
        cfg5 = [];
        cfg5.parameter = 'tscore';
        cfg5.layout = layout;
        
        cfg5.ylim = cfg.gpowcorr.freq(f).freq;
        cfg5.zlim = [-4 4];
        cfg5.style = 'straight';
        cfg5.marker = 'off';
        cfg5.comment = 'xlim';
        cfg5.commentpos = 'title';
        
        %-no topoplot if the data contains NaN
        i_freq1 = nearest(gpowcorr{f}.freq, cfg.gpowcorr.freq(f).freq(1));
        onedat = squeeze(gpowcorr{t}.tscore(1, i_freq1, :)); % take one example, lowest frequency
        cfg5.xlim = gpowcorr{t}.time(~isnan(onedat));
        
        ft_topoplotER(cfg5, gplot);
        %--------%
        %-----------------%
        
        %-----------------%
        %-save and link
        pngname = sprintf('gpowcorr_topo_%s_%s', cfg.gpowcorr.chan(c).name, cfg.gpowcorr.freq(f).name);
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
    for c = 1:numel(cfg.gpowcorr.chan)
      
      %-loop over freq
      for f = 1:numel(cfg.gpowcorr.freq)
        
        %-----------------%
        %-figure
        h = figure;
        set(h, 'Renderer', 'painters')
        
        %--------%
        %-plot
        cfg4 = [];
        cfg4.channel = cfg.gpowcorr.chan(c).chan;
        cfg4.parameter = 'tscore';
        cfg4.zlim = cfg.gpowcorr.freq(f).freq;
        ft_singleplotER(cfg4, gpowcorr{:});
        
        legend('cond1', 'cond2')
        ylabel(cfg4.parameter)
        
        title([comp ' ' cfg.gpowcorr.chan(c).name ' ' cfg.gpowcorr.freq(f).name], 'Interpreter', 'none')
        %--------%
        %-----------------%
        
        %-----------------%
        %-save and link
        pngname = sprintf('gpowcorr_val_%s_%s_%s', comp, cfg.gpowcorr.chan(c).name, cfg.gpowcorr.freq(f).name);
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
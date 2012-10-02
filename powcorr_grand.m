function powcorr_grand(info, opt)
%POWCORR_GRAND grand power average
% 1) read single subject-data and create gpowcorr in info.dpow
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
%  .powcorr.source: read virtual electrode data (logical)
%
%-Statistics
%  .gpowcorr.comp: cells within cell (e.g. {{'cond1' 'cond2'} {'cond1'} {'cond2'}})
%        but you cannot have more than 2 conditions (it's always a t-test).
%   If empty, not statistics.
%   If stats,
%     .gpowcorr.stat.time: time limit for statistics (two scalars)
%     .gpowcorr.stat.freq: freq limit for statistics (two scalars)
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
%  [info.dpow 'powcorr_SUBJ_COND'] 'powcorr_subj': power correlation for single-subject
%
% OUT
%  [info.dpow 'powcorr_COND'] 'powcorr': power correlation for all subjects
%  [info.dpow 'powcorrpeak_COMP'] 'powcorr_peak': significant peaks in the POWCORR for the comparison
%
% FIGURES (saved in info.log and, if not empty, cfg.rslt)
%  gpowcorr_tfr_COMP_COND: time-frequency plot powcorr, for each comparison, for one channel group
%  gpowcorr_val_CHAN_FREQ: singleplot powcorr, all conditions, for one channel group, one frequency
%  gpowcorr_topo_COMP_FREQ: topoplot powcorr for each condition and frequency, over time
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

%-----------------------------------------------%
%-read the data
%---------------------------%
%-loop over conditions
for k = 1:numel(cfg.powcorr.cond)
  cond     = cfg.powcorr.cond{k};
  condname = regexprep(cond, '*', '');
  
  %-----------------%
  %-powcorr over subj
  [outtmp data] = load_subj(info, 'powcorr', cond);
  output = [output outtmp];
  if isempty(data); continue; end
  %-----------------%
  
  %-----------------%
  %-average across subjects
  cfg1 = [];
  cfg1.keepindividual = 'yes';
  gpowcorrall = ft_freqgrandaverage(cfg1, data{:});
  
  cfg2 = [];
  cfg2.variance = 'yes';
  powcorr = ft_freqdescriptives(cfg2, gpowcorrall);
  powcorr.tscore =  powcorr.powspctrm ./ powcorr.powspctrmsem;
  powcorr.cfg = []; % remove cfg
  %-----------------%
  
  %-----------------%
  %-save
  save([info.dpow 'powcorr_' condname], 'powcorr')
  %-----------------%
  
end
clear gpowcorr gpowcorrall
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
      output = sprintf('%s\n   COMPARISON %s\n', output, cond);
      
      %-------%
      %-powcorr over subj
      [outtmp data] = load_subj(info, 'powcorr', cond);
      output = [output outtmp];
      if isempty(data); continue; end
      
      cfg1 = [];
      cfg1.keepindividual = 'yes';
      gpowcorrall1 = ft_freqgrandaverage(cfg1, data{:});
      
      cfg2 = [];
      cfg2.variance = 'yes';
      powcorr{1} = ft_freqdescriptives(cfg2, gpowcorrall1);
      powcorr{1}.tscore =  powcorr{1}.powspctrm ./ powcorr{1}.powspctrmsem;
      powcorr{1}.cfg = []; % remove cfg
      %-------%
      
      %-------%
      %-data to plot
      gplot = powcorr{1};
      %-------%
      
      [powcorr_peak stat outtmp] = report_cluster(cfg, gpowcorrall1);
      %-----------------%
      
    else
      
      %-----------------%
      %-compare two conditions
      cond1 = cfg.gpowcorr.comp{t}{1};
      cond2 = cfg.gpowcorr.comp{t}{2};
      comp = [regexprep(cond1, '*', '') '_' regexprep(cond2, '*', '')];
      output = sprintf('%s\n   COMPARISON %s vs %s\n', output, cond1, cond2);
      
      %-------%
      %-powcorr over subj
      [outtmp data1 data2] = load_subj(info, 'powcorr', cfg.gpowcorr.comp{t});
      output = [output outtmp];
      if isempty(data1) || isempty(data2); continue; end
      
      cfg1 = [];
      cfg1.keepindividual = 'yes';
      gpowcorrall1 = ft_freqgrandaverage(cfg1, data1{:});
      gpowcorrall2 = ft_freqgrandaverage(cfg1, data2{:});
      
      cfg2 = [];
      cfg2.variance = 'yes';
      powcorr{1} = ft_freqdescriptives(cfg2, gpowcorrall1);
      powcorr{1}.tscore =  powcorr{1}.powspctrm ./ powcorr{1}.powspctrmsem;
      powcorr{1}.cfg = []; % remove cfg

      powcorr{2} = ft_freqdescriptives(cfg2, gpowcorrall2);
      powcorr{2}.tscore =  powcorr{2}.powspctrm ./ powcorr{2}.powspctrmsem;
      powcorr{2}.cfg = []; % remove cfg
      %-------%
      
      %-------%
      %-data for plot
      gplot = powcorr{2};
      gplot.tscore = powcorr{2}.tscore - powcorr{1}.tscore;
      %-------%
      
      [powcorr_peak stat outtmp] = reportcluster(cfg, gpowcorrall1, gpowcorrall2);
      %-----------------%
      
    end
    
    save([info.dpow 'powcorr_peak_' comp], 'powcorr_peak', 'stat')
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
      saveas(gcf, [info.log filesep pngname '.png'])
      close(gcf); drawnow
      
      [~, logfile] = fileparts(info.log);
      system(['ln ' info.log filesep pngname '.png ' cfg.rslt pngname '_' logfile '.png']);
      %-----------------%
      
    end
    %---------------------------%
    
    %---------------------------%
    %-topoplotTFR (loop over tests)
    if haslay
      
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
        i_freq1 = nearest(powcorr{f}.freq, cfg.gpowcorr.freq(f).freq(1));
        onedat = squeeze(powcorr{t}.tscore(1, i_freq1, :)); % take one example, lowest frequency
        cfg5.xlim = powcorr{t}.time(~isnan(onedat));
        
        ft_topoplotER(cfg5, gplot);
        %--------%
        %-----------------%
        
        %-----------------%
        %-save and link
        pngname = sprintf('gpowcorr_topo_%s_%s', cfg.gpowcorr.chan(c).name, cfg.gpowcorr.freq(f).name);
        saveas(gcf, [info.log filesep pngname '.png'])
        close(gcf); drawnow
        
        [~, logfile] = fileparts(info.log);
        system(['ln ' info.log filesep pngname '.png ' cfg.rslt pngname '_' logfile '.png']);
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
        ft_singleplotER(cfg4, powcorr{:});
        
        legend('cond1', 'cond2')
        ylabel(cfg4.parameter)
        
        title([comp ' ' cfg.gpowcorr.chan(c).name ' ' cfg.gpowcorr.freq(f).name], 'Interpreter', 'none')
        %--------%
        %-----------------%
        
        %-----------------%
        %-save and link
        pngname = sprintf('gpowcorr_val_%s_%s_%s', comp, cfg.gpowcorr.chan(c).name, cfg.gpowcorr.freq(f).name);
        saveas(gcf, [info.log filesep pngname '.png'])
        close(gcf); drawnow
        
        [~, logfile] = fileparts(info.log);
        system(['ln ' info.log filesep pngname '.png ' cfg.rslt pngname '_' logfile '.png']);
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
fid = fopen([info.log '.txt'], 'a');
fwrite(fid, output);
fclose(fid);
%-----------------%
%---------------------------%

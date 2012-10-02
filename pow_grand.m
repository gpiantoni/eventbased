function pow_grand(info, opt)
%POW_GRAND power-analysis over subject
% 1) read single subject-data and create gpow in info.dpow
% 2) do statistics for condition indicated by opt.comp
% 3) plot the topoplot over time, frequency and singleplot for some electrodes
%
% INFO
%  .log: name of the file and directory with analysis log
%  .dpow: directory with POW data
%  .sens.layout: file with layout. It should be a mat containing 'layout'
%                If empty, it does not plot topo.
%  .rslt: directory images are saved into
%
% CFG.OPT
%  .cond*: cell with conditions (e.g. {'*cond1' '*cond2'})'
%  .comp*: comparisons to test (cell within cell, e.g. {{'cond1' 'cond2'} {'cond1'} {'cond2'}})
%        but you cannot have more than 2 conditions (it's always a t-test). If empty, not statistics and no plots
%
%  .plot.chan(1).name: 'name_of_channels'
%  .plot.chan(1).chan: cell with labels of channels of interest
%  .plot.freq(1).name: 'name_of_frequency'
%  .plot.freq(1).freq: two scalars with the frequency limits
%
% TODO: Baseline correction at the single-subject level
%
% IN
%  [info.dpow 'pow_SUBJ_COND'] 'pow_subj': timelock analysis for single-subject
%
% OUT
%  [info.dpow 'pow_COND'] 'pow': power analysis for all subjects
%  [info.dpow 'pow_peak_COMP'] 'pow_peak': significant peaks in the POW for the comparison
%
% FIGURES (saved in info.log and, if not empty, info.rslt)
%  gpow_tfr_COMP_CHAN: time-frequency plot POW, for each condition, for one channel group
%  gpow_val_CHAN_FREQ: singleplot POW, all conditions, for one channel group, one frequency % TODO: does it re-write files
%  gpow_topo_COMP_FREQ: topoplot POW for each frequency, over time
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
if ~isfield(opt.plot, 'freq'); opt.plot.freq = []; end
%---------------------------%

%-----------------------------------------------%
%-read the data
%---------------------------%
%-loop over conditions
for k = 1:numel(opt.cond)
  cond = opt.cond{k};
  condname = regexprep(cond, '*', '');
  
  %-----------------%
  %-pow over subj
  [outtmp data] = load_subj(info, 'pow', cond);
  output = [output outtmp];
  if isempty(data); continue; end
  %-----------------%
  
  %-----------------%
  %-average across subjects
  cfg = [];
  cfg.keepindividual = 'yes';
  gpowall = ft_freqgrandaverage(cfg, data{:});
  
  cfg = [];
  cfg.variance = 'yes';
  pow = ft_freqdescriptives(cfg, gpowall);
  pow.tscore =  pow.powspctrm ./ pow.powspctrmsem;
  pow.cfg = []; % remove cfg
  %-----------------%
  
  %-----------------%
  %-save
  save([info.dpow 'pow_' condname], 'pow')
  %-----------------%
  
end
clear pow powall
%---------------------------%
%-----------------------------------------------%

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
    clear pow
    
    %---------------------------%
    %-statistics for effects of interest
    if numel(opt.comp{t}) == 1
      
      %-----------------%
      %-compare against zero
      cond = opt.comp{t}{1};
      comp = regexprep(cond, '*', '');
      output = sprintf('%s\n   COMPARISON %s\n', output, cond);
      
      %-------%
      %-pow over subj
      [outtmp data] = load_subj(info, 'pow', cond);
      output = [output outtmp];
      if isempty(data); continue; end
      
      cfg = [];
      pow{1} = ft_freqgrandaverage(cfg, data{:});
      cfg.keepindividual = 'yes';
      gpowall1 = ft_freqgrandaverage(cfg, data{:});
      %-------%
      
      %-------%
      %-data to plot
      gplot = pow{1};
      %-------%

      [pow_peak stat outtmp] = report_cluster(opt, gpowall1);
      %-----------------%
      
    else
      
      %-----------------%
      %-compare two conditions
      cond1 = opt.comp{t}{1};
      cond2 = opt.comp{t}{2};
      comp = [regexprep(cond1, '*', '') '_' regexprep(cond2, '*', '')];
      output = sprintf('%s\n   COMPARISON %s vs %s\n', output, cond1, cond2);
      
      %-------%
      %-pow over subj
      [outtmp data1 data2] = load_subj(info, 'pow', opt.comp{t});
      output = [output outtmp];
      if isempty(data1) || isempty(data2); continue; end
      
      cfg = [];
      pow{1} = ft_freqgrandaverage(cfg, data1{:});
      pow{2} = ft_freqgrandaverage(cfg, data2{:});
      cfg.keepindividual = 'yes';
      gpowall1 = ft_freqgrandaverage(cfg, data1{:});
      gpowall2 = ft_freqgrandaverage(cfg, data2{:});
      %-------%
      
      %-------%
      %-data for plot
      gplot = pow{2};
      if 1
        gplot.powspctrm = log(pow{2}.powspctrm ./ pow{1}.powspctrm);
      else % TODO: with baseline correction, take the difference
        gplot.powspctrm = pow{2}.powspctrm - pow{1}.powspctrm;
      end
      %-------%
      
      [pow_peak stat outtmp] = report_cluster(opt, gpowall1, gpowall2);
      %-----------------%
      
    end
    
    save([info.dpow 'pow_peak_' comp], 'pow_peak', 'stat')
    output = [output outtmp];
    %---------------------------%
    
    %---------------------------%
    %-loop over channels
    for c = 1:numel(opt.plot.chan)
      
      %-----------------%
      %-figure
      h = figure;
      set(h, 'Renderer', 'painters')
      
      %--------%
      %-options
      cfg = [];
      cfg.channel = opt.plot.chan(c).chan;
      cfg.zlim = 'maxabs';
      cfg.parameter = 'powspctrm';
      ft_singleplotTFR(cfg, gplot);
      colorbar
      
      title([comp ' ' opt.plot.chan(c).name], 'Interpreter', 'none')
      %--------%
      %-----------------%
      
      %-----------------%
      %-save and link
      pngname = sprintf('gpow_tfr_%s_%s', comp, opt.plot.chan(c).name);
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
      
      %-loop over freq
      for f = 1:numel(opt.plot.freq)
        
        %-----------------%
        %-figure
        h = figure;
        set(h, 'Renderer', 'painters')
        
        %--------%
        %-options
        cfg = [];
        cfg.parameter = 'powspctrm';
        cfg.layout = layout;
        
        cfg.ylim = opt.plot.freq(f).freq;
        
        %-color scaling specific to each frequency band
        i_freq1 = nearest(pow{1}.freq, opt.plot.freq(f).freq(1));
        i_freq2 = nearest(pow{1}.freq, opt.plot.freq(f).freq(2));
        powspctrm = gplot.powspctrm(:, i_freq1:i_freq2, :);
        cfg.zlim = [-1 1] * max(powspctrm(:));
        
        cfg.style = 'straight';
        cfg.marker = 'off';
        cfg.comment = 'xlim';
        cfg.commentpos = 'title';
        
        %-no topoplot if the data contains NaN
        onedat = squeeze(pow{1}.powspctrm(1, i_freq1, :)); % take one example, lowest frequency
        cfg.xlim = pow{1}.time(~isnan(onedat));
        
        ft_topoplotER(cfg, gplot);
        %--------%
        
        %--------%
        colorbar
        %--------%
        %-----------------%
        
        %-----------------%
        %-save and link
        pngname = sprintf('gpow_topo_%s_%s', condname, opt.plot.freq(f).name);
        saveas(gcf, [info.log filesep pngname '.png'])
        close(gcf); drawnow
        
        [~, logfile] = fileparts(info.log);
        system(['ln ' info.log filesep pngname '.png ' info.rslt pngname '_' logfile '.png']);
        %-----------------%
        
      end
      
    end
    %---------------------------%
    
    %---------------------------%
    %-singleplotER (multiple conditions at once)
    for c = 1:numel(opt.plot.chan)
      
      %-loop over freq
      for f = 1:numel(opt.plot.freq)
        
        %-----------------%
        %-figure
        h = figure;
        set(h, 'Renderer', 'painters')
        
        %--------%
        %-plot
        cfg4 = [];
        cfg4.channel = opt.plot.chan(c).chan;
        cfg4.parameter = 'powspctrm';
        cfg4.zlim = opt.plot.freq(f).freq;
        ft_singleplotER(cfg4, pow{:});
        
        legend('cond1', 'cond2')
        ylabel(cfg4.parameter)
        
        title([comp ' ' opt.plot.chan(c).name ' ' opt.plot.freq(f).name], 'Interpreter', 'none')
        %--------%
        %-----------------%
        
        %-----------------%
        %-save and link
        pngname = sprintf('gpow_val_%s_%s_%s', comp, opt.plot.chan(c).name, opt.plot.freq(f).name);
        saveas(gcf, [info.log filesep pngname '.png'])
        close(gcf); drawnow
        
        [~, logfile] = fileparts(info.log);
        system(['ln ' info.log filesep pngname '.png ' info.rslt pngname '_' logfile '.png']);
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

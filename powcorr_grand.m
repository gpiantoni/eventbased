function powcorr_grand(info, opt)
%POWCORR_GRAND power-trial correlation over subjects
% 1) read single subject-data and create gpowcorr in info.dpow
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
%  [info.dpow 'powcorr_SUBJ_COND'] 'powcorr_subj': power correlation for single-subject
%
% OUT
%  [info.dpow 'powcorr_COND'] 'powcorr': power correlation for all subjects
%  [info.dpow 'powcorrpeak_COMP'] 'powcorr_peak': significant peaks in the POWCORR for the comparison
%
% FIGURES (saved in info.log and, if not empty, info.rslt)
%  gpowcorr_tfr_COMP_COND: time-frequency plot powcorr, for each comparison, for one channel group
%  gpowcorr_val_CHAN_FREQ: singleplot powcorr, all conditions, for one channel group, one frequency
%  gpowcorr_topo_COMP_FREQ: topoplot powcorr for each condition and frequency, over time
%
% * indicates obligatory parameter
%
% Part of EVENTBASED group-analysis
% see also ERP_SUBJ, ERP_GRAND, 
% ERPSOURCE_SUBJ, ERPSOURCE_GRAND, ERPSTAT_SUBJ, ERPSTAT_GRAND,
% POW_SUBJ, POW_GRAND, POW_GRP, POWCORR_SUBJ, POWCORR_GRAND,
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
  [outtmp data] = load_subj(info, 'powcorr', cond);
  output = [output outtmp];
  if isempty(data); continue; end
  %-----------------%
  
  %-----------------%
  %-average across subjects
  cfg = [];
  cfg.keepindividual = 'yes';
  gpowcorrall = ft_freqgrandaverage(cfg, data{:});
  
  cfg = [];
  cfg.variance = 'yes';
  powcorr = ft_freqdescriptives(cfg, gpowcorrall);
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
    clear powcorr
    
    %---------------------------%
    %-statistics for effects of interest
    if numel(opt.comp{t}) == 1
      
      %-----------------%
      %-compare against zero
      cond = opt.comp{t}{1};
      comp = regexprep(cond, '*', '');
      output = sprintf('%s\n   COMPARISON %s\n', output, cond);
      
      %-------%
      %-powcorr over subj
      [outtmp data] = load_subj(info, 'powcorr', cond);
      output = [output outtmp];
      if isempty(data); continue; end
      
      cfg = [];
      cfg.keepindividual = 'yes';
      gpowcorrall1 = ft_freqgrandaverage(cfg, data{:});
      
      cfg = [];
      cfg.variance = 'yes';
      powcorr{1} = ft_freqdescriptives(cfg, gpowcorrall1);
      powcorr{1}.tscore =  powcorr{1}.powspctrm ./ powcorr{1}.powspctrmsem;
      powcorr{1}.cfg = []; % remove cfg
      %-------%
      
      %-------%
      %-data to plot
      gplot = powcorr{1};
      %-------%
      
      [powcorr_peak stat outtmp] = report_cluster(opt, gpowcorrall1);
      %-----------------%
      
    else
      
      %-----------------%
      %-compare two conditions
      cond1 = opt.comp{t}{1};
      cond2 = opt.comp{t}{2};
      comp = [regexprep(cond1, '*', '') '_' regexprep(cond2, '*', '')];
      output = sprintf('%s\n   COMPARISON %s vs %s\n', output, cond1, cond2);
      
      %-------%
      %-powcorr over subj
      [outtmp data1 data2] = load_subj(info, 'powcorr', opt.comp{t});
      output = [output outtmp];
      if isempty(data1) || isempty(data2); continue; end
      
      cfg = [];
      cfg.keepindividual = 'yes';
      gpowcorrall1 = ft_freqgrandaverage(cfg, data1{:});
      gpowcorrall2 = ft_freqgrandaverage(cfg, data2{:});
      
      cfg = [];
      cfg.variance = 'yes';
      powcorr{1} = ft_freqdescriptives(cfg, gpowcorrall1);
      powcorr{1}.tscore =  powcorr{1}.powspctrm ./ powcorr{1}.powspctrmsem;
      powcorr{1}.cfg = []; % remove cfg

      powcorr{2} = ft_freqdescriptives(cfg, gpowcorrall2);
      powcorr{2}.tscore =  powcorr{2}.powspctrm ./ powcorr{2}.powspctrmsem;
      powcorr{2}.cfg = []; % remove cfg
      %-------%
      
      %-------%
      %-data for plot
      gplot = powcorr{2};
      gplot.tscore = powcorr{2}.tscore - powcorr{1}.tscore;
      %-------%
      
      [powcorr_peak stat outtmp] = report_cluster(opt, gpowcorrall1, gpowcorrall2);
      %-----------------%
      
    end
    
    save([info.dpow 'powcorr_peak_' comp], 'powcorr_peak', 'stat')
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
      cfg.zlim = [-4 4];
      cfg.parameter = 'tscore';
      if numel(gplot.time) > 1 % real TFR
        ft_singleplotTFR(cfg, gplot);
        colorbar
      else
        gplot.dimord = 'chan_freq'; % so it plots frequency on vertical axis
        ft_singleplotER(cfg, gplot);
        gplot.dimord = 'chan_freq_time';
      end
      
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
        cfg.parameter = 'tscore';
        cfg.layout = layout;
        
        cfg.ylim = opt.plot.freq(f).freq;
        cfg.zlim = [-4 4];
        cfg.style = 'straight';
        cfg.marker = 'off';
        cfg.comment = 'xlim';
        cfg.commentpos = 'title';
        
        %-no topoplot if the data contains NaN
        i_freq1 = nearest(powcorr{f}.freq, opt.plot.freq(f).freq(1));
        onedat = squeeze(powcorr{t}.tscore(1, i_freq1, :)); % take one example, lowest frequency
        cfg.xlim = powcorr{t}.time(~isnan(onedat));
        if numel(cfg.xlim) == 1; cfg.xlim = [1 1] * cfg.xlim; end
        
        ft_topoplotER(cfg, gplot);
        %--------%
        %-----------------%
        
        %-----------------%
        %-save and link
        pngname = sprintf('gpowcorr_topo_%s_%s', comp, opt.plot.freq(f).name);
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
        %-reassign correct dimord to pow, for the last plot
        for i = 1:numel(powcorr)
          if numel(powcorr{i}.time) == 1
            powcorr{i}.dimord = powcorr{i}.dimord(1:end-5);
          end
        end
        %-----------------%
        
        %-----------------%
        %-figure
        h = figure;
        set(h, 'Renderer', 'painters')
        
        %--------%
        %-plot
        cfg = [];
        cfg.channel = opt.plot.chan(c).chan;
        cfg.parameter = 'tscore';
        cfg.zlim = opt.plot.freq(f).freq;
        ft_singleplotER(cfg, powcorr{:});
        
        legend('cond1', 'cond2')
        ylabel(cfg.parameter)
        
        title([comp ' ' opt.plot.chan(c).name ' ' opt.plot.freq(f).name], 'Interpreter', 'none')
        %--------%
        %-----------------%
        
        %-----------------%
        %-save and link
        pngname = sprintf('gpowcorr_val_%s_%s_%s', comp, opt.plot.chan(c).name, opt.plot.freq(f).name);
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

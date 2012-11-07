function pow_grp(info, opt)
%POW_GRP opt.grp(g).nameare power between two groups of subjects
%
% INFO
%  .log: name of the file and directory with analysis log
%  .dpow: directory with POW data
%  .sens.layout: file with layout. It should be a mat containing 'layout'
%                If empty, it does not plot topo.
%  .rslt: directory images are saved into
%
% CFG.OPT
%  .grp*: structure with fields
%    .name: name of the two groups (string)
%    .grp1: index of the subjects in group 1
%    .grp2: index of the subjects in group 2
%    .cond: condition to test (a string, not a cell)
%
%  .bl: Baseline correction at the single-subject level. If empty, no baseline.
%  .bl.baseline: two scalars with baseline windows
%  .bl.baselinetype: type of baseline ('relchange')
%  .compstyle: 'logratio' or 'diff' (when comp contains two conditions,
%              take log(x2/x1) or x2 - x1)
%
%  .plot.chan(1).name: 'name_of_channels'
%  .plot.chan(1).chan: cell with labels of channels of interest
%  .plot.freq(1).name: 'name_of_frequency'
%  .plot.freq(1).freq: two scalars with the frequency limits
%
% IN
%  [info.dpow 'pow_SUBJ_COND'] 'pow_subj': power analysis for single-subject
%
% OUT
%  [info.dpow 'pow_peak_grp_GROUPNAME'] 'pow_peak': significant peaks in
%                         the POW for the opt.grp(g).namearison between the two groups
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
%-defaults
if ~isfield(opt, 'compstyle'); opt.compstyle = 'logratio'; end
if ~isfield(opt, 'plot'); opt.plot = []; end
if ~isfield(opt.plot, 'chan'); opt.plot.chan = []; end
if ~isfield(opt.plot, 'freq'); opt.plot.freq = []; end
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

%---------------------------%
%-loop over groups
for g = 1:numel(opt.grp)
  
  clear pow
  
  %-----------------%
  %-opt.grp(g).nameare two conditions
  condname = regexprep(opt.grp(g).cond, '*', '');
  output = sprintf('%s\n   COMPARISON %s between groups %s\n', output, opt.grp(g).cond, opt.grp(g).name);
  
  %-------%
  %-load the data
  tmpinfo = info;
  tmpinfo.subjall = opt.grp(g).grp1;
  [outtmp data1] = load_subj(tmpinfo, 'pow', opt.grp(g).cond);
  output = [output outtmp];
  if isempty(data1); continue; end
  
  tmpinfo = info;
  tmpinfo.subjall = opt.grp(g).grp2;
  [outtmp data2] = load_subj(tmpinfo, 'pow', opt.grp(g).cond);
  output = [output outtmp];
  if isempty(data2); continue; end
  %-------%
  
  %-------%
  %-baseline correction
  if isfield(opt, 'bl') && ~isempty(opt.bl)
    for i = 1:numel(data1)
      cfg = [];
      cfg.baseline = opt.bl.baseline;
      cfg.baselinetype = opt.bl.baselinetype;
      data1{i} = ft_freqbaseline(cfg, data1{i});
      
      if strcmp(opt.bl.baselinetype, 'relative')
        data1{i}.powspctrm = 10*log10(data1{i}.powspctrm);
      end
    end
    
    for i = 1:numel(data2)
      cfg = [];
      cfg.baseline = opt.bl.baseline;
      cfg.baselinetype = opt.bl.baselinetype;
      data2{i} = ft_freqbaseline(cfg, data2{i});
      
      if strcmp(opt.bl.baselinetype, 'relative')
        data2{i}.powspctrm = 10*log10(data2{i}.powspctrm);
      end
    end
  end
  %-------%
  
  %-------%
  %-pow over subj
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
  if strcmp(opt.compstyle, 'logratio')
    gplot.powspctrm = log(pow{2}.powspctrm ./ pow{1}.powspctrm);
  else
    gplot.powspctrm = pow{2}.powspctrm - pow{1}.powspctrm;
  end
  %-------%
  
  [pow_peak stat outtmp] = report_cluster(opt, gpowall1, gpowall2, false);
  %-----------------%
  
  save([info.dpow 'pow_peak_grp_' opt.grp(g).name], 'pow_peak', 'stat', 'gplot')
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
    if numel(gplot.time) > 1 % real TFR
      ft_singleplotTFR(cfg, gplot);
      colorbar
    else
      gplot.dimord = 'chan_freq'; % so it plots frequency on vertical axis
      ft_singleplotER(cfg, gplot);
      gplot.dimord = 'chan_freq_time';
    end
    
    title([opt.grp(g).name ' ' opt.plot.chan(c).name], 'Interpreter', 'none')
    %--------%
    %-----------------%
    
    %-----------------%
    %-save and link
    pngname = sprintf('gpow_grp_tfr_%s_%s', opt.grp(g).name, opt.plot.chan(c).name);
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
      if numel(cfg.xlim) == 1; cfg.xlim = [1 1] * cfg.xlim; end
      
      ft_topoplotER(cfg, gplot);
      %--------%
      
      %--------%
      colorbar
      %--------%
      %-----------------%
      
      %-----------------%
      %-save and link
      pngname = sprintf('gpow_grp_topo_%s_%s', opt.grp(g).name, opt.plot.freq(f).name);
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
      for i = 1:numel(pow)
        if numel(pow{i}.time) == 1
          pow{i}.dimord = pow{i}.dimord(1:end-5);
        end
      end
      %-----------------%
      
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
      
      title([opt.grp(g).name ' ' opt.plot.chan(c).name ' ' opt.plot.freq(f).name], 'Interpreter', 'none')
      %--------%
      %-----------------%
      
      %-----------------%
      %-save and link
      pngname = sprintf('gpow_grp_val_%s_%s_%s', opt.grp(g).name, opt.plot.chan(c).name, opt.plot.freq(f).name);
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
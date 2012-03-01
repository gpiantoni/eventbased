function pow_grand(cfg)
%POW_GRAND grand power average

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
  inputfile = sprintf('pow_*_%s.mat', condname);
  
  allsub = dir([cfg.dpow inputfile]);
  
  if isempty(allsub)
    outtmp = sprintf('%s does not match any file\n', condname);
    output = [output outtmp];
    continue
  end
  %-----------------%
  
  %-----------------%
  %-POW over subj
  spcell = @(name) sprintf('%s%s', cfg.dpow, name);
  allname = cellfun(spcell, {allsub.name}, 'uni', 0);
  
  cfg1 = [];
  cfg1.inputfile = allname;
  cfg1.keepindividual = 'yes';
  gfreq{k} = ft_freqgrandaverage(cfg1);
  
  cfg2 = [];
  cfg2.variance = 'yes';
  gpow{k} = ft_freqdescriptives(cfg2, gfreq{k});
  gpow{k}.tscore =  gpow{k}.powspctrm ./ gpow{k}.powspctrmsem;
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
load(cfg.sens.layout, 'layout');

if ~isempty(gpow)
  
  %---------------------------%
  %-statistics for main effects
  if ~isempty(cfg.poweffect)
    [powpeak outtmp] = reportcluster(gfreq{cfg.poweffect}, cfg);
    
    save([cfg.dpow cfg.proj '_powpeak'], 'powpeak')
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
      cfg5.xlim = gpow{t}.time;
      cfg5.ylim = cfg.gpow.freq{f};
      cfg5.zlim = 'maxabs';
      cfg5.style = 'straight';
      cfg5.marker = 'off';
      cfg5.comment = 'xlim';
      cfg5.commentpos = 'title';
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
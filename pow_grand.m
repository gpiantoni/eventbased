function pow_grand(cfg)
%POW_GRAND grand power average

mversion = 8;
%08 12/02/08 nicer plots and link to results in cfg.rslt
%07 12/02/06 deal with cases when gerp is empty
%06 12/02/03 renamed to pow_grand
%05 12/01/10 create powpeak when power is almost significant, to be used in gosdpowsource
%04 12/01/10 calculate t-statistic and sem as well
%03 12/01/10 singleplotTFR for cfg.test
%02 11/09/27 make png
%01 11/09/27 created from granderp

%---------------------------%
%-start log
output = sprintf('%s (v%02.f) started at %s on %s\n', ...
  mfilename,  mversion, datestr(now, 'HH:MM:SS'), datestr(now, 'dd-mmm-yy'));
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
  [powpeak outtmp] = reportcluster(gfreq{cfg.poweffect}, cfg);
  
  save([cfg.dpow cfg.proj '_powpeak'], 'powpeak')
  output = [output outtmp];
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
      cfg3.zlim = [-4 4];
      cfg3.parameter = 'tscore';
      ft_singleplotTFR(cfg3, gpow{t});
      colorbar
      
      title([cfg.test{t} ' ' cfg.gpow.chan(c).name])
      %--------%
      
      %--------%
      %-save and link
      condname = regexprep(cfg.test{t}, '*', '');
      pngname = sprintf('_gpow_TFR_c%02.f_%s.png', c, condname);
      saveas(gcf, [cfg.log filesep cfg.proj pngname])
      close(gcf); drawnow
      
      [~, logfile] = fileparts(cfg.log);
      system(['ln ' cfg.log filesep cfg.proj pngname ' ' cfg.rslt logfile pngname]);
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
      cfg4.parameter = 'tscore';
      cfg4.zlim = cfg.gpow.freq{f};
      ft_singleplotER(cfg4, gpow{:});
      
      legend(cfg.test)
      ylabel(cfg4.parameter)
      
      freqname = sprintf('f%02.f-%02.f', cfg.gpow.freq{f}(1), cfg.gpow.freq{f}(2));
      title([cfg.gpow.chan(c).name, ' ' freqname])
      %--------%
      
      %--------%
      %-save and link
      pngname = sprintf('_gpow_ERP_c%02.f_%s.png', c, freqname);
      saveas(gcf, [cfg.log filesep cfg.proj pngname])
      close(gcf); drawnow
      
      [~, logfile] = fileparts(cfg.log);
      system(['ln ' cfg.log filesep cfg.proj pngname ' ' cfg.rslt logfile pngname]);
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
      cfg5.parameter = 'tscore';
      cfg5.layout = layout;
      cfg5.xlim = gpow{t}.time;
      cfg5.ylim = cfg.gpow.freq{f};
      cfg5.zlim = [-4 4];
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

      pngname = sprintf('_gpow_topo_%s_%s.png', condname, freqname);
      saveas(gcf, [cfg.log filesep cfg.proj pngname])
      close(gcf); drawnow
      
      [~, logfile] = fileparts(cfg.log);
      system(['ln ' cfg.log filesep cfg.proj pngname ' ' cfg.rslt logfile pngname]);
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
outtmp = sprintf('%s (v%02.f) ended at %s on %s after %s\n\n', ...
  mfilename, mversion, datestr(now, 'HH:MM:SS'), datestr(now, 'dd-mmm-yy'), ...
  datestr( datenum(0, 0, 0, 0, 0, toc_t), 'HH:MM:SS'));
output = [output outtmp];

%-----------------%
fprintf(output)
fid = fopen([cfg.log '.txt'], 'a');
fwrite(fid, output);
fclose(fid);
%-----------------%
%---------------------------%
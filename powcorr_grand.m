function powcorr_grand(cfg)
%POWCORR_GRAND grand power average for powcorr
% 1) read single subject-data and create gpowcorr in cfg.dpow
% 2) do statistics for condition indicated by cfg.powcorreffect, to create powcorrpeak
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
%  .dpow: directory to save POWCORR data
%  .powcorreffect: effect of interest to create powcorrpeak. If empty, no stats.
%
%  .gpowcorr.chan(1).name = 'name of channel group';
%  .gpowcorr.chan(1).chan =  cell with labels of channels of interest
%
% OUT
%  [cfg.dpow 'COND_grandpowcorr']: power analysis for all subjects
%  [cfg.dpow 'COND_powcorrpeak']: significant peaks in the POWCORR
%
% FIGURES
%  gpowcorr_tfr_c01_COND: time-frequency plot POWCORR, for each condition, for one channel group
%  gpowcorr_val_c01_FREQ: singleplot POWCORR, all conditions, for one channel group, one frequency
%  gpowcorr_topo_COND_FREQ: topoplot POWCORR for each condition and frequency, over time
%
% Part of EVENTBASED group-analysis
% see also ERP_SUBJ, ERP_GRAND, ERPSOURCE_SUBJ, ERPSOURCE_GRAND,
% POW_SUBJ, POW_GRAND, POWSOURCE_SUBJ, POWSOURCE_GRAND,
% POWCORR_SUBJ, POWCORR_GRAND,
% CONN_SUBJ, CONN_GRAND, CONN_STAT

% practically identical to pow_grand but not powpeaks
% pow -> powcorr; cfg.gpow -> cfg.gpowcorr

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
  subjfile = @(s) sprintf('%spowcorr_%02.f_%s.mat', cfg.dpow, s, condname);
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
  %-POWCORR over subj
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
save([cfg.dpow cfg.proj '_grandpowcorr'], 'gpow')
%-----------------%
%---------------------------%

%-----------------------------------------------%
%-stats and plots
if ~isempty(cfg.sens.layout)
  load(cfg.sens.layout, 'layout');
end

if ~isempty(gpow)
  
  %---------------------------%
  %-statistics for main effects
  if ~isempty(cfg.poweffect)
    [powcorrpeak outtmp] = reportcluster(gfreq{cfg.powcorreffect}, cfg);
    outtmp = sprintf('%sattention! .wndw is adjusted incorrectly based on cfg.pow, not cfg.powcorr\n', output);
    
    save([cfg.dpow cfg.proj '_powcorrpeak'], 'powcorrpeak')
    output = [output outtmp];
  end
  %---------------------------%
  
  %---------------------------%
  %-singleplotTFR
  for t = 1:numel(cfg.test)
    
    %-----------------%
    %-loop over channels
    for c = 1:numel(cfg.gpowcorr.chan)
      
      %--------%
      %-figure
      h = figure;
      set(h, 'Renderer', 'painters')
      %--------%
      
      %--------%
      %-options
      cfg3 = [];
      cfg3.channel = cfg.gpowcorr.chan(c).chan;
      cfg3.zlim = 'maxabs';
      cfg3.parameter = 'tscore';
      ft_singleplotTFR(cfg3, gpow{t});
      colorbar
      
      title([cfg.test{t} ' ' cfg.gpowcorr.chan(c).name])
      %--------%
      
      %--------%
      %-save and link
      condname = regexprep(cfg.test{t}, '*', '');
      pngname = sprintf('gpowcorr_tfr_c%02.f_%s', c, condname);
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
  for c = 1:numel(cfg.gpowcorr.chan)
    
    %-----------------%
    %-loop over freq
    for f = 1:numel(cfg.gpowcorr.freq)
      
      %--------%
      %-figure
      h = figure;
      set(h, 'Renderer', 'painters')
      %--------%
      
      %--------%
      %-plot
      cfg4 = [];
      cfg4.channel = cfg.gpowcorr.chan(c).chan;
      cfg4.parameter = 'tscore';
      cfg4.zlim = cfg.gpowcorr.freq{f};
      ft_singleplotER(cfg4, gpow{:});
      
      legend(cfg.test)
      ylabel(cfg4.parameter)
      
      freqname = sprintf('f%02.f-%02.f', cfg.gpowcorr.freq{f}(1), cfg.gpowcorr.freq{f}(2));
      title([cfg.gpowcorr.chan(c).name, ' ' freqname])
      %--------%
      
      %--------%
      %-save and link
      pngname = sprintf('gpowcorr_val_c%02.f_%s', c, freqname);
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
      for f = 1:numel(cfg.gpowcorr.freq)
        
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
        cfg5.ylim = cfg.gpowcorr.freq{f};
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
        freqname = sprintf('f%02.f-%02.f', cfg.gpowcorr.freq{f}(1), cfg.gpowcorr.freq{f}(2));
        
        pngname = sprintf('gpowcorr_topo_%s_%s', condname, freqname);
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
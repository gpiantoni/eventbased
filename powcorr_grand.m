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
%  .powcorreffect: index of interest to create powcorrpeak, can be a row vector. If empty, no stats.
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

% practically identical to pow_grand
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
gpowcorr = [];
for k = 1:numel(cfg.powcorr.cond) % DOC: CFG.POWCORR.COND
  cond     = cfg.powcorr.cond{k};
  condname = regexprep(cond, '*', '');
  
  %-----------------%
  %-file for each cond
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
  gpowcorr{k} = ft_freqdescriptives(cfg2, gfreq{k});
  gpowcorr{k}.tscore =  gpowcorr{k}.powspctrm ./ gpowcorr{k}.powspctrmsem;
  %-----------------%
  
end

%-----------------%
%-save
save([cfg.dpow cfg.cond '_grandpowcorr'], 'gpowcorr')
%-----------------%
%---------------------------%

%-----------------------------------------------%
%-stats and plots
if ~isempty(cfg.sens.layout)
  load(cfg.sens.layout, 'layout');
end

if ~isempty(gpowcorr)
  
  %-------------------------------------%
  %-loop over statistics conditions
  for t = 1:numel(cfg.gpowcorr.stat.cond) % DOC: cfg.gpowcorr.stat.cond as cell
    
    %---------------------------%
    %-statistics for effects of interest
    if numel(cfg.gpowcorr.stat.cond{t}) == 1
      
      %-----------------%
      %-compare against baseline
      cond = cfg.gpowcorr.stat.cond{t}{1};
      i_cond = strfind(cfg.powcorr.cond, cond);
      condname = regexprep(cond, '*', '');
      
      [powcorrpeak outtmp] = reportcluster(cfg, gfreq{i_cond});
      %-----------------%
      
      %-----------------%
      %-data for plot
      gplot = gpowcorr{i_cond};
      gtime{1} = gplot; % to plot freq fluctuations over time
      %-----------------%
      
    else
      
      %-----------------%
      %-compare two conditions
      cond1 = cfg.gpowcorr.stat.cond{t}{1};
      cond2 = cfg.gpowcorr.stat.cond{t}{2};
      i_cond1 = strfind(cfg.powcorr.cond, cond1);
      i_cond2 = strfind(cfg.powcorr.cond, cond2);
      condname = [regexprep(cond1, '*', '') '_' regexprep(cond2, '*', '')];
      
      [powcorrpeak outtmp] = reportcluster(cfg, gfreq{i_cond1}, gfreq{i_cond2});
      %-----------------%
      
      %-----------------%
      %-data for plot
      gplot = gpowcorr{i_cond1};
      gplot.tscore = log(gpowcorr{i_cond1}.tscore ./ gpowcorr{i_cond2}.tscore);
      
      gtime{1} = gpowcorr{i_cond1}; % to plot freq fluctuations over time
      gtime{2} = gpowcorr{i_cond2}; % to plot freq fluctuations over time
      %-----------------%
      
    end
    
    save([cfg.dpow cfg.cond condname '_powcorrpeak'], 'powcorrpeak')
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
      
      title([condname ' ' cfg.gpowcorr.chan(c).name])
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
        
        cfg5.ylim = cfg.gpowcorr.freq{f};
        cfg5.zlim = [-4 4];
        cfg5.style = 'straight';
        cfg5.marker = 'off';
        cfg5.comment = 'xlim';
        cfg5.commentpos = 'title';
        
        %-no topoplot if the data contains NaN
        onedat = squeeze(gpowcorr{t}.tscore(1, cfg.gpowcorr.freq{f}(1), :)); % take one example, lowest frequency
        cfg5.xlim = gpowcorr{t}.time(~isnan(onedat));
        
        ft_topoplotER(cfg5, gplot);
        %--------%
        %-----------------%
        
        %-----------------%
        %-save and link
        freqname = sprintf('f%02.f-%02.f', cfg.gpowcorr.freq{f}(1), cfg.gpowcorr.freq{f}(2));
        
        pngname = sprintf('gpowcorr_topo_%s_%s', condname, freqname);
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
        cfg4.zlim = cfg.gpow.freq{f};
        ft_singleplotER(cfg4, gtime{:});
        
        legend('cond1', 'cond2')
        ylabel(cfg4.parameter)
        
        freqname = sprintf('f%02.f-%02.f', cfg.gpowcorr.freq{f}(1), cfg.gpowcorr.freq{f}(2));
        title([condname ' ' cfg.gpowcorr.chan(c).name ' ' freqname])
        %--------%
        %-----------------%
        
        %-----------------%
        %-save and link
        pngname = sprintf('gpowcorr_val_%s_c%02.f_%s', condname, c, freqname);
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
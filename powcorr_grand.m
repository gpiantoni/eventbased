function powcorr_grand(cfg)
%POWCORR_GRAND grand power average
% 1) read single subject-data and create gpowcorr in cfg.dpow
% 2) do statistics for condition indicated by cfg.gpowcorr.stat.cond
% 3) plot the topoplot over time, frequency and singleplot for some electrodes
%
% CFG
%-Average
%  .nick: NICKNAME to save files specific to each NICKNAME
%  .log: name of the file and directory with analysis log
%  .subjall: index of the number of subjects
%
%  .dpow: directory with ERP data
%  .powcorr.cond: conditions to make averages
%
%-Statistics
%  .gpowcorr.stat.cond: cells within cell (e.g. {{'cond1' 'cond2'} {'cond1'} {'cond2'}})
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
%  [cfg.derp 'erp_SUBJ_COND']: timelock analysis for single-subject
%
% OUT
%  [cfg.dpow 'COND_grandpow']: power analysis for all subjects
%  [cfg.dpow 'COND_powpeak']: significant peaks in the POW
%
% FIGURES (saved in cfg.log and, if not empty, cfg.rslt)
%  gpow_tfr_CHAN_COND: time-frequency plot POW, for each condition, for one channel group
%  gpow_val_CHAN_FREQ: singleplot POW, all conditions, for one channel group, one frequency
%  gpow_topo_COND_FREQ: topoplot POW for each condition and frequency, over time
%
% Part of EVENTBASED group-analysis
% see also ERP_SUBJ, ERP_GRAND, ERPSOURCE_SUBJ, ERPSOURCE_GRAND, 
% POW_SUBJ, POW_GRAND, POWSOURCE_SUBJ, POWSOURCE_GRAND, 
% POWCORR_SUBJ, POWCORR_GRAND, POWSTAT_SUBJ, POWSTAT_GRAND, 
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
gpowcorr = [];
for k = 1:numel(cfg.powcorr.cond)
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
save([cfg.dpow cfg.nick '_grandpowcorr'], 'gpowcorr')
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
    
    save([cfg.dpow cfg.nick condname '_powcorrpeak'], 'powcorrpeak')
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
      pngname = sprintf('gpow_tfr_c%02.f_%s', cfg.gpowcorr.chan(c).name, condname);
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
        onedat = squeeze(gpowcorr{t}.tscore(1, cfg.gpowcorr.freq(f).freq(1), :)); % take one example, lowest frequency
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
        ft_singleplotER(cfg4, gtime{:});
        
        legend('cond1', 'cond2')
        ylabel(cfg4.parameter)
        
        title([condname ' ' cfg.gpowcorr.chan(c).name ' ' cfg.gpowcorr.freq(f).name])
        %--------%
        %-----------------%
        
        %-----------------%
        %-save and link
        pngname = sprintf('gpowcorr_val_%s_%s_%s', condname, cfg.gpowcorr.chan(c).name, cfg.gpowcorr.freq(f).name);
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
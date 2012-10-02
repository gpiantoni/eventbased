function conn_stat(info, opt)
%CONN_STAT statistics on connectivity analysis
% Attention: this does not handle one single condition case 
% 
% CFG
%  .cond: name to be used in projects/PROJNAME/subjects/0001/MOD/CONDNAME/
%  .dcon: directory with connectivity data
%  .conn.method: method used for connectivity
%  .log: name of the file and directory with analysis log
%  .rslt: directory images are saved into
% 
%  .conn.toi: vector with time points to run connectivity on
%  .statconn.bl.baseline: two scalars with baseline windows (if empty, no baseline)
%  .statconn.bl.baselinetype: type of baseline ('relchange' 'relative' 'absolute')
%
%  .statconn.ttest2: two index values of cfg.conn.test to compare directly
%  .statconn.time: cell with two scalars. Each cell gives the time window to test the t-test between conditions
% 
% OUT
%  [info.log connsum]: values from the plots to report in csv
% 
% FIGURES
%  gtrs_LABEL1_LABEL2_CONNMETHOD: connectivity over time, with one subplot per frequency
%
% Part of EVENTBASED single-subject
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

if isfield(cfg.gconn, 'comp')
  
  %-------------------------------------%
  %-loop over statistics conditions
  for t = 1:numel(cfg.gconn.comp)
    
    %---------------------------%
    %-statistics for effects of interest
    if numel(cfg.gpow.comp{t}) == 1
      % TODO
      
    else
      
      %-----------------%
      %-compare two conditions
      cond1 = cfg.gconn.comp{t}{1};
      cond2 = cfg.gconn.comp{t}{2};
      comp = [regexprep(cond1, '*', '') '_' regexprep(cond2, '*', '')];
      output = sprintf('%s\n   COMPARISON %s vs %s\n', output, cond1, cond2);
      
      %-------%
      %-conn over subj
      [outtmp data1 data2] = load_subj(cfg, 'conn', cfg.gconn.comp{t});
      output = [output outtmp];
      if isempty(data1) || isempty(data2); continue; end
      
      cfg1 = [];
      pow{1} = ft_freqgrandaverage(cfg1, data1{:});
      pow{2} = ft_freqgrandaverage(cfg1, data2{:});
      cfg1.keepindividual = 'yes';
      gpowall1 = ft_freqgrandaverage(cfg1, data1{:});
      gpowall2 = ft_freqgrandaverage(cfg1, data2{:});
      %-------%
      
      %-------%
      %-data for plot
      gplot = pow{2};
      if isempty(cfg.pow.bl)
        gplot.powspctrm = log(pow{2}.powspctrm ./ pow{1}.powspctrm);
      else % with baseline correction, take the difference
        gplot.powspctrm = pow{2}.powspctrm - pow{1}.powspctrm;
      end
      %-------%
      
      [pow_peak outtmp] = reportcluster(cfg, gpowall1, gpowall2);
      %-----------------%
      
    end
    
    save([cfg.dpow 'pow_peak_' comp], 'pow_peak')
    output = [output outtmp];
    %---------------------------%
  end
end

%---------------------------%
%-load data
load([cfg.dcon cfg.cond '_' cfg.conn.method '_grandconn'], 'gconn')
gshort = mean(mean(mean(mean(gconn.mat,3),4),5),6);
symm = all(all(gshort - gshort' > eps(10))); % check if matrix is symmetric

if symm
  output = sprintf('%s%s should be symmetric\n', output, cfg.conn.method);
else
  output = sprintf('%s%s should be asymmetric\n', output, cfg.conn.method);
end
%---------------------------%

%-------------------------------------%
%-loop over labels
sem = @(x) nanstd(x,[],3) / sqrt(size(x,3));
connsum = []; % summary of connectivity (for gosd2csv)
cnt = 0;

for chan1 = 1:numel(gconn.label)
  
  %-------%
  %-look in both directions if asymmetrical
  if symm
    nextchan = chan1 + 1;
  else
    nextchan = 1;
  end
  %-------%
  
  for chan2 = nextchan:numel(gconn.label)
    if chan1 ~= chan2 % for asymm, use all but this combination
      
      h = figure;
      nplot = numel(gconn.freq);
      nyplot = ceil(sqrt(nplot));
      nxplot = ceil(nplot./nyplot);
      
      %-----------------%
      %-subplot for frequency
      for f = 1:numel(gconn.freq)
        
        %-------%
        %-plot
        subplot(nxplot,nyplot,f)
        hold on
        
        x = squeeze(gconn.mat(chan1, chan2, :, f, :, :));
        
        if ~isempty(cfg.statconn.bl.baseline)
          x = performNormalization(cfg.conn.toi, x, cfg.statconn.bl.baseline, cfg.statconn.bl.baselinetype);
        end
        
        conntime = nanmean(x, 3); % connectivity over time
        % plot(cfg.conn.toi, conntime)
        xax = repmat(cfg.conn.toi, [numel(cfg.conn.test) 1]);
        errorbar(xax', conntime, sem(x))
        %-------%
        
        %-------%
        %-info summary
        cnt = cnt + 1;
        connsum(cnt).chan1 = gconn.label{chan1};
        connsum(cnt).chan2 = gconn.label{chan2};
        connsum(cnt).freq  = sprintf('% 3.f-% 3.f', gconn.freq{f}(1), gconn.freq{f}(2));
        connsum(cnt).cond1 = regexprep(cfg.conn.test{cfg.statconn.ttest2(1)}, '*', '');
        connsum(cnt).cond2 = regexprep(cfg.conn.test{cfg.statconn.ttest2(2)}, '*', '');
           
        for t = 1:numel(cfg.statconn.time)
          t1 = nearest(cfg.conn.toi, cfg.statconn.time{t}(1));
          t2 = nearest(cfg.conn.toi, cfg.statconn.time{t}(2));
          
          x1 = squeeze(x(t1:t2, cfg.statconn.ttest2(1), :))';
          x2 = squeeze(x(t1:t2, cfg.statconn.ttest2(2), :))';
          
          [~, p] = ttest(x1, x2);
          
          [~, minp] = min(p);
          connsum(cnt).time(t).time = cfg.conn.toi(t1 + minp - 1);
          connsum(cnt).time(t).pval = p(minp);
          
          %-assess significance: 0 and -1 is not significant, 1 and 2 are significant
          if p(minp) > .1; sig = -1; elseif p(minp) > .05; sig = 0; elseif p(minp) > .01; sig = 1; else sig = 2; end
          connsum(cnt).time(t).sign = sig;
          
          connsum(cnt).time(t).cond1 = mean(x1(:,minp));
          connsum(cnt).time(t).cond2 = mean(x2(:,minp));
        end
        %-------%
        
        %-------%
        title([connsum(cnt).chan1 ' ' connsum(cnt).chan2 ' ' connsum(cnt).freq])
        xlim(cfg.conn.toi([1 end]))
        xlabel('time (s)')
        ylabel(cfg.conn.method);
        legend(cfg.conn.test{:}, 'Location', 'NorthWest')
        %-------%
        
      end
      %-----------------%
      
      %--------%
      %-save and link
      pngname = sprintf('gtrs_%s_%s_%s', gconn.label{chan1}, gconn.label{chan2}, cfg.conn.method);
      saveas(gcf, [info.log filesep pngname '.png'])
      close(gcf); drawnow
      
      [~, logfile] = fileparts(info.log);
      system(['ln ' info.log filesep pngname '.png ' cfg.rslt pngname '_' logfile '.png']);
      %--------%
      
    end
  end
end

%-------------------------------------%

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

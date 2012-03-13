function conn_stat(cfg)
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
%  .statconn.time: cell with two scalars. Each cell gives the time window to test the t-test between conditions
% 
% OUT
%  [cfg.log connsum]: values from the plots to report in csv
% 
% FIGURES
%  gtrs_CONNMETHOD_LABEL1_LABEL2: connectivity over time, with one subplot per frequency
%
% Part of EVENTBASED group-analysis
% see also ERP_SUBJ, ERP_GRAND, ERPSOURCE_SUBJ, ERPSOURCE_GRAND, 
% POW_SUBJ, POW_GRAND, POWSOURCE_SUBJ, POWSOURCE_GRAND, 
% POWCORR_SUBJ, POWCORR_GRAND,
% CONN_SUBJ, CONN_GRAND, CONN_STAT

%---------------------------%
%-start log
output = sprintf('%s started at %s on %s\n', ...
  mfilename,  datestr(now, 'HH:MM:SS'), datestr(now, 'dd-mmm-yy'));
tic_t = tic;
%---------------------------%

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
sem = @(x) std(x,[],3) / sqrt(size(x,3));
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
        
        conntime = mean(x, 3); % connectivity over time
        % plot(cfg.conn.toi, conntime)
        xax = repmat(cfg.conn.toi, [numel(cfg.test) 1]);
        errorbar(xax', conntime, sem(x))
        %-------%
        
        %-------%
        %-info summary
        cnt = cnt + 1;
        connsum(cnt).chan1 = gconn.label{chan1};
        connsum(cnt).chan2 = gconn.label{chan2};
        connsum(cnt).freq  = sprintf('% 3.f-% 3.f', gconn.freq{f}(1), gconn.freq{f}(2));
        connsum(cnt).cond1 = regexprep(cfg.test{cfg.statconn.ttest2(1)}, '*', '');
        connsum(cnt).cond2 = regexprep(cfg.test{cfg.statconn.ttest2(2)}, '*', '');
           
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
        legend(cfg.test{:})
        %-------%
        
      end
      %-----------------%
      
      %--------%
      %-save and link
      pngname = sprintf('gtrs_%s_%s_%s', cfg.conn.method, gconn.label{chan1}, gconn.label{chan2});
      saveas(gcf, [cfg.log filesep pngname '.png'])
      close(gcf); drawnow
      
      [~, logfile] = fileparts(cfg.log);
      system(['ln ' cfg.log filesep pngname '.png ' cfg.rslt pngname '_' logfile '.png']);
      %--------%
      
    end
  end
end

save([cfg.log filesep 'connsum'], 'connsum')
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
fid = fopen([cfg.log '.txt'], 'a');
fwrite(fid, output);
fclose(fid);
%-----------------%
%---------------------------%

%-------------------------------------%
%-performNormalization (copied from ft_freqbaseline)
function data = performNormalization(timeVec, data, baseline, baselinetype)

baselineTimes = (timeVec >= baseline(1) & timeVec <= baseline(2));

% compute mean of time/frequency quantity in the baseline interval,
% ignoring NaNs, and replicate this over time dimension
meanVals = repmat(nanmean(data(baselineTimes,:,:), 1), [size(data, 1) 1 1]);

if (strcmp(baselinetype, 'absolute'))
  data = data - meanVals;
elseif (strcmp(baselinetype, 'relative'))
  data = data ./ meanVals;
elseif (strcmp(baselinetype, 'relchange'))
  data = (data - meanVals) ./ meanVals;
else
  error('unsupported method for baseline normalization: %s', baselinetype);
end
%-------------------------------------%

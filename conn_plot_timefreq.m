function conn_plot_timefreq(info, opt)
%CONN_PLOT_TIMEFREQ plot connectivity analysis when contains frequency info
% (XXX to test XXX)
% INFO
%  .log: name of the file and directory with analysis log
%  .dcon: directory with CONN data
%  .sens.layout: file with layout. It should be a mat containing 'layout'
%                If empty, it does not plot topo.
%  .rslt: directory images are saved into
%
% CFG.OPT
%  .conn.method*: (if .conn.type == 'cca') 'gc'; (if .conn.type == 'ft')
%                'coh', 'csd', 'plv', 'powcorr', 'amplcorr', 'ppc', 'psi'
%                (symmetric) or 'granger', 'dtf', 'pdc' (asymmetric)
%  .comp*: comparisons to test (cell within cell, e.g. {{'cond1' 'cond2'} {'cond1'} {'cond2'}})
%        but you cannot have more than 2 conditions (it's always a t-test). If empty, not statistics and no plots
%  .conn.toi*: vector with time points to run connectivity on
%
%  Baseline correction at the single-subject level:
%  .conn.bl: if empty, no baseline. Otherwise:
%  .conn.bl.baseline: two scalars with baseline windows
%  .conn.bl.baselinetype: type of baseline ('relchange')
%
% IN
%  [info.dcon COND_CONNMETHOD_GRANDCONN]: a matrix with all connectivity measures (chan X chan X time X freq X test X subj)
%
% * indicates obligatory parameter
%
% Part of EVENTBASED group-analysis
% see also ERP_SUBJ, ERP_GRAND,
% ERPSOURCE_SUBJ, ERPSOURCE_GRAND, ERPSTAT_SUBJ, ERPSTAT_GRAND,
% POW_SUBJ, POW_GRAND, POW_GRP, POWCORR_SUBJ, POWCORR_GRAND,
% POWSOURCE_SUBJ, POWSOURCE_GRAND, POWSTAT_SUBJ, POWSTAT_GRAND, SOURCE_SUBJ, 
% CONN_SUBJ, CONN_GRAND, CONN_PLOT_TIME, CONN_PLOT_TIMEFREQ, CONN_STAT

%---------------------------%
%-start log
output = sprintf('%s began at %s on %s\n', ...
  mfilename,  datestr(now, 'HH:MM:SS'), datestr(now, 'dd-mmm-yy'));
tic_t = tic;
%---------------------------%

%-----------------------------------------------%
%-compare conditions
sem = @(x) std(x,[],3) / sqrt(size(x,3));

if isfield(opt, 'comp')
  
  %-------------------------------------%
  %-loop over statistics conditions
  for t = 1:numel(opt.comp)
    
    %---------------------------%
    %-load the data
    if numel(opt.comp{t}) == 1
      
      %-----------------%
      %-one condition
      cond = regexprep(opt.comp{t}{1}, '*', '');
      comp = cond;
      output = sprintf('%s\n   COMPARISON %s\n', output, cond);
      load([info.dcon 'conn_' cond], 'conn')
      %-----------------%
      
    else
      
      %-----------------%
      %-two conditions
      cond1 = regexprep(opt.comp{t}{1}, '*', '');
      cond2 = regexprep(opt.comp{t}{2}, '*', '');
      comp = [cond1 '_' cond2];
      output = sprintf('%s\n   COMPARISON %s vs %s\n', output, cond1, cond2);
      
      load([info.dcon 'conn_' cond1], 'conn')
      conn1 = conn;
      load([info.dcon 'conn_' cond2], 'conn')
      conn2 = conn;
      
      conn.mat = log(conn1.mat ./ conn2.mat);
      clear conn1 conn2
      %-----------------%
      
    end
    %---------------------------%
    
    %---------------------------%
    %-load data
    gshort = mean(mean(mean(mean(conn.mat,3),4),5),6);
    symm = all(all(gshort - gshort' < eps(10))); % check if matrix is symmetric
    
    if symm
      output = sprintf('%s%s should be symmetric\n', output, opt.conn.method);
    else
      output = sprintf('%s%s should be asymmetric\n', output, opt.conn.method);
    end
    %---------------------------%
    
    %-------------------------------------%
    %-loop over labels
    nchan = numel(conn.label);
    
    h = figure('vis', 'off')
    for chan1 = 1:nchan
      
      %-----------------%
      %-look in both directions if asymmetrical
      if symm
        nextchan = chan1 + 1;
      else
        nextchan = 1;
      end
      %-----------------%
      
      for chan2 = nextchan:nchan
        if chan1 ~= chan2 % for asymm, use all but this combination
          
          subplot(nchan, nchan, (chan1 - 1) * nchan + chan2)
          dat = shiftdim(conn.mat(chan1, chan2, :, :, :), 2);
          
          if isfield(opt.conn, 'bl') && isfield(opt.conn.bl, 'baseline') ...
              && ~isempty(opt.conn.bl.baseline)
            dat = performNormalization(opt.conn.toi, dat, opt.conn.bl.baseline, opt.conn.bl.baselinetype);
          end
          
          %-----------------%
          %-plot image
          dat = mean(dat,3);
          if find(dat(:) < 0)
            clim = [-1 1] * max(abs(dat(:)));
          else
            clim = [ 0 1] * max(abs(dat(:)));
          end
          
          imagesc(cfg.conn.toi, 1:numel(conn.freq), dat', clim)
          
          axis xy
          colorbar
          freq2str = @(x)sprintf('%1.f-%1.f', x(1), x(2));
          set(gca, 'ytick', 1:numel(conn.freq), 'yticklabel', cellfun(freq2str, conn.freq, 'uni', 0))
          ylabel('frequency bands')
          %-----------------%
          
          %-----------------%
          %-title
          title([conn.label{chan1} ' -> ' conn.label{chan2}])
          if numel(conn.time) > 1
            xlim(conn.time([1 end]))
            xlabel('time (s)')
          end
          %-----------------%
          
        end % chan1 ~= chan2
      end % chan2
    end % chan1
    %-------------------------------------%
    
    %-----------------%
    %-save and link
    pngname = sprintf('gtrs_timefreq_%s_%s', opt.conn.method, comp);
    saveas(h, [info.log filesep pngname '.png'])
    close(h); drawnow
    
    [~, logfile] = fileparts(info.log);
    system(['ln ' info.log filesep pngname '.png ' info.rslt pngname '_' logfile '.png']);
    %-----------------%
    
  end % numel(gcomp)
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

%-------------------------------------%
%-performNormalization (copied from ft_freqbaseline)
function data = performNormalization(timeVec, data, baseline, baselinetype)

baselineTimes = (timeVec + eps(10) >= baseline(1) & timeVec - eps(10) <= baseline(2));

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

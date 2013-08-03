function conn_plot_time(info, opt)
%CONN_PLOT_TIME plot connectivity analysis over time only
%
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
%          you can have as many conditions as you want.
%  .conn.toi*: vector with time points to run connectivity on
%  .stat: if not empty, plot the significant/not-significant parts
%       .pvalue: pvalue for significance (default: 0.025)
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

%---------------------------%
%-check input
stat = false;
if isfield(opt, 'stat') && ~isempty(opt.stat)
  stat = true;
  if ~isfield(opt.stat, 'pvalue')
    opt.stat.pvalue = 0.05 / 2;
  end
end
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
    comp = '';
    output = sprintf('%s\n   COMPARISON', output);
    for c = 1:numel(opt.comp{t})
      
      %-----------------%
      %-conditions
      cond = regexprep(opt.comp{t}{c}, '*', '');
      comp = [comp '_' cond];
      output = sprintf('%s %s', output, cond);
      load([info.dcon 'conn_' cond], 'conn')
      
      if c == 1
        connall = conn;
      else
        connall.mat = cat(6, connall.mat, conn.mat); % CONN.MAT has dimension: label X label X time X freq X subj
      end
      %-----------------%
      
    end
    conn = connall;
    output = [output sprintf('\n')];
    %---------------------------%
    
    %---------------------------%
    %-
    if size(conn.mat,4) > 1
      error('You cannot use CONN_PLOT_TIME for data with multiple frequency.')
      %TODO: implement average over frequency or multiple plots per
      %frequency with CONN_PLOT_TIME
    end
    %---------------------------%
    
    %---------------------------%
    %-test if symmetric
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
    
    h = figure('vis', 'off');
    h_subplot = [];
    maxy = 0;
    
    if stat
      df = size(conn.mat,5) - 1; % always tested against zero
      t_thresh = abs(tinv(opt.stat.pvalue, df));
      
      Sh = figure('vis', 'off');
      Sh_subplot = [];
      Smaxy = 0;
    end
    
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
          
          figure(h)
          h_subplot = [h_subplot; subplot(nchan, nchan, (chan1 - 1) * nchan + chan2)];
          hold on
          dat = shiftdim(conn.mat(chan1, chan2, :, :, :, :), 2);
          
          if isfield(opt.conn, 'bl') && isfield(opt.conn.bl, 'baseline') ...
              && ~isempty(opt.conn.bl.baseline)
            for c = 1:size(dat, 4)
              dat(:,:,:,c) = performNormalization(opt.conn.toi, dat(:,:,:,c), opt.conn.bl.baseline, opt.conn.bl.baselinetype);
            end
          end
          
          %-----------------%
          %-only one frequency, plot line with error bar
          conntime = repmat(conn.time', 1, size(dat,4));
          m_dat = permute(mean(dat,3), [1 4 2 3]);
          sem_dat = permute(sem(dat), [1 4 2 3]);
          errorbar(conntime, m_dat, sem_dat);
          ylabel(opt.conn.method);
          
          maxy = max(maxy, nanmax(abs(m_dat(:))));
          %-----------------%

          %-----------------%
          %-title
          title([conn.label{chan1} ' -> ' conn.label{chan2}])
          if numel(conn.time) > 1
            xlim(conn.time([1 end]))
            xlabel('time (s)')
          end
          %-----------------%

          %-----------------%
          if stat
            figure(Sh)
            Sh_subplot = [Sh_subplot; subplot(nchan, nchan, (chan1 - 1) * nchan + chan2)];
          
            signplot = sign(m_dat) .* (abs(m_dat ./ sem_dat) > t_thresh);
            plot(conntime, signplot, 'o')
            title([conn.label{chan1} ' -> ' conn.label{chan2}])

          end
          %-----------------%
          
        end % chan1 ~= chan2
      end % chan2
    end % chan1
    
    set(h_subplot, 'ylim', [-1 1] * maxy)
    
    if stat
      set(Sh_subplot, 'ylim', [-1.5 1.5])
    end
    %-------------------------------------%
    
    %-----------------%
    %-save and link
    pngname = sprintf('gtrs_time_%s_%s', opt.conn.method, comp);
    saveas(h, [info.log filesep pngname '.png'])
    close(h); drawnow
    
    [~, logfile] = fileparts(info.log);
    system(['ln ' info.log filesep pngname '.png ' info.rslt pngname '_' logfile '.png']);
    %-----------------%
    
    %-----------------%
    %-save and link
    if stat
      pngname = sprintf('gtrs_time_%s_%s_pvalue', opt.conn.method, comp);
      saveas(Sh, [info.log filesep pngname '.png'])
      close(Sh); drawnow
      
      [~, logfile] = fileparts(info.log);
      system(['ln ' info.log filesep pngname '.png ' info.rslt pngname '_' logfile '.png']);
    end
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

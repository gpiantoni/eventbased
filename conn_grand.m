function conn_grand(cfg)
%CONN_GRAND connectivity analysis across subjects
%
% CFG
%  .cond: name to be used in projects/PROJNAME/subjects/0001/MOD/CONDNAME/
%  .dcon: directory with connectivity data
%  .log: name of the file and directory with analysis log
%  .conn.test: a cell with the condition defined by redef.
%              It can, but need not, be identical to cfg.test
%  .conn.method: method used for connectivity
%  .conn.toi: vector with time points to run connectivity on
%  .gconn.freq:
%    if 'all': it takes all the frequency in the data (one value)
%    if 'any': it takes each frequency in the data (can be alot alot)
%    if two scalar: it takes each frequency between the extremes ([8 12], means each frequency between 8 and 12, so 8 9 10 11 12, five values)
%    if a cell with scalar: it takes the average between the two limits ({[8 12]}, means average of all the frequencies between 8 and 12, one value)
%
% OUT
%  [cfg.dcon COND_CONNMETHOD_GRANDCONN]: a matrix with all connectivity measures (chan X chan X time X freq X test X subj)
%
% Part of EVENTBASED group-analysis
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

%---------------------------%
%-preparation
%-----------------%
%-define name of spectrum
switch cfg.conn.method
  
  %-symmetric
  case 'coh'
    spctrm = 'cohspctrm';
  case 'csd'
    spctrm = 'crsspctrm';
  case 'plv'
    spctrm = 'plvspctrm';
  case 'powcorr'
    spctrm = 'powcorrspctrm';
  case 'amplcorr'
    spctrm = 'amplcorrspctrm';
  case 'ppc'
    spctrm = 'ppcspctrm';
  case 'psi'
    spctrm = 'psispctrm';
    
    %-asymmetric
  case 'granger'
    spctrm = 'grangerspctrm';
  case 'dtf'
    spctrm = 'dtfspctrm';
  case 'pdc'
    spctrm = 'pdcspctrm';
    
  case 'gc'
    %-cca toolbox
    spctrm = 'gc';
end
%-----------------%

%-----------------%
%-freq info
%-------%
%-load one example
condname = regexprep(cfg.conn.cond{1}, '*', ''); % it does not matter
outputfile = sprintf('conn_*_%s.mat', condname);

allsub = dir([cfg.dcon outputfile]);
load([cfg.dcon allsub(1).name], 'stat')
%-------%

if ischar(cfg.gconn.freq) && strcmp(cfg.gconn.freq, 'all')
  connfreq = [stat.freq(1) stat.freq(end)];
  
elseif ischar(cfg.gconn.freq) && strcmp(cfg.gconn.freq, 'any')
  connfreq = [stat.freq' stat.freq'];
  
elseif isnumeric(cfg.gconn.freq)
  freq1 = nearest(stat.freq, cfg.gconn.freq(1));
  freq2 = nearest(stat.freq, cfg.gconn.freq(2));
  connfreq = [stat.freq(freq1:freq2)' stat.freq(freq1:freq2)'];
  
elseif iscell(cfg.gconn.freq)
  connfreq = cell2mat(cfg.gconn.freq');
  
end
%-----------------%
%---------------------------%

%-----------------------------------------------%
%-read the data
%---------------------------%
%-loop over conditions
for k = 1:numel(cfg.conn.cond)
  
  cond     = cfg.conn.cond{k};
  condname = regexprep(cond, '*', '');
  
  %-----------------%
  %-file for each cond
  [outtmp data] = load_subj(cfg, 'conn', cond);
  output = [output outtmp];
  if isempty(data); continue; end
  %-----------------%
  
  %---------------------------%
  %-prepare conn structure
  % MAT is chan X chan X time X freq X subj
  conn.label = stat.label;
  conn.time = cfg.conn.toi;% conn.time = data{1}.time;
  conn.freq = mat2cell(connfreq, size(connfreq,1), [2]);
  conn.mat = nan(numel(conn.label), numel(conn.label), numel(conn.time), size(connfreq,1), numel(cfg.subjall));
  %---------------------------%
  
  %---------------------------%
  %-loop over frequencies
  for f = 1:size(connfreq,1)
    
    %--------%
    %-freq of interest
    foi = nearest(data{1}.freq, connfreq(f,1)) : nearest(data{1}.freq, connfreq(f,2));
    
    if isempty(foi)
      error(['no frequencies are selected'])
    end
    %--------%
    
    %-----------------%
    %-loop over subj
    % data{i} has dimord "chan_chan_freq_time"
    for i = 1:numel(data)
      conn.mat(:, :, :, f, i) = permute(mean( data{i}.(spctrm)(:,:, foi, :), 3), [1 2 4 3 5 6]);
    end
    %-----------------%
    
  end
  %---------------------------%
  
  %-----------------%
  %-save
  save([cfg.dcon 'conn_' condname], 'conn')
  clear conn
  %-----------------%
  
end
%-----------------------------------------------%

%-----------------------------------------------%
%-compare conditions
if isfield(cfg.gconn, 'comp')
  
  %-------------------------------------%
  %-loop over statistics conditions
  for t = 1:numel(cfg.gconn.comp)
    
    %---------------------------%
    %-compare two conditions
    %-----------------%
    %-load data
    cond1 = regexprep(cfg.gconn.comp{t}{1}, '*', '');
    cond2 = regexprep(cfg.gconn.comp{t}{2}, '*', '');
    comp = [cond1 '_' cond2];
    output = sprintf('%s\n   COMPARISON %s vs %s\n', output, cond1, cond2);
    
    load([cfg.dcon 'conn_' cond1], 'conn')
    conn1 = conn;
    load([cfg.dcon 'conn_' cond2], 'conn')
    conn2 = conn;
    
    clear conn
    %-----------------%
    %---------------------------%
    
    %---------------------------%
    %-load data
    gshort = mean(mean(mean(mean(conn1.mat,3),4),5),6);
    symm = all(all(gshort - gshort' < eps(10))); % check if matrix is symmetric
    
    if symm
      output = sprintf('%s%s should be symmetric\n', output, cfg.conn.method);
    else
      output = sprintf('%s%s should be asymmetric\n', output, cfg.conn.method);
    end
    %---------------------------%
    
    %-------------------------------------%
    %-loop over labels
    sem = @(x) std(x,[],2) / sqrt(size(x,2));
    cnt = 0;
    
    for chan1 = 1:numel(conn1.label)
      
      %-------%
      %-look in both directions if asymmetrical
      if symm
        nextchan = chan1 + 1;
      else
        nextchan = 1;
      end
      %-------%
      
      for chan2 = nextchan:numel(conn1.label)
        if chan1 ~= chan2 % for asymm, use all but this combination
          
          h = figure;
          nplot = numel(conn1.freq);
          nyplot = ceil(sqrt(nplot));
          nxplot = ceil(nplot./nyplot);
          
          %-----------------%
          %-subplot for frequency
          for f = 1:numel(conn1.freq)
            
            %-------%
            %-plot
            subplot(nxplot,nyplot,f)
            hold on
            
            x1 = permute(conn1.mat(chan1, chan2, :, f, :), [3 5 1 2 4]); % 1 2 4 are singleton dim
            x2 = permute(conn2.mat(chan1, chan2, :, f, :), [3 5 1 2 4]); % 1 2 4 are singleton dim
            
            if isfield(cfg.gconn, 'log') && cfg.gconn.log % DOC
              % the log can be justified. GC follows an F-distribution, which is the ratio of two chi-square
              % If you deal with ratios, it's better to take the log
              x1 = log(x1);
              x2 = log(x2);
            end
            
            if isfield(cfg.gconn, 'bl') && isfield(cfg.gconn.bl, 'baseline') ...
                && ~isempty(cfg.gconn.bl.baseline)
              x1 = performNormalization(cfg.conn.toi, x1, cfg.gconn.bl.baseline, cfg.gconn.bl.baselinetype);
              x2 = performNormalization(cfg.conn.toi, x2, cfg.gconn.bl.baseline, cfg.gconn.bl.baselinetype);
            end
            
            errorbar(conn1.time, mean(x1, 2), sem(x1)); hold on
            errorbar(conn2.time, mean(x2, 2), sem(x2), 'r')
            %-------%
            
            %-------%
            %-info summary
            %-TODO: useful for export2csv, see conn_stat old revisions
            %-------%
            
            %-------%
            title_freq  = sprintf('% 3.f-% 3.f', conn1.freq{f}(1), conn1.freq{f}(2));
            title([conn1.label{chan1} ' -> ' conn1.label{chan2} ' ' title_freq])
            if numel(conn1.time) > 1
              xlim(conn1.time([1 end]))
              xlabel('time (s)')
            end
            ylabel(cfg.conn.method);
            legend(cond1, cond2, 'Location', 'NorthWest')
            %-------%
            
          end
          %-----------------%
          
          %--------%
          %-save and link
          pngname = sprintf('gtrs_%s_%s_%s', conn1.label{chan1}, conn1.label{chan2}, cfg.conn.method);
          saveas(gcf, [cfg.log filesep pngname '.png'])
          close(gcf); drawnow
          
          [~, logfile] = fileparts(cfg.log);
          system(['ln ' cfg.log filesep pngname '.png ' cfg.rslt pngname '_' logfile '.png']);
          %--------%
          
        end % chan1 ~= chan2
      end % chan2
    end % chan1 
    %-------------------------------------%
    
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
fid = fopen([cfg.log '.txt'], 'a');
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

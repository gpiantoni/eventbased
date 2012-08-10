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
%  Baseline correction at the single-subject level:
%  .conn.bl: if empty, no baseline. Otherwise:
%  .conn.bl.baseline: two scalars with baseline windows
%  .conn.bl.baselinetype: type of baseline ('relchange')
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
  conn.time = cfg.conn.toi;
  conn.freq = mat2cell(connfreq, ones(size(connfreq,1),1), [2]);
  conn.mat = nan(numel(conn.label), numel(conn.label), numel(conn.time), size(conn.freq,1), numel(cfg.subjall));
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
sem = @(x) std(x,[],3) / sqrt(size(x,3));

if isfield(cfg.gconn, 'comp')
  
  %-------------------------------------%
  %-loop over statistics conditions
  for t = 1:numel(cfg.gconn.comp)
    
    %---------------------------%
    %-load the data
    if numel(cfg.gconn.comp{t}) == 1
      
      %-----------------%
      %-one condition
      cond = regexprep(cfg.gconn.comp{t}{1}, '*', '');
      comp = cond;
      output = sprintf('%s\n   COMPARISON %s\n', output, cond);
      load([cfg.dcon 'conn_' cond], 'conn')
      %-----------------%
      
    else
      
      %-----------------%
      %-two conditions
      cond1 = regexprep(cfg.gconn.comp{t}{1}, '*', '');
      cond2 = regexprep(cfg.gconn.comp{t}{2}, '*', '');
      comp = [cond1 '_' cond2];
      output = sprintf('%s\n   COMPARISON %s vs %s\n', output, cond1, cond2);
      
      load([cfg.dcon 'conn_' cond1], 'conn')
      conn1 = conn;
      load([cfg.dcon 'conn_' cond2], 'conn')
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
      output = sprintf('%s%s should be symmetric\n', output, cfg.conn.method);
    else
      output = sprintf('%s%s should be asymmetric\n', output, cfg.conn.method);
    end
    %---------------------------%
    
    %-------------------------------------%
    %-loop over labels
    nchan = numel(conn.label);
    
    figure
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
          
          if isfield(cfg.gconn, 'bl') && isfield(cfg.gconn.bl, 'baseline') ...
              && ~isempty(cfg.gconn.bl.baseline)
            dat = performNormalization(cfg.conn.toi, dat, cfg.gconn.bl.baseline, cfg.gconn.bl.baselinetype);
          end
          
          if size(dat, 2) == 1

            %-----------------%
            %-only one frequency, plot line with error bar
            errorbar(conn.time, mean(dat, 3), sem(dat));
            ylabel(cfg.conn.method);
            %-----------------%
            
          else
            
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
            
          end
          
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
    pngname = sprintf('gtrs_%s_%s', cfg.conn.method, comp);
    saveas(gcf, [cfg.log filesep pngname '.png'])
    close(gcf); drawnow
    
    [~, logfile] = fileparts(cfg.log);
    system(['ln ' cfg.log filesep pngname '.png ' cfg.rslt pngname '_' logfile '.png']);
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

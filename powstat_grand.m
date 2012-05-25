function powstat_grand(cfg)
%POWSTAT_GRAND group-analysis of POW source data
%
% CFG
%-Average
%  .nick: NICK to save files specific to each NICK
%  .log: name of the file and directory with analysis log
%  .subjall: index of the number of subjects
%
%  .dpow: directory with POW data
%
%-Statistics
%  .powstat.comp: cells within cell (e.g. {{'cond1' 'cond2'} {'cond1'} {'cond2'}})
%         but you cannot have more than 2 conditions (it's always a t-test)
%  .powstat.nai: use NAI calculated from baseline or not (logical)
%
%  .powsource.peaks: how to speficy peaks to analyze, 'manual' or 'pow_peak'
%          (peaks from grandpow) or 'powcorr_peak' (peaks from grandpowcorr)
%    if 'manual'
%      .powsource.pow_peak(1).name: string ('name_of_the_time_window')
%      .powsource.pow_peak(1).time: scalar (center of the time window in s)
%      .powsource.pow_peak(1).wndw: scalar (length of the time window in s)
%      .powsource.pow_peak(1).freq = 10; % center of the frequency
%      .powsource.pow_peak(1).band = 4; % width of the frequency band
%    if 'pow_peak' or 'powcorr_peak'
%      .powsource.refcomp: cell of string(s) of the comparison whose peaks
%                     will be localized (one of the cells of cfg.gpow.comp or cfg.gpowcorr.comp)
%
%  .powstat.clusterstatistics: 'maxsize' or 'max'
%  .powstat.clusteralpha: level to select sensors (default 0.05)
%                           it can be a string in format '5%' to take top 5 voxels and put them in a cluster.
%  .powstat.maxvox: max number of significant voxels to be used in soupeak
%  .powstat.param: 'pow' or 'nai'
%  .powstat.clusterthr: threshold to report clusters in output
%
%  .rslt: directory images are saved into
%
% IN
%  [cfg.dpow 'powstat_SUBJ_COND']: source data for period of interest and baseline for each subject
%
% OUT
%  [cfg.dpow 'powstat_COND']: source analysis for all subject
%  [cfg.dpow 'powsoupeak_COND']: significant source peaks in POW
%
% FIGURES
%  gpow_peak_COND_POWPEAK: 3d plot of the source for one peak
%
% Part of EVENTBASED group-analysis
% see also ERP_SUBJ, ERP_GRAND,
% ERPSOURCE_SUBJ, ERPSOURCE_GRAND, ERPSTAT_SUBJ, ERPSTAT_GRAND,
% POW_SUBJ, POW_GRAND, POWCORR_SUBJ, POWCORR_GRAND,
% POWSOURCE_SUBJ, POWSOURCE_GRAND, POWSTAT_SUBJ, POWSTAT_GRAND,
% CONN_SUBJ, CONN_GRAND, CONN_STAT

%---------------------------%
%-start log
output = sprintf('%s began at %s on %s\n', ...
  mfilename,  datestr(now, 'HH:MM:SS'), datestr(now, 'dd-mmm-yy'));
tic_t = tic;
%---------------------------%

pow_peak = getpeak(cfg, 'pow');

%-------------------------------------%
%-loop over statistics conditions
for t = 1:numel(cfg.powstat.comp) % DOC: cfg.powstat.comp
  
  %---------------------------%
  %-statistics for effects of interest
  if numel(cfg.powstat.comp{t}) == 1
    
    %-----------------%
    %-compare against baseline
    cond = cfg.powstat.comp{t}{1};
    comp = regexprep(cond, '*', '');
    output = sprintf('%s\n   COMPARISON %s\n', output, cond);
    
    %-------%
    %-pow over subj
    [outtmp data] = load_subj(cfg, 'powstat', cond);
    output = [output outtmp];
    if isempty(data); continue; end
    %-------%
    
    sousubj1 = squeeze(data(:,2,:)); % period of interest
    sousubj2 = squeeze(data(:,1,:)); % baseline
    param = 'pow';
    %-----------------%
    
  else
    
    %-----------------%
    %-compare two conditions (cond1 - cond2)
    cond1 = cfg.powstat.comp{t}{1};
    cond2 = cfg.powstat.comp{t}{2};
    comp = [regexprep(cond1, '*', '') '_' regexprep(cond2, '*', '')];
    output = sprintf('%s\n   COMPARISON %s vs %s\n', output, cond1, cond2);
    
    %-------%
    %-pow over subj
    [outtmp data1] = load_subj(cfg, 'powstat', cond1);
    output = [output outtmp];
    if isempty(data1); continue; end
    %-------%
    
    %-------%
    %-pow over subj
    [outtmp data2] = load_subj(cfg, 'powstat', cond2);
    output = [output outtmp];
    if isempty(data2); continue; end
    %-------%
    param = 'pow';
    
    if isfield(cfg.powstat, 'nai') && cfg.powstat.nai
      for i1 = 1:size(data1,1)
        for i3 = 1:size(data1,3)
          data1{i1,2,i3}.avg.nai = data1{i1,2,i3}.avg.pow ./ data1{i1,1,i3}.avg.pow;
          data1{i1,2,i3}.avg = rmfield(data1{i1,2,i3}.avg, 'pow');
          
          data2{i1,2,i3}.avg.nai = data2{i1,2,i3}.avg.pow ./ data2{i1,1,i3}.avg.pow;
          data2{i1,2,i3}.avg = rmfield(data2{i1,2,i3}.avg, 'pow');
        end
      end
      param = 'nai';
    end
    
    sousubj1 = squeeze(data1(:,2,:)); % cond1
    sousubj2 = squeeze(data2(:,2,:)); % cond2
    %-----------------%
    
  end
  clear data*
  %---------------------------%
  
  %---------------------------%
  %-loop over peaks
  soupeak = [];
  for p = 1:numel(pow_peak)
    output = sprintf('%s\n%s:\n', output, pow_peak(p).name);

    %-----------------%
    %-grand average
    cfg1 = [];
    cfg1.keepindividual = 'yes';
    cfg1.parameter = param;
    soustat1 = ft_sourcegrandaverage(cfg1, sousubj1{:,p});
    soustat2 = ft_sourcegrandaverage(cfg1, sousubj2{:,p});
    %-----------------%    
    
    %-----------------%
    %-do stats and figure
    h = figure;
    [soupos powstat{p} outtmp] = reportsource(cfg.powstat, soustat1, soustat2); % DOC: cfg.powstat
    soupeak(p).pos = soupos;
    soupeak(p).center = mean(soupos,1);
    soupeak(p).name = pow_peak(p).name;
    output = [output outtmp];
    
    %--------%
    pngname = sprintf('gpowstat_%s_%s', comp, pow_peak(p).name);
    saveas(gcf, [cfg.log filesep pngname '.png'])
    close(gcf); drawnow
    
    [~, logfile] = fileparts(cfg.log);
    system(['ln ' cfg.log filesep pngname '.png ' cfg.rslt pngname '_' logfile '.png']);
    %--------%
    %-----------------%
    
  end
  %---------------------------%
  
  %-----------------%
  %-save
  save([cfg.dpow 'powstatpeak_' condname], 'soupeak')
  
  for p = 1:numel(powstat)
    powstat{p}.cfg = []; % this is huge
  end
  save([cfg.dpow 'powstat_' condname], 'powstat', '-v7.3')
  %-----------------%
  
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
fid = fopen([cfg.log '.txt'], 'a');
fwrite(fid, output);
fclose(fid);
%-----------------%
%---------------------------%
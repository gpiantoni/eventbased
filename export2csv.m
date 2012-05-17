function export2csv(cfg)
%EXPORT2CSV write results in csv file
% It will write the cfg, erppeak (and soupeak), powpeak (and soupeak), and
% results from connectivity analysis.
%
% You add columns with specific results. You need to specify:
%   .export2csv.extrainfo = 'functionname'
% where 'functionname' is as:
%   [output] = functionname(cfg)
% and output should contain no new lines but commas to separate values
%
% Part of EVENTBASED group-analysis
% see also ERP_SUBJ, ERP_GRAND, ERPSOURCE_SUBJ, ERPSOURCE_GRAND,
% POW_SUBJ, POW_GRAND, POWSOURCE_SUBJ, POWSOURCE_GRAND,
% POWCORR_SUBJ, POWCORR_GRAND, POWSTAT_SUBJ, POWSTAT_GRAND,
% CONN_SUBJ, CONN_GRAND, CONN_STAT

mincol = 10; % minimum columns
mincolcfg = 150; % minimum columns for the whole cfg
npeaks = 2; % max # of peaks to plot

%---------------------------%
%-get log file name and steps which were run
[~, logfile] = fileparts(cfg.log);
output = [logfile ','];
output = [output sprintf(' %s', cfg.step{cfg.run}) ','];
%---------------------------%

%-------------------------------------%
%-loop over steps
%---------------------------%
%-POW and POWSOURCE
%-----------------%
%-CFG
%-------%
%-pow
output = [output realign(struct2log(cfg.pow, 'csv'), mincol)];
%-------%

%-------%
%-pow_grand
if isfield(cfg.gpow, 'stat') && isfield(cfg.gpow.stat, 'time')&& ~isempty(cfg.gpow.stat.time)
  output = [output sprintf('[%4.2f %4.2f],', cfg.gpow.stat.time(1), cfg.gpow.stat.time(2))];
else
  output = [output ' ,'];
end

if isfield(cfg.gpow, 'stat') && isfield(cfg.gpow.stat, 'freq') && ~isempty(cfg.gpow.stat.freq)
  output = [output sprintf('[%4.2f %4.2f],', cfg.gpow.stat.freq(1), cfg.gpow.stat.freq(2))];
else
  output = [output ' ,'];
end
%-------%
%-----------------%

%-----------------%
%-loop over comparisons
for t = 1:numel(cfg.gpow.comp)

  %-------%
  %-load powpeak
  cond1 = cfg.gpow.comp{t}{1};
  comp = regexprep(cond1, '*', '');
  if numel(cfg.gpow.comp{t}) == 2
    cond2 = cfg.gpow.comp{t}{2};
    comp = [comp '_' regexprep(cond2, '*', '')];
  end
  output = [output comp ','];
  
  load([cfg.dpow 'powpeak_' comp], 'powpeak')
  %-------%
  
  %-------%
  %-loop over first peaks
  
  for p = 1:npeaks
    if  p <= numel(powpeak)
      
      output = [output sprintf('p%d, %s, %1.3f, %1.3f, %1.3f, %1.3f, ', ...
        p, powpeak(p).name, powpeak(p).time, powpeak(p).wndw, powpeak(p).freq, powpeak(p).band)];
      output = [output sprintf('%1.4f,', powpeak(p).pval)]; % TODO: if manually specified, pvalue does not exist
      
    else
      output = [output ' , , , , , , ,'];
    end
  end
  %-------%
  
end
%-----------------%
%---------------------------%
%-------------------------------------%

%-------------------------------------%
%-read extra, dataset-specific info
if isfield(cfg, 'export2csv') && isfield(cfg.export2csv, 'extrainfo') && ...
    ~isempty(cfg.export2csv.extrainfo)
  outtmp = feval(cfg.export2csv.extrainfo, cfg);
  output = [output outtmp];
end
%-------------------------------------%

%-------------------------------------%
%-write to file
fid = fopen(cfg.csvf, 'a+');
fwrite(fid, output);
fprintf(fid, '\n');
fclose(fid);
%-------------------------------------%

%-------------------------------------%
%-add commas at the end so that columns are realigned
function csvlog = realign(csvlog, mincol)

ncol = mincol - numel(strfind(csvlog, ','));
if ncol <= 0
  return
else
  csvlog = [csvlog repmat(',', [1 ncol])];
end
%-------------------------------------%

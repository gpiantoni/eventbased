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
% TODO:
%   - powcorr
%   - testing
%
% Part of EVENTBASED

mincol = 10; % minimum columns
mincolcfg = 150; % minimum columns for the whole cfg

%-----------------%
%-get log file name and steps which were run
[~, logfile] = fileparts(cfg.log);
output = [logfile ','];
output = [output sprintf(' %s', cfg.step{cfg.run}) ','];
%-----------------%

%-------------------------------------%
%-loop over steps
if isfield(cfg, 'vol')
  output = [output sprintf('%s,', cfg.vol.type)];
else
  output = [output ','];
end

for st = 1:numel(cfg.step)
  output = [output cfg.step{st} ','];
  
  %---------------------------%
  %-preprocessing
  switch cfg.step{st}
    case 'seldata'
      output = [output sprintf('%1.f,', numel(cfg.seldata))];
      
    case 'gclean'
      output = [output sprintf('%1f,', cfg.gtool.bad_samples.MADs)];
      output = [output sprintf('%1f,', cfg.gtool.bad_channels.MADs)];
      output = [output sprintf('%1f,', cfg.gtool.eog.correction)];
      output = [output sprintf('%1f,', cfg.gtool.emg.correction)];
      
    case 'preproc'
      output = [output realign(struct2log(cfg.preproc, 'csv'), mincol)];
      
    case 'redef'
      output = [output realign(struct2log(cfg.preproc, 'csv'), mincol)];
  end
  %---------------------------%
  
  %---------------------------%
  %-analysis
  switch cfg.step{st}
    case 'erp_subj'
      output = [output realign(struct2log(cfg.erp, 'csv'), mincol)];
    case 'erp_grand'
      output = [output sprintf(' %f', cfg.erpeffect) ','];
      
    case 'erpsource_subj'
      output = [output realign(struct2log(cfg.erpsource.erp, 'csv'), mincol)];
      output = [output sprintf('%s,', cfg.erpsource.areas)];
      output = [output sprintf('%f,', cfg.erpsource.bline)];
      output = [output sprintf('%s,', cfg.erpsource.lambda)];
      output = [output sprintf('%s,', cfg.erpsource.powmethod)];
    case 'erpsource_grand'
      output = [output sprintf('%s,', cfg.erpsource.clusterstatistics)];
      output = [output sprintf('%f,', cfg.erpsource.clusteralpha)];
      output = [output sprintf('%f,', cfg.erpsource.maxvox)];
      
    case 'pow_subj'
      output = [output realign(struct2log(cfg.pow, 'csv'), mincol)];
    case 'pow_grand'
      output = [output sprintf(' %f', cfg.poweffect) ','];
      
    case 'powsource_subj'
      output = [output sprintf('%s,', cfg.powsource.areas)];
      output = [output sprintf('%f,', cfg.powsource.bline)];
      output = [output sprintf('%s,', cfg.powsource.lambda)];
      output = [output sprintf('%s,', cfg.powsource.powmethod)];
    case 'powsource_grand'
      output = [output sprintf('%s,', cfg.powsource.clusterstatistics)];
      output = [output sprintf('%f,', cfg.powsource.clusteralpha)];
      output = [output sprintf('%f,', cfg.powsource.maxvox)];
      
    case 'powcorr_subj'
      output = [output realign(struct2log(cfg.powcorr, 'csv'), mincol)];
    case 'powcorr_grand'
      output = [output sprintf('%f,', cfg.powcorreffect)];
      
    case 'conn_subj'
      output = [output realign(struct2log(cfg.conn, 'csv'), 2*mincol)];
    case 'conn_grand'
      output = [output struct2log(cfg.gconn, 'csv')];
    case 'conn_grand'
      output = [output realign(struct2log(cfg.statconn, 'csv'), mincol)];
      
  end
  %---------------------------%
  
end

output = realign(output, mincolcfg);
%-------------------------------------%

%-------------------------------------%
%-main findings
%---------------------------%
%-ERP
if numel(dir(cfg.derp)) > 2 % dir is not empty
  
  %-----------------%
  %-erppeak and soupeak
  %-------%
  %-get erppeak
  if isfield(cfg, 'erpsource') && ... % in case it's not specified at all
      strcmp(cfg.erpsource.areas, 'manual')
    erppeak = {cfg.erpsource.erppeak};
    
  else % if strcmp(cfg.erpsource.areas, 'erppeak')
    
    for p = cfg.erpeffect
      condname = regexprep(cfg.test{p}, '*', '');
      if exist([cfg.derp cfg.cond condname '_erppeak.mat'], 'file')
        load([cfg.derp cfg.cond condname '_erppeak'], 'erppeak')
        erppeakall{p} = erppeak;
      end
    end
    
    erppeak = erppeakall;
  end
  %-------%
  
  %-------%
  %-load ERP source peaks
  sou = false;
  
  if exist([cfg.derp cfg.cond '_soupeak.mat'], 'file')
    load([cfg.derp cfg.cond '_soupeak'], 'soupeak')
    
    if numel({erppeak{1}.name}) == numel({soupeak.name}) && ...
        all(strcmp({erppeak{1}.name}, {soupeak.name}))
      sou = true;
    end
  end
  %-------%
  %-----------------%
  
  %-----------------%
  %-loop over peaks (sorted by time)
  for p = 1:numel(erppeak)
    [~, peaktime] = sort([erppeak{p}.time]);
    
    for i = peaktime
      output = [output sprintf('p%d, %s, %1.3f, %1.3f,%1.3f,', ...
        p, erppeak{p}(i).name, erppeak{p}(i).pval, erppeak{p}(i).time, erppeak{p}(i).wndw)];
      if strcmp(cfg.erpsource.areas, 'erppeak')
        output = [output sprintf('%1.3f,', erppeak{p}(i).pval)];
      else
        output = [output ','];
      end
      
      if sou
        output = [output sprintf('[%1.2f %1.2f %1.2f], %1.f,', ...
          soupeak(i).center(1), soupeak(i).center(2), soupeak(i).center(3), size(soupeak(i).pos,1))];
      else
        output = [output ',,'];
      end
    end
  end
  %-----------------%
  
end
output = realign(output, mincol*4); % to be tested
%-----------------%
%---------------------------%

%---------------------------%
%-POW
if numel(dir(cfg.dpow)) > 2 % dir is not empty
  
  %-----------------%
  %-powpeak and soupeak
  %-------%
  %-get powpeak
  if isfield(cfg, 'powsource') && ... % in case it's not specified at all
      strcmp(cfg.powsource.areas, 'manual')
    powpeak = {cfg.powsource.powpeak};
    
  else % if strcmp(cfg.powsource.areas, 'powpeak')
    
    for p = cfg.poweffect
      condname = regexprep(cfg.test{p}, '*', '');
      if exist([cfg.dpow cfg.cond condname '_powpeak.mat'], 'file')
        load([cfg.dpow cfg.cond condname '_powpeak'], 'powpeak')
        powpeakall{p} = powpeak;
      end
    end
    
    powpeak = powpeakall;
  end
  %-------%
  
  %-------%
  %-load POW source peaks
  sou = false;
  
  if exist([cfg.dpow cfg.cond '_soupeak.mat'], 'file')
    load([cfg.dpow cfg.cond '_soupeak'], 'soupeak')
    
    if numel({powpeak{1}.name}) == numel({soupeak.name}) && ... % only the first poweffect has source
        all(strcmp({powpeak{1}.name}, {soupeak.name}))
      sou = true;
    end
  end
  %-------%
  %-----------------%
  
  %-----------------%
  %-loop over peaks (sorted by time)
  for p = 1:numel(powpeak)
    [~, peaktime] = sort([powpeak{p}.time]);
    
    for i = peaktime
      output = [output sprintf('p%d, %s,%1.3f,%1.3f,%1.3f,%1.3f,%1.3f,', ...
        p, powpeak{p}(i).name, powpeak{p}(i).pval, powpeak{p}(i).time, powpeak{p}(i).wndw, powpeak{p}(i).freq, powpeak{p}(i).band)];
      if isfield(cfg, 'powsource') && ... % in case it's not specified at all
          strcmp(cfg.powsource.areas, 'powpeak')
        output = [output sprintf('%1.4f,', powpeak{p}(i).pval)];
      else
        output = [output ','];
      end
      
      if sou && p == 1 % only the first poweffect has source
        output = [output sprintf('[%1.2f %1.2f %1.2f], %1.f,', ...
          soupeak(i).center(1), soupeak(i).center(2), soupeak(i).center(3), size(soupeak(i).pos,1))];
      else
        output = [output ',,'];
      end
    end
  end
  %-----------------%
  
end
output = realign(output, mincol*4);
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
%-connectivity analysis (and write)
% one row per channel combination/ or one row for the whole analysis if no connectivity
fid = fopen(cfg.csvf, 'a+');

if exist([cfg.log filesep 'connsum.mat'], 'file')
  
  %---------------------------%
  %-loop over chan-chan
  load([cfg.log filesep 'connsum'], 'connsum')
  
  for i = 1:numel(connsum) % channel combinations
    
    fwrite(fid, output);
    fprintf(fid, '%s, %s, %s, %s, %s,', ...
      connsum(i).cond1, connsum(i).cond2, connsum(i).freq, connsum(i).chan1, connsum(i).chan2);
    
    %-----------------%
    %-loop over time of interest (beginning, middle, end of the trial)
    for t = 1:numel(connsum(i).time)
      fprintf(fid, '%1.2f, %1.f, %1.3f, %1.4f, %1.4f,', ...
        connsum(i).time(t).time, connsum(i).time(t).sign, connsum(i).time(t).pval, connsum(i).time(t).cond1, connsum(i).time(t).cond2);
    end
    %-----------------%
    
    fprintf(fid, '\n');
  end
  %---------------------------%
  
else
  
  fwrite(fid, output);
  fprintf(fid, '\n');
  
end
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

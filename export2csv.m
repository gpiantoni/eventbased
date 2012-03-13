function export2csv(cfg)
%EXPORT2CSV write results in csv file
% 
% CFG
%  
%




mincol = 10; % minimum columns for erp and pow, to keep them aligned

%-------------------------------------%
%-simplified cfg
[~, logfile] = fileparts(cfg.log);
output = [logfile ','];

%-----------------%
%-trial definition
output = [output struct2log(cfg.preproc, 'csv')];
output = [output struct2log(cfg.redef, 'csv')];

output = [output sprintf('%s,', cfg.voltype)];
output = [output ',,,,']; % leave some extras in case we add options
%-----------------%
%-------------------------------------%

%-------------------------------------%
%-main findings
%---------------------------%
%-ERP
%-----------------%
%-output parameters
output = [output 'ERP,'];
output = [output sprintf('%s,', cfg.test{cfg.erpeffect})];

output = [output sprintf('%1.f,', cfg.erp.preproc.baselinewindow(1))];
output = [output sprintf('%1.f,', cfg.erp.preproc.baselinewindow(2))];

if strcmp(cfg.erp.preproc.lpfilter, 'no')
  output = [output 'no,'];
else
  output = [output sprintf('%1.f,', cfg.erp.preproc.lpfreq)];
end

output = [output ',,'];
output = [output 'ERPsource,'];

output = [output sprintf('%f,', cfg.gerp.bline(1))];
output = [output sprintf('%f,', cfg.gerp.bline(2))];

output = [output sprintf('%s,', cfg.erpsource.lambda)];
output = [output sprintf('%s,', cfg.erpsource.powmethod)];
output = [output sprintf('%f,', cfg.erpsource.bline)];

output = [output sprintf('%s,', cfg.erpsource.areas )];
output = [output ',,,,'];
%-----------------%

%-----------------%
%-use predefined or ERP-peaks for areas of interest
erp = false;

if strcmp(cfg.erpsource.areas, 'manual')
  erppeak = cfg.erpsource.erppeak;
  erp = true;
  
elseif strcmp(cfg.erpsource.areas, 'erppeak')
  if exist([cfg.derp cfg.proj '_erppeak.mat'], 'file')
    load([cfg.derp cfg.proj '_erppeak'], 'erppeak')
    erp = true;
  end
end

%-------%
%-load ERP source peaks
sou = false;

if exist([cfg.derp cfg.proj '_soupeak.mat'], 'file')
  load([cfg.derp cfg.proj '_soupeak'], 'soupeak')
  
  if numel({erppeak.name}) == numel({soupeak.name}) && ...
      all(strcmp({erppeak.name}, {soupeak.name}))
    sou = true;
  end
end
%-------%
%-----------------%

if erp
  
  %-----------------%
  %-loop over peaks (sorted by time)
  [~, peaktime] = sort([erppeak.time]);
  
  for i = peaktime
    output = [output sprintf('%s,%1.3f,%1.3f,', ...
      erppeak(i).name, erppeak(i).time, erppeak(i).wndw)];
    if strcmp(cfg.erpsource.areas, 'erppeak')
      output = [output sprintf('%1.3f,', erppeak(i).pval)];
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
  %-----------------%
  
end

%-----------------%
%-add extra column for realignment
if erp
  extracol = (mincol - numel(erppeak)) * 6; % there are six cells in ERP
else
  extracol = mincol * 6;
end
output = [output sprintf(repmat(' ,', 1, extracol))];
%-----------------%
%---------------------------%

%---------------------------%
%-POW
%-----------------%
%-output
output = [output 'POW,'];
output = [output sprintf('%s,', cfg.test{cfg.poweffect})];

output = [output sprintf('%s,', cfg.pow.method)];
output = [output sprintf('%s,', cfg.pow.taper)];

output = [output sprintf('%f,', cfg.pow.foi(1))];
output = [output sprintf('%f,', cfg.pow.foi(end))];

output = [output sprintf('%f,', cfg.pow.t_ftimwin(1))];
output = [output sprintf('%f,', cfg.pow.t_ftimwin(end))];

output = [output sprintf('%f,', cfg.pow.toi(1))];
output = [output sprintf('%f,', cfg.pow.toi(end))];

if ~isempty(cfg.pow.bl.baseline)
  output = [output sprintf('%s,', cfg.pow.bl.baselinetype)];
  output = [output sprintf('%f,', cfg.pow.bl.baseline(1))];
  output = [output sprintf('%f,', cfg.pow.bl.baseline(end))];
else
  output = [output ',,,'];
end

output = [output ',,'];
output = [output 'POWsource,'];

output = [output sprintf('%s,', cfg.powsource.lambda)];
output = [output sprintf('%s,', cfg.powsource.powmethod)];
output = [output sprintf('%f,', cfg.powsource.bline)];

output = [output sprintf('%s,', cfg.powsource.areas )];
output = [output ',,,,'];
%-----------------%

%-----------------%
%-use predefined or POW-peaks for areas of interest
pow = false;

if strcmp(cfg.powsource.areas, 'manual')
  powpeak = cfg.powsource.powpeak;
  pow = true;
  
elseif strcmp(cfg.powsource.areas, 'powpeak')
  if exist([cfg.dpow cfg.proj '_powpeak.mat'], 'file')
    load([cfg.dpow cfg.proj '_powpeak'], 'powpeak')
    pow = true;
  end
end

%-------%
%-load POW source peaks
sou = false;

if exist([cfg.dpow cfg.proj '_soupeak.mat'], 'file')
  load([cfg.dpow cfg.proj '_soupeak'], 'soupeak')
  
  if numel({powpeak.name}) == numel({soupeak.name}) && ... 
    all(strcmp({powpeak.name}, {soupeak.name}))
    sou = true;
  end
end
%-------%
%-----------------%

if pow
  
  %-----------------%
  %-loop over peaks (sorted by time)
  [~, peaktime] = sort([powpeak.time]);
  
  for i = peaktime
    output = [output sprintf('%s,%1.3f,%1.3f,%1.3f,%1.3f,', ...
      powpeak(i).name, powpeak(i).time, powpeak(i).wndw, powpeak(i).freq, powpeak(i).band)];
    if strcmp(cfg.powsource.areas, 'powpeak')
      output = [output sprintf('%1.4f,', powpeak(i).pval)];
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
  %-----------------%
  
end

%-----------------%
%-add extra column for realignment
if pow
  extracol = (mincol - numel(powpeak)) * 8; % there are six cells in pow
else
  extracol = mincol * 8;
end
output = [output sprintf(repmat(' ,', 1, extracol))];
%-----------------%
%---------------------------%
%-------------------------------------%

%-------------------------------------%
%-read extra, dataset-specific info
if isfield(cfg, 'export2csv') && isfield(cfg.export2csv, 'extrainfo')
  outtmp = feval(cfg.export2csv.extrainfo, cfg);
  output = [output outtmp];
end
%-------------------------------------%

%-------------------------------------%
%-connectivity analysis (and write)
% one row per channel combination/ or one row for the whole analysis if no connectivity
fid = fopen(cfg.csvf, 'a+');

if exist([cfg.log filesep 'connsum.mat'], 'file')
 
  %-----------------%
  %-output
  output = [output 'CONN,'];
  
  output = [output sprintf('%s,', cfg.conn.areas)];
  
  switch cfg.conn.areas
    case {'erppeak' 'powpeak'}
      output = [output sprintf('%s,', cfg.conn.fixedmom)];
    otherwise
      output = [output ','];
  end
  
  output = [output sprintf('%f,', cfg.conn.toi(1))];
  output = [output sprintf('%f,', cfg.conn.toi(end))];
  
  output = [output sprintf('%f,', cfg.conn.t_ftimwin)];
  
  if ~isempty(cfg.statconn.bl.baseline)
    output = [output sprintf('%s,', cfg.statconn.bl.baselinetype)];
    output = [output sprintf('%f,', cfg.statconn.bl.baseline(1))];
    output = [output sprintf('%f,', cfg.statconn.bl.baseline(end))];
  else
    output = [output ',,,'];
  end
  
  output = [output sprintf('%s,', cfg.conn.type)];
  
  switch cfg.conn.type
    case 'ft'
      output = [output sprintf('%s,', cfg.conn.mvar)];
      output = [output sprintf('%s,', cfg.conn.toolbox)];
      
      output = [output sprintf('%s,', cfg.conn.freq)];
      output = [output sprintf('%s,', cfg.conn.method)];
      
      
    case 'cca'
      output = [output sprintf('%s,', cfg.conn.method)];
      output = [output ',,,'];
  end
  
  if strcmp(cfg.conn.type, 'ft') && strcmp(cfg.conn.mvar, 'ft')
    output = [output ','];
  else
    output = [output sprintf('%1.f,', cfg.conn.order)];
  end
  
  output = [output ',,,,']; % leave some extras in case we add options
  %-----------------%
  
  load([cfg.log filesep 'connsum'], 'connsum')
  
  %---------------------------%
  %-loop over chan-chan
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

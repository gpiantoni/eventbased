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
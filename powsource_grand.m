function powsource_grand(cfg)
%POWSOURCE_GRAND grand pow source average

mversion = 6;
%06 12/02/14 save source stat for average
%05 12/02/03 renamed to powsource_grand
%04 12/01/31 plot image with most significant voxels
%03 12/01/11 find peaks with pos and areas for sources
%02 12/01/10 similar to grandpow, it reports clusters
%01 12/01/10 created from grandpow

%---------------------------%
%-start log
output = sprintf('%s (v%02.f) started at %s on %s\n', ...
  mfilename,  mversion, datestr(now, 'HH:MM:SS'), datestr(now, 'dd-mmm-yy'));
tic_t = tic;
%---------------------------%

%---------------------------%
%-loop over conditions
for e = 1:numel(cfg.poweffect)
  k = cfg.erpeffect(e);
  
  %-----------------%
  %-file for each cond
  condname = regexprep(cfg.test{k}, '*', '');
  outputfile = sprintf('powsource_*_%s.mat', condname);
  
  allsub = dir([cfg.dpow outputfile]);
  %-----------------%
  
  %-----------------%
  %-read data
  for s = 1:numel(allsub)
    load([cfg.dpow allsub(s).name]);
    spre(s,:) = souPre;
    sall(s,:) = source;
    clear source
  end
  %-----------------%
  
  %-----------------%
  %-powsource over subj: loop over areas
  for a = 1:size(sall,2) % this is powpeak, but implicit
    cfg1 = [];
    cfg1.keepindividual = 'yes';
    gpowsouPre{k,a} = ft_sourcegrandaverage(cfg1, spre{:,a});
    % gpowsouPre{k,a}.avg.tscore =  gpowsouPre{k,a}.avg.pow ./ (sqrt(gpowsouPre{k,a}.var.pow) / sqrt(numel(cfg.subjall)));
    gpowsource{k,a} = ft_sourcegrandaverage(cfg1, sall{:,a});
    % gpowsource{k,a}.avg.tscore =  gpowsource{k,a}.avg.pow ./ (sqrt(gpowsource{k,a}.var.pow) / sqrt(numel(cfg.subjall)));
  end
  %-----------------%
  
  clear sall spre
end
%---------------------------%

%---------------------------%
%-statistics for main effects
%-----------------%
%-use predefined or power-peaks for areas of interest
if strcmp(cfg.powsource.areas, 'manual')
  powpeak = cfg.powsource.powpeak;
elseif strcmp(cfg.powsource.areas, 'powpeak')
  load([cfg.dpow cfg.proj '_powpeak'], 'powpeak')
end
%-----------------%

soupeak = [];
for p = 1:numel(powpeak)
  output = sprintf('%s\n%s:\n', output, powpeak(p).name);
  
  h = figure;
  [soupos powstat{p} outtmp] = reportsource(gpowsource{cfg.poweffect, p}, gpowsouPre{cfg.poweffect, p});
  soupeak(p).pos = soupos;
  soupeak(p).center = mean(soupos,1);
  soupeak(p).name = powpeak(p).name;
  output = [output outtmp];
  
  %--------%
  pngname = sprintf('gpowpeak_%1.f_%s', cfg.poweffect, powpeak(p).name);
  saveas(gcf, [cfg.log filesep pngname '.png'])
  close(gcf); drawnow
  
  [~, logfile] = fileparts(cfg.log);
  system(['ln ' cfg.log filesep pngname '.png ' cfg.rslt pngname '_' logfile '.png']);
  %--------%

end

save([cfg.dpow cfg.proj '_soupeak'], 'soupeak')

%-----------------%
%-save
for p = 1:numel(powstat)
  powstat{p} = rmfield(powstat{p}, ...
    {'mask', 'ref', ...
    'posclusters', 'posclusterslabelmat', 'posdistribution', ...
    'negclusters', 'negclusterslabelmat', 'negdistribution'});
  powstat{p}.cfg = []; % this is huge
end
save([cfg.dpow cfg.proj '_grandpowsource'], 'powstat', '-v7.3')
%-----------------%
%---------------------------%

%---------------------------%
%-end log
toc_t = toc(tic_t);
outtmp = sprintf('%s (v%02.f) ended at %s on %s after %s\n\n', ...
  mfilename, mversion, datestr(now, 'HH:MM:SS'), datestr(now, 'dd-mmm-yy'), ...
  datestr( datenum(0, 0, 0, 0, 0, toc_t), 'HH:MM:SS'));
output = [output outtmp];

%-----------------%
fprintf(output)
fid = fopen([cfg.log '.txt'], 'a');
fwrite(fid, output);
fclose(fid);
%-----------------%
%---------------------------%
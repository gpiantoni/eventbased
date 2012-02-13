function powsource_grand(cfg)
%POWSOURCE_GRAND grand pow source average

mversion = 5;
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

%-----------------%
%-save (it's too big)
[s1 s2] = size(gpowsource);

for i1 = 1:s1
  for i2 = 1:s2
    souPre{i1,i2} = rmfield(gpowsouPre{i1,i2}, 'trial');
    source{i1,i2} = rmfield(gpowsource{i1,i2}, 'trial');
  end
end
save([cfg.dpow cfg.proj '_grandpowsource'], 'source', 'souPre', '-v7.3')
%-----------------%
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
  [soupos outtmp] = reportsource(gpowsource{cfg.poweffect, p}, gpowsouPre{cfg.poweffect, p});
  soupeak(p).pos = soupos;
  soupeak(p).center = mean(soupos,1);
  soupeak(p).name = powpeak(p).name;
  output = [output outtmp];  
  
  saveas(h, [cfg.log filesep 'pow_' powpeak(p).name '.png'])
  close(h)
end

save([cfg.dpow cfg.proj '_soupeak'], 'soupeak')
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
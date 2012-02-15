function erpsource_grand(cfg)
%ERPSOURCE_GRAND grand erp source average

mversion = 5;
%05 12/02/14 save source stat for average
%04 12/02/03 renamed to erpsource_grand
%03 12/01/31 plot image with most significant voxels
%02 12/01/12 compare against baseline
%01 12/01/11 created from grandpowsource

%---------------------------%
%-start log
output = sprintf('%s (v%02.f) started at %s on %s\n', ...
  mfilename,  mversion, datestr(now, 'HH:MM:SS'), datestr(now, 'dd-mmm-yy'));
tic_t = tic;
%---------------------------%

%---------------------------%
%-loop over conditions
for e = 1:numel(cfg.erpeffect)
  k = cfg.erpeffect(e);
  
  %-----------------%
  %-file for each cond
  condname = regexprep(cfg.test{k}, '*', '');
  outputfile = sprintf('erpsource_*_%s.mat', condname);
  
  allsub = dir([cfg.derp outputfile]);
  %-----------------%
  
  %-----------------%
  %-read data
  for s = 1:numel(allsub)
    load([cfg.derp allsub(s).name]);
    spre(s,:) = souPre;
    sall(s,:) = source;
    clear source souPre
  end
  %-----------------%
  
  %-----------------%
  %-powsource over subj: loop over areas
  for a = 1:size(sall,2) % this is erppeak, but implicit
    cfg1 = [];
    cfg1.keepindividual = 'yes';
    cfg1.parameter = 'pow'; % instead of nai
    gerpsouPre{k,a} = ft_sourcegrandaverage(cfg1, spre{:,a});
    gerpsource{k,a} = ft_sourcegrandaverage(cfg1, sall{:,a});
    % gerpsource{k,a}.avg.tscore =  gerpsource{k,a}.avg.nai ./ (sqrt(gerpsource{k,a}.var.nai) / sqrt(numel(cfg.subjall))); % nai in lcmv instead of avg in pow
  end
  %-----------------%
  
  clear sall spre
end
%---------------------------%

%---------------------------%
%-statistics for main effects
%-----------------%
%-use predefined or erp-peaks for areas of interest
if strcmp(cfg.erpsource.areas, 'manual')
  erppeak = cfg.erpsource.erppeak;
elseif strcmp(cfg.erpsource.areas, 'erppeak')
  load([cfg.derp cfg.proj '_erppeak'], 'erppeak')
end
%-----------------%

soupeak = [];
for p = 1:numel(erppeak)
  output = sprintf('%s\n%s:\n', output, erppeak(p).name);
  h = figure;
  [soupos erpstat{p} outtmp] = reportsource(gerpsource{cfg.erpeffect, p}, gerpsouPre{cfg.erpeffect,p});
  soupeak(p).pos = soupos;
  soupeak(p).center = mean(soupos,1);
  soupeak(p).name = erppeak(p).name;
  output = [output outtmp];  
  
  %--------%
  pngname = sprintf('gerppeak_%1.f_%s', cfg.poweffect, powpeak(p).name);
  saveas(gcf, [cfg.log filesep pngname '.png'])
  close(gcf); drawnow
  
  [~, logfile] = fileparts(cfg.log);
  system(['ln ' cfg.log filesep pngname '.png ' cfg.rslt pngname '_' logfile '.png']);
  %--------%
      
end

save([cfg.derp cfg.proj '_soupeak'], 'soupeak')

%-----------------%
%-save
for p = 1:numel(erpstat)
  erpstat{p} = rmfield(erpstat{p}, ...
    {'mask', 'ref', ...
    'posclusters', 'posclusterslabelmat', 'posdistribution', ...
    'negclusters', 'negclusterslabelmat', 'negdistribution'});
  erpstat{p}.cfg = []; % this is huge
end
save([cfg.derp cfg.proj '_granderpsource'], 'erpstat', '-v7.3')
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
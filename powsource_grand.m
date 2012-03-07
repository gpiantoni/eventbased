function powsource_grand(cfg)
%POWSOURCE_GRAND group-analysis of POW source data
%
% CFG
%  .cond: name to be used to save powsource_PROJNAME and figures
%  .test: a cell with the condition defined by redef. 
%
%  .dpow: directory to save POW data
%  .poweffect: effect of interest to create powpeak. If empty, no stats.
%
% Options from reportsource:
%  .powsource.clusterstatistics: 'maxsize' or 'max'
%  .powsource.clusteralpha: level to select sensors (default 0.05)
%  .powsource.maxvox: max number of significant voxels to be used in soupeak
%
% OUT
%  [cfg.dpow 'COND_grandpowsource']: source analysis for all subject
%  [cfg.dpow 'COND_soupeak']: significant source peaks in the POW
%
% FIGURES
%  gpowpeak_POWEFFECT_POWPEAKNAME: 3d plot of the source for one peak
%
% Part of EVENTBASED group-analysis
% see also ERP_SUBJ, ERP_GRAND, ERPSOURCE_SUBJ, ERPSOURCE_GRAND, 
% POW_SUBJ, POW_GRAND, POWSOURCE_SUBJ, POWSOURCE_GRAND, 
% POWCORR_SUBJ, POWCORR_SUBJ,
% CONN_SUBJ, CONN_GRAND, CONN_STAT

%---------------------------%
%-start log
output = sprintf('%s started at %s on %s\n', ...
  mfilename,  datestr(now, 'HH:MM:SS'), datestr(now, 'dd-mmm-yy'));
tic_t = tic;
%---------------------------%

%---------------------------%
%-loop over conditions
for e = 1:numel(cfg.poweffect)
  k = cfg.poweffect(e);
  
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
    gpowsource{k,a} = ft_sourcegrandaverage(cfg1, sall{:,a});
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
  [soupos powstat{p} outtmp] = reportsource(cfg.powsource, gpowsource{cfg.poweffect, p}, gpowsouPre{cfg.poweffect, p});
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
  powstat{p}.cfg = []; % this is huge
end
save([cfg.dpow cfg.proj '_grandpowsource'], 'powstat', '-v7.3')
%-----------------%
%---------------------------%

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
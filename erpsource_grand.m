function erpsource_grand(cfg)
%ERPSOURCE_GRAND group-analysis of ERP source data
%
% CFG
%  .cond: name to be used to save erpsource_PROJNAME and figures
%  .test: a cell with the condition defined by redef. 
%  .log: name of the file and directory with analysis log
%  .rslt: directory images are saved into
%
%  .derp: directory to save ERP data
%  .erpeffect: effect of interest to create erppeak. If empty, no stats.
%
% Options from reportsource:
%  .erpsource.clusterstatistics: 'maxsize' or 'max'
%  .erpsource.clusteralpha: level to select sensors (default 0.05)
%  .erpsource.maxvox: max number of significant voxels to be used in soupeak
%
% OUT
%  [cfg.derp 'COND_granderpsource']: source analysis for all subject
%  [cfg.derp 'COND_soupeak']: significant source peaks in the ERP
%
% FIGURES
%  gerppeak_ERPEFFECT_ERPPEAKNAME: 3d plot of the source for one peak
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
  end
  %-----------------%
  
  clear sall spre
end
%---------------------------%

%---------------------------%
%-statistics for main effects
%-----------------%
%-use predefined or erppeaks for areas of interest
if strcmp(cfg.erpsource.areas, 'manual')
  erppeak = cfg.erpsource.erppeak;
elseif strcmp(cfg.erpsource.areas, 'erppeak')
  load([cfg.derp cfg.cond '_erppeak'], 'erppeak')
end
%-----------------%

soupeak = [];
for p = 1:numel(erppeak)
  output = sprintf('%s\n%s:\n', output, erppeak(p).name);
  h = figure;
  [soupos erpstat{p} outtmp] = reportsource(cfg.erpsource, gerpsource{cfg.erpeffect, p}, gerpsouPre{cfg.erpeffect,p});
  soupeak(p).pos = soupos;
  soupeak(p).center = mean(soupos,1);
  soupeak(p).name = erppeak(p).name;
  output = [output outtmp];  
  
  %--------%
  pngname = sprintf('gerppeak_%1.f_%s', cfg.erpeffect, erppeak(p).name);
  saveas(gcf, [cfg.log filesep pngname '.png'])
  close(gcf); drawnow
  
  [~, logfile] = fileparts(cfg.log);
  system(['ln ' cfg.log filesep pngname '.png ' cfg.rslt pngname '_' logfile '.png']);
  %--------%
      
end

save([cfg.derp cfg.cond '_soupeak'], 'soupeak')

%-----------------%
%-save
for p = 1:numel(erpstat)
  erpstat{p}.cfg = []; % this is huge
end
save([cfg.derp cfg.cond '_granderpsource'], 'erpstat', '-v7.3')
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
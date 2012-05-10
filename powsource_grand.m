function powsource_grand(cfg)
%POWSOURCE_GRAND group-analysis of POW source data
%
% CFG
%  .cond: name to be used to save powsource_PROJNAME and figures
%  .test: a cell with the condition defined by redef. 
%  .log: name of the file and directory with analysis log
%  .rslt: directory images are saved into
%
%  .dpow: directory to save POW data
%  .poweffect: index of interest to create powpeak, can be a row vector,
%              but it only uses the first one
%
% Options from reportsource:
%  .powsource.clusterstatistics: 'maxsize' or 'max'
%  .powsource.clusteralpha: level to select sensors (default 0.05)
%                           it can be a string in format '5%' to take top 5 voxels and put them in a cluster.
%  .powsource.maxvox: max number of significant voxels to be used in soupeak
%
% Options if you want to create significance mask
%  .powsource.nifti: directory and initial part of the name where you want to save the masks
%  .dti.ref: template for mask ('/usr/share/data/fsl-mni152-templates/MNI152_T1_1mm_brain.nii.gz')
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
% POWCORR_SUBJ, POWCORR_GRAND,
% CONN_SUBJ, CONN_GRAND, CONN_STAT

%---------------------------%
%-start log
output = sprintf('%s started at %s on %s\n', ...
  mfilename,  datestr(now, 'HH:MM:SS'), datestr(now, 'dd-mmm-yy'));
tic_t = tic;
%---------------------------%

%---------------------------%
%-loop over conditions
for k = 1:numel(cfg.powsource.cond) % DOC: CFG.POWSOURCE.COND
  cond     = cfg.powsource.cond{k};
  condname = regexprep(cond, '*', '');
  
  %-----------------%
  %-file for each cond
  subjfile = @(s) sprintf('%spowsource_%02.f_%s.mat', cfg.dpow, s, condname);
  allname = cellfun(subjfile, num2cell(cfg.subjall), 'uni', 0);
  
  allfiles = true(1, numel(allname));
  for i = 1:numel(allname)
    if ~exist(allname{i}, 'file')
      output = [output sprintf('%s does not exist\n', allname{i})];
      allfiles(i) = false;
    end
  end
  allname = allname(allfiles);
  %-----------------%
  
  %-----------------%
  %-read data
  for s = 1:numel(allname)
    load(allname{s});
    spre(s,:) = souPre;
    sall(s,:) = source;
    clear source souPre
  end
  %-----------------%
  
  %-----------------%
  %-powsource over subj: loop over areas
  for a = 1:size(sall,2) % this is powpeak, but implicit
    cfg1 = [];
    cfg1.keepindividual = 'yes';
    gpowsouPre{p,a} = ft_sourcegrandaverage(cfg1, spre{:,a});
    gpowsource{p,a} = ft_sourcegrandaverage(cfg1, sall{:,a});
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
  peakname = regexprep(cfg.pow.refcond, '*', ''); % DOC: CFG.POW.REFCOND
  load([cfg.dpow cfg.cond peakname '_powpeak'], 'powpeak')
  
elseif strcmp(cfg.powsource.areas, 'powcorrpeak')
  peakname = regexprep(cfg.powcorr.refcond, '*', ''); % DOC: CFG.POWCORR.REFCOND
  load([cfg.dpow cfg.cond peakname '_powcorrpeak'], 'powcorrpeak')
  powpeak = powcorrpeak;
  
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
  
  %--------%
  if isfield(cfg.powsource, 'nifti') && ~isempty(cfg.powsource.nifti)
    
    dtimri = ft_read_mri(cfg.dti.ref);
    
    cfg1 = [];
    cfg1.parameter = 'image';
    souinterp = ft_sourceinterpolate(cfg1, powstat{p}, dtimri);
    
    cfg1 = [];
    cfg1.parameter = 'image';
    cfg1.filename = [cfg.powsource.nifti soupeak(p).name];
    ft_sourcewrite(cfg1, souinterp);
    gzip([cfg.powsource.nifti soupeak(p).name '.nii'])
    delete([cfg.powsource.nifti soupeak(p).name '.nii'])
  end
  %--------%

end

save([cfg.dpow cfg.cond '_soupeak'], 'soupeak')

%-----------------%
%-save
for p = 1:numel(powstat)
  powstat{p}.cfg = []; % this is huge
end
save([cfg.dpow cfg.cond '_grandpowsource'], 'powstat', '-v7.3')
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
function source_subj(cfg, subj)
%SOURCE_SUBJ create virtual electrode, to be used for connectivity analysis
%
% CFG
%  .data: path of /data1/projects/PROJ/subjects/
%  .rec: REC in /data1/projects/PROJ/recordings/REC/
%  .nick: NICK in /data1/projects/PROJ/subjects/0001/MOD/NICK/
%  .mod: modality, MOD in /data1/projects/PROJ/subjects/0001/MOD/NICK/
%  .endname: includes preprocessing steps (e.g. '_seldata_gclean_redef')
%
%  .log: name of the file and directory to save log
%  .dsou: directory for connectivity data
%  .sou.cond: cell with conditions (e.g. {'*cond1' '*cond2'})'
%
%-ROI parameters
%  .source.areas: 'channel', 'erp', 'dip', 'erppeak' or 'powpeak'
%
%    if 'channel'
%      .source.chan: a struct with
%        .name: 'name of group elec'
%        .chan: cell with electrode labels for each group
%
%    if 'erp'
%      .source.dip: a struct with
%        .name: 'name of dipole'
%        .time: time window of the ERP activity of interest (two scalars)
%      .derp: directory with ERP data
%      .source.refcond: condition with ERP used for reference topography (string)
%
%    if 'dip' (use dipoles of interest):
%      .source.beamformer: 'erp' or 'pow' (time-domain, lcmv, or freq-domain, dics)
%                        you need to have run 'erpsource_subj' or 'powsource_subj' respectively
%                        don't forget to keepfilter
%      .source.refcond: string of the condition used for source location
%      .source.fixedmom: logical (use the same moment for source or change it every time)
%      .source.dip(1).name: name of the dipole
%      .source.dip(1).pos: position of the dipole (X x 3, it'll do PCA on it)
%      (this function only works with beamforming, NOT with simple inversion of leadfield)
%
%    if 'erppeak' (use beamformer to construct virtual electrode):
%      .derp: directory with ERP data (you need 'erpsource_subj' with keepfilder)
%      .source.refcond: string of the condition used for source location
%      .source.fixedmom: logical (use the same moment for source or change it every time)
%
%    if 'powpeak' (use beamformer to construct virtual electrode):
%      Because this function returns the activity at the source in
%      time-domain, we need to use the beamformer in time-domain, i.e. LCMV
%      .dpow: directory with POW data
%      .derp: you need 'erpsource_subj' with keepfilder
%      .source.refcond: string of the condition used for source location
%      .source.fixedmom: logical (use the same moment for source or change it every time)
%
% IN:
%  data in /PROJ/subjects/SUBJ/MOD/NICK/
%  if .source.areas == 'erppeak' or ('dip' and .source.beamformer == 'erp')
%     [cfg.derp 'erpsource_SUBJ_COND']: source data for period of interest for each subject
%     [cfg.derp 'NICK_COND_soupeak']: significant source peaks in the ERP
%  if .source.areas == 'powpeak'  or ('dip' and .source.beamformer == 'pow')
%     [cfg.derp 'erpsource_SUBJ_COND']: source data for period of interest for each subject
%     [cfg.dpow 'NICK_COND_soupeak']: significant source peaks in the POW
%
% OUT
%  [cfg.dsou 'conn_SUBJ_COND']: virtual electrode for each channel
%
% Part of EVENTBASED single-subject
% see also ERP_SUBJ, ERP_GRAND,
% ERPSOURCE_SUBJ, ERPSOURCE_GRAND, ERPSTAT_SUBJ, ERPSTAT_GRAND,
% POW_SUBJ, POW_GRAND, POWCORR_SUBJ, POWCORR_GRAND,
% POWSOURCE_SUBJ, POWSOURCE_GRAND, POWSTAT_SUBJ, POWSTAT_GRAND,
% SOURCE_SUBJ, CONN_SUBJ, CONN_GRAND, CONN_STAT

%---------------------------%
%-start log
output = sprintf('%s (%04d) began at %s on %s\n', ...
  mfilename, subj, datestr(now, 'HH:MM:SS'), datestr(now, 'dd-mmm-yy'));
tic_t = tic;
%---------------------------%

%---------------------------%
%-prepare montage
switch cfg.source.areas
  
  case 'channel'
    [mont outtmp] = prepare_montage(cfg);
    
  case 'erp'
    condname = regexprep(cfg.source.refcond, '*', '');
    load([cfg.derp 'erp_' condname], 'erp')
    
    [mont outtmp] = prepare_montage(cfg, erp);
    
  case 'dip'
    condname = regexprep(cfg.source.refcond, '*', '');
    sourcename = sprintf('%ssource_s_A', cfg.source.beamformer);
    load(sprintf('%s%ssource_%04d_%s', cfg.derp, cfg.source.beamformer, subj, condname), sourcename) % source of interest
    
    [mont outtmp] = prepare_montage(cfg, eval(sourcename), cfg.source.dip);
    
  case 'erppeak'
    condname = regexprep(cfg.source.refcond, '*', '');
    load(sprintf('%serpsource_%04d_%s', cfg.derp, subj, condname), 'erpsource_s_A') % source of interest
    load(sprintf('%serpsource_peak_%s', cfg.derp, condname), 'erpsource_peak') % peaks in ERP
    
    [mont outtmp] = prepare_montage(cfg, erpsource_s_A, erpsource_peak);
    
  case 'powpeak'
    condname = regexprep(cfg.source.refcond, '*', '');
    load(sprintf('%serpsource_%04d_%s', cfg.derp, subj, condname), 'erpsource_s_A') % source of interest
    load(sprintf('%spowsource_peak_%s', cfg.dpow, condname), 'powsource_peak') % peaks in POW
    
    [mont outtmp] = prepare_montage(cfg, erpsource_s_A, powsource_peak);
    
end
output = [output outtmp];
%---------------------------%

%-------------------------------------%
%-loop over conditions
for k = 1:numel(cfg.source.cond)
  cond     = cfg.source.cond{k};
  condname = regexprep(cond, '*', '');
  
  %---------------------------%
  %-read data
  [data badchan] = load_data(cfg, subj, cond);
  if isempty(data)
    output = sprintf('%sCould not find any file for condition %s\n', ...
      output, cond);
    continue
  end
  
  outputfile = sprintf('source_%04d_%s', subj, condname);
  %---------------------------%
  
  %---------------------------%
  %-apply montage (if using two-step procedure)
  
  data = ft_apply_montage(data, mont, 'feedback', 'none');
  
  if strcmp(cfg.source.areas, 'dip')
    data = pcadata(data, cfg.source.dip, cfg.source.fixedmom);
  end
  if strcmp(cfg.source.areas, 'erppeak')
    data = pcadata(data, erpsource_peak, cfg.source.fixedmom);
  end
  if strcmp(cfg.source.areas, 'powpeak')
    data = pcadata(data, powsource_peak, cfg.source.fixedmom);
  end
  
  source = ft_checkdata(data, 'hassampleinfo', 'yes'); % recreate sampleinfo which got lost (necessary for selfromraw)
  %---------------------------%
  
  save([cfg.dsou outputfile], 'source')
  
end
%-------------------------------------%

%---------------------------%
%-end log
toc_t = toc(tic_t);
outtmp = sprintf('%s (%04d) ended at %s on %s after %s\n\n', ...
  mfilename, subj, datestr(now, 'HH:MM:SS'), datestr(now, 'dd-mmm-yy'), ...
  datestr( datenum(0, 0, 0, 0, 0, toc_t), 'HH:MM:SS'));
output = [output outtmp];

%-----------------%
fprintf(output)
fid = fopen([cfg.log '.txt'], 'a');
fwrite(fid, output);
fclose(fid);
%-----------------%
%---------------------------%

%---------------------------------------------------------%
%-SUBFUNCTION: pcadata
%---------------------------------------------------------%
function [data] = pcadata(data, soupeak, fixedmom)
%PCADATA simplify data using pca on each region of interest
% keep only the first component

%-------------------------------------%
%-loop over regions of interest
trial = [];
for i1 = 1:numel(soupeak)
  
  %-----------------%
  %-find channels belonging to region of interest
  newname = sprintf('%s', soupeak(i1).name);
  label{i1,1} = newname;
  ivox = ~cellfun(@isempty, strfind(data.label, newname));
  %-----------------%
  
  if ~fixedmom
    %-----------------%
    %-the moment is different in each trial
    for t = 1:numel(data.trial)
      [~, pcatrl] = princomp(data.trial{t}(ivox,:)');
      trial{t}(i1,:) = pcatrl(:,1)';
    end
    %-----------------%
    
  else
    
    %-----------------%
    %-get coefficient for all the trials, then apply it to the data
    alltrl = [data.trial{:}];
    [coef] = princomp(alltrl(ivox,:)');
    
    for t = 1:numel(data.trial)
      [pcatrl] = data.trial{t}(ivox,:)' * coef;
      trial{t}(i1,:) = pcatrl(:,1)';
    end
    %-----------------%
  end
  
end
%-------------------------------------%

data.trial = trial;
data.label = label;
%---------------------------------------------------------%

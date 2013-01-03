function source_subj(info, opt, subj)
%SOURCE_SUBJ create virtual electrode, to be used for connectivity analysis
%
% INFO

% CFG.OPT
%  .roi: 'channel', 'erp', 'dip', 'erppeak' or 'powpeak'
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

% Check that when you use freesurfer, you need to check which dipoles are
% actually used
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

% IN:
%  data in /PROJ/subjects/SUBJ/MOD/NICK/
%  if .source.areas == 'erppeak' or ('dip' and .source.beamformer == 'erp')
%     [info.derp 'erpsource_SUBJ_COND']: source data for period of interest for each subject
%     [info.derp 'NICK_COND_soupeak']: significant source peaks in the ERP
%  if .source.areas == 'powpeak'  or ('dip' and .source.beamformer == 'pow')
%     [info.derp 'erpsource_SUBJ_COND']: source data for period of interest for each subject
%     [info.dpow 'NICK_COND_soupeak']: significant source peaks in the POW
%
% OUT
%  [cfg.dsou 'conn_SUBJ_COND']: virtual electrode for each channel
%
% Part of EVENTBASED single-subject
% see also ERP_SUBJ, ERP_GRAND,
% ERPSOURCE_SUBJ, ERPSOURCE_GRAND, ERPSTAT_SUBJ, ERPSTAT_GRAND,
% POW_SUBJ, POW_GRAND, POW_GRP, POWCORR_SUBJ, POWCORR_GRAND,
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
switch opt.roi 
  
  case 'channel'
    [mont outtmp] = prepare_montage(opt);
    
  case 'erp'
    condname = regexprep(cfg.source.refcond, '*', '');
    load([info.derp 'erp_' condname], 'erp')
    
    [mont outtmp] = prepare_montage(cfg, erp);
    
  case 'dip'
    condname = regexprep(cfg.source.refcond, '*', '');
    sourcename = sprintf('%ssource_s_A', cfg.source.beamformer);
    load(sprintf('%s%ssource_%04d_%s', info.derp, cfg.source.beamformer, subj, condname), sourcename) % source of interest
    
    [mont outtmp] = prepare_montage(cfg, eval(sourcename), cfg.source.dip);
    
  case 'erppeak'
    condname = regexprep(cfg.source.refcond, '*', '');
    load(sprintf('%serpsource_%04d_%s', info.derp, subj, condname), 'erpsource_s_A') % source of interest
    load(sprintf('%serpsource_peak_%s', info.derp, condname), 'erpsource_peak') % peaks in ERP
    
    [mont outtmp] = prepare_montage(cfg, erpsource_s_A, erpsource_peak);
    
  case 'powpeak'
    condname = regexprep(cfg.source.refcond, '*', '');
    load(sprintf('%serpsource_%04d_%s', info.derp, subj, condname), 'erpsource_s_A') % source of interest
    load(sprintf('%spowsource_peak_%s', info.dpow, condname), 'powsource_peak') % peaks in POW
    
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
fid = fopen([info.log '.txt'], 'a');
fwrite(fid, output);
fclose(fid);
%-----------------%
%---------------------------%

%---------------------------------------------------------%
%-SUBFUNCTION FOR MONTAGE---------------------------------%
%---------------------------------------------------------%
%-----------------------------------------------%
%-PREPARE_MONTAGE_CHANNEL-----------------------%
%-----------------------------------------------%
function [mont output] = prepare_montage_channel(cfg)

%-----------------%
%-rename
grpchan = cfg.source.chan;
label = cfg.seldata.label;
%-----------------%

%-----------------%
%-prepare TRA
nnew = numel(grpchan);
nold = numel(cfg.seldata.label);
tra = zeros(nnew, nold);

output = sprintf('%d sensor groups:', nnew);

for i = 1:nnew
  ism = ismember(label, grpchan(i).chan)';
  tra(i,:) = ism / numel(find(ism));
  output = [output sprintf(' %s (%d elec),', grpchan(i).name, numel(find(ism)))];
end
%-----------------%

%-----------------%
%-prepare MONT
mont.labelorg = label;
mont.labelnew = {grpchan.name}';
mont.tra = tra;
%-----------------%
%-----------------------------------------------%


%-----------------------------------------------%
%-PREPARE_MONTAGE_PEAK--------------------------%
%-----------------------------------------------%
function [mont output] = prepare_montage_peak(source, peak)

%---------------------------%
%-check input
if numel(source) == 1 % called from cfg.source.areas = 'dip'
  source = repmat(source, size(peak));
end

if numel(source) ~= numel(peak)
  error(sprintf('the number of sources (%1d) should be identical to the number of significant peaks (%1d)', numel(source), numel(peak)))
end

output = '';
%---------------------------%

%-------------------------------------%
%-loop over sources
%-----------------%
%-alloc mont
nvox = size(cat(1,peak.pos),1) * 3;
nchan = size(source{1}.avg.filter{source{1}.inside(1)},2);
tra = NaN(nvox, nchan);
%-----------------%

cnt = 1;
for i = 1:numel(source)
  %-----------------%
  %-only inside dipole
  sou = source{i};
  sou.pos = source{i}.pos(source{i}.inside,:);
  sou.avg.filter = source{i}.avg.filter(source{i}.inside);
  %-----------------%
  
  [~, isou, ipeak] = intersect(sou.pos, peak(i).pos, 'rows');
  
  %-----------------%
  %-output
  maxpos = max(peak(i).pos(ipeak,:));
  meanpos = mean(peak(i).pos(ipeak,:));
  minpos = min(peak(i).pos(ipeak,:));
  
  outtmp = sprintf(['%s (defined at% 5d locations) has% 5d dipoles:\n' ...
  '                         x = [% 6.1f % 6.1f % 6.1f]\n', ...
  '                         y = [% 6.1f % 6.1f % 6.1f]\n', ...
  '                         z = [% 6.1f % 6.1f % 6.1f]\n'], ...
  peak(i).name, size(peak(i).pos,1), size(ipeak,1), ...
  maxpos(1), meanpos(1), minpos(1), ...
  maxpos(2), meanpos(2), minpos(2), ...
  maxpos(3), meanpos(3), minpos(3));
  output = [output outtmp];
  %-----------------%
  
  %-----------------%
  %-are all the voxels in the brain
  if size(ipeak,1) ~= size(peak(i).pos,1)
    outtmp = sprintf('warning: source %1d has % 3d voxels in the brain out of % 3d\n', i, size(ipeak,1), size(peak(i).pos,1));
    output = [output outtmp];
  end
  %-----------------%
  
  %-----------------%
  %-per voxel
  for v = 1:numel(isou)
    tra(cnt + (0:2),:) = sou.avg.filter{isou(v)};
    
    %-------%
    %-per moment
    labelnew{cnt,1} = sprintf('%s_%04d_a', peak(i).name, v);
    labelnew{cnt+1,1} = sprintf('%s_%04d_b', peak(i).name, v);
    labelnew{cnt+2,1} = sprintf('%s_%04d_c', peak(i).name, v);
    %-------%
    
    cnt = cnt + 3;
  end
  %-----------------%
  
end

%-----------------%
%-clean up the tra
tra = tra(1:cnt-1,:);
%-----------------%
%-------------------------------------%

mont.labelorg = source{1}.cfg.channel;
mont.labelnew = labelnew;
mont.tra = tra;
%-----------------------------------------------%
%---------------------------------------------------------%

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

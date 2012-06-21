function [mont output] = prepare_montage(cfg, data, peak)
%PREPARE_MONTAGE create montage as spatial filters to be used for sources
%
% .conn.areas
%-CHANNEL: average over channels
%   .conn.chan: a struct with
%      .name: 'name of group elec'
%      .chan: cell with electrode labels for each group
%   .seldata.label: labels of electrodes in the data
%
%-ERP: use averaged topography
%   .conn.dip: a struct with
%      .name: 'name of dipole'
%      .time: time window of the ERP activity of interest (two scalars)
%   erp: data with erp for condition of interest
%
%-ERPPEAK, POWPEAK: convert from beamforming source into montage
%   source: subject-specific source, you should have kept the real filters
%   soupeak: peaks for the condition of interest
%
% This function returns "mont" which can be used with ft_apply_montage
%
% TODO: ICA to create spatial filters

switch cfg.conn.areas
  
  case 'channel'
    [mont output] = prepare_montage_channel(cfg);
    
  case 'erp'
    [mont output] = prepare_montage_erp(cfg, data);
    
  case {'erppeak' 'powpeak'}
    [mont output] = prepare_montage_peak(data, peak);
 
end

%-----------------------------------------------%
%-PREPARE_MONTAGE_CHANNEL-----------------------%
%-----------------------------------------------%
function [mont output] = prepare_montage_channel(cfg)

%-----------------%
%-rename
grpchan = cfg.conn.chan;
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
%-PREPARE_MONTAGE_ERP---------------------------%
%-----------------------------------------------%
function [mont output] = prepare_montage_erp(cfg, erp)

%-----------------%
%-rename
dip = cfg.conn.dip;
%-----------------%

%-----------------%
%-prepare TRA
nnew = numel(dip);
nold = numel(erp.label);
tra = zeros(nnew, nold);

output = sprintf('%d sensor groups:', nnew);

for i = 1:numel(dip)
  dipavg = ft_selectdata(erp, 'toilim', dip(i).time, 'avgovertime', 'yes');
  tra(i, :) = dipavg.avg';
end
%-----------------%

%-----------------%
%-prepare MONT
mont.labelorg = erp.label;
mont.labelnew = {dip.name}';
mont.tra = tra;
%-----------------%

%-----------------------------------------------%
%-PREPARE_MONTAGE_PEAK--------------------------%
%-----------------------------------------------%
function [mont output] = prepare_montage_peak(source, peak)

%---------------------------%
%-check input
if numel(source) ~= numel(peak)
  error(sprintf('the number of sources (%1.f) should be identical to the number of significant peaks (%1.f)', numel(source), numel(peak)))
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
  [~, isou, ipeak] = intersect(source{i}.pos, peak(i).pos, 'rows');
  
  %-----------------%
  %-are all the voxels in the brain
  if size(ipeak,1) ~= size(peak(i).pos,1)
    outtmp = sprintf('warning: source %1.f has % 3.f voxels in the brain out of % 3.f\n', i, size(ipeak,1), size(peak(i).pos,1));
    output = [output outtmp];
  end
  %-----------------%
  
  %-----------------%
  %-per voxel
  for v = 1:numel(isou)
    tra(cnt + (0:2),:) = source{i}.avg.filter{isou(v)};
    
    %-------%
    %-per moment
    labelnew{cnt,1} = sprintf('%s_%04.f_a', peak(i).name, v);
    labelnew{cnt+1,1} = sprintf('%s_%04.f_b', peak(i).name, v);
    labelnew{cnt+2,1} = sprintf('%s_%04.f_c', peak(i).name, v);
    %-------%
    
    cnt = cnt + 3;
  end
  %-----------------%
  
end
%-------------------------------------%

mont.labelorg = source{1}.cfg.channel;
mont.labelnew = labelnew;
mont.tra = tra;
%-----------------------------------------------%
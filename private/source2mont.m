function [mont output] = source2mont(source, soupeak)
%SOURCE2MONT convert from beamforming source into montage
% TODO: another way to do this (at least for lcmv) is to keep the mom (see
% gosderpsource) and then use that moment. So you only have one timeseries
% per voxel (you can run pca on that afterwards).
% At the moment, it does not keep the moment because it uses too much space

%02 12/02/11 get orig channel labels from source{1} (it need not use all channels)
%01 12/01/11 created

%---------------------------%
%-check input
if numel(source) ~= numel(soupeak)
  error(sprintf('the number of sources (%1.f) should be identical to the number of significant peaks (%1.f)', numel(source), numel(soupeak)))
end

output = '';
%---------------------------%

%-------------------------------------%
%-loop over sources
%-----------------%
%-alloc mont
nvox = size(cat(1,soupeak.pos),1) * 3;
nchan = size(source{1}.avg.filter{source{1}.inside(1)},2);
tra = NaN(nvox, nchan);
%-----------------%

cnt = 1;
for i = 1:numel(source)
  [~, isou, ipeak] = intersect(source{i}.pos, soupeak(i).pos, 'rows');
  
  %-----------------%
  %-are all the voxels in the brain
  if size(ipeak,1) ~= size(soupeak(i).pos,1)
    outtmp = sprintf('warning: source %1.f has % 3.f voxels in the brain out of % 3.f\n', i, size(ipeak,1), size(soupeak(i).pos,1));
    output = [output outtmp];
  end
  %-----------------%
  
  %-----------------%
  %-per voxel
  for v = 1:numel(isou)
    tra(cnt + (0:2),:) = source{i}.avg.filter{isou(v)};
    
    %-------%
    %-per moment
    labelnew{cnt,1} = sprintf('%s_%04.f_a', soupeak(i).name, v);
    labelnew{cnt+1,1} = sprintf('%s_%04.f_b', soupeak(i).name, v);
    labelnew{cnt+2,1} = sprintf('%s_%04.f_c', soupeak(i).name, v);
    %-------%
    
    cnt = cnt + 3;
  end
  %-----------------%
  
end
%-------------------------------------%

mont.labelorg = source{1}.cfg.channel;
mont.labelnew = labelnew;
mont.tra = tra;

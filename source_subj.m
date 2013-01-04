function source_subj(info, opt, subj)
%SOURCE_SUBJ create virtual electrode, to be used for connectivity analysis
%
% INFO
%  .log: name of the file and directory to save log
%  .dsou: directory with SOURCE data
%
% CFG.OPT
%  .cond: cell with conditions (e.g. {'*cond1' '*cond2'})
%  .mont.type: type of montage to compute source activity
%  if 'channel'
%    .mont.chan: a struct with
%      .name: 'name of group elec'
%      .chan: cell with electrode labels for each group
%    .label: labels of electrodes in the data
%
%  if 'roi'
%    .roi: see GET_ROI
%
%    .mont.forward: forward model
%    if 'leadfield'
%      TODO
%    if 'lcmv'
%      .mont.refcond: reference condition, in which you estimated
%                     erpsource_subj with keepfilter (as cell)
%    .fixedmom: logical (use the same moment for source at every trial)
%
% IN:
%  data in /PROJ/subjects/SUBJ/MOD/NICK/
%
% OUT
%  [cfg.dsou 'conn_SUBJ_COND'] 'source': virtual electrode for each channel
%
% Part of EVENTBASED group-analysis
% see also ERP_SUBJ, ERP_GRAND, 
% ERPSOURCE_SUBJ, ERPSOURCE_GRAND, ERPSTAT_SUBJ, ERPSTAT_GRAND,
% POW_SUBJ, POW_GRAND, POW_GRP, POWCORR_SUBJ, POWCORR_GRAND,
% POWSOURCE_SUBJ, POWSOURCE_GRAND, POWSTAT_SUBJ, POWSTAT_GRAND,
% SOURCE_SUBJ, CONN_SUBJ, CONN_GRAND, CONN_STAT,
% R_GRAND

%---------------------------%
%-start log
output = sprintf('%s (%04d) began at %s on %s\n', ...
  mfilename, subj, datestr(now, 'HH:MM:SS'), datestr(now, 'dd-mmm-yy'));
tic_t = tic;
%---------------------------%

if ~isfield(opt, 'fixedmom'); opt.fixedmom = false; end

%---------------------------%
%-prepare montage
switch opt.mont.type
  
  case 'channel'
    [mont outtmp] = prepare_montage_channel(opt);
        
  case 'roi'
    
    roi = get_roi(info, opt.mont.roi);
    
    switch opt.mont.forward
      case 'lcmv'
        info.subjall = subj;
        [~, fwd] = load_subj(info, 'erpsource', opt.mont.refcond{1});
        
    end
    [mont outtmp] = prepare_montage_roi(roi, fwd{1});
    
end
output = [output outtmp];
%---------------------------%

%-------------------------------------%
%-loop over conditions
for k = 1:numel(opt.cond)
  cond     = opt.cond{k};
  condname = regexprep(cond, '*', '');
  
  %---------------------------%
  %-read data
  [data] = load_data(info, subj, cond);
  if isempty(data)
    output = sprintf('%sCould not find any file for condition %s\n', ...
      output, cond);
    continue
  end
  
  outputfile = sprintf('source_%04d_%s', subj, condname);
  %---------------------------%
  
  %---------------------------%
  %-apply montage
  source = ft_apply_montage(data, mont, 'feedback', 'none');
  if strcmp(opt.mont.type, 'roi')
    source = pcadata(source, roi, opt.fixedmom);
  end
  %---------------------------%
  source.cfg = [];
  
  save([info.dsou outputfile], 'source')
  
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
grpchan = cfg.mont.chan;
label = cfg.label;
%-----------------%

%-----------------%
%-prepare TRA
nnew = numel(grpchan);
nold = numel(label);
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
function [mont output] = prepare_montage_roi(roi, fwd)

%-------------------------------------%
%-forward model
%---------------------------%
%-based on LCMV
%-----------------%
%-only inside dipole
pos = fwd.pos(fwd.inside,:);
filter = fwd.avg.filter(fwd.inside);
label = fwd.cfg.channel;
%-----------------%
%---------------------------%
%-------------------------------------%

%-------------------------------------%
%-loop over sources
%-----------------%
%-alloc mont
nvox = size(cat(1, roi.pos),1) * 3;
nchan = size(fwd.avg.filter{fwd.inside(1)},2);
tra = NaN(nvox, nchan);
%-----------------%

output = '';
cnt = 1;

for i = 1:numel(roi)
  
  [~, isou, iroi] = intersect(pos, roi(i).pos, 'rows');
  
  %-----------------%
  %-ERROR if regions of interest is not in the brain
  if isempty(isou)
    output = [output sprintf('ROI %s is not in the brain, check the position of the dipoles\n', roi(i).name)];
    continue
  end
  %-----------------%
  
  %-----------------%
  %-output
  maxpos = max(roi(i).pos(iroi,:), [], 1);
  meanpos = mean(roi(i).pos(iroi,:), 1);
  minpos = min(roi(i).pos(iroi,:), [], 1);
  
  outtmp = sprintf(['%s (defined at% 5d locations) has% 5d dipoles:\n' ...
  '                         x = [% 6.1f % 6.1f % 6.1f]\n', ...
  '                         y = [% 6.1f % 6.1f % 6.1f]\n', ...
  '                         z = [% 6.1f % 6.1f % 6.1f]\n'], ...
  roi(i).name, size(roi(i).pos,1), size(iroi,1), ...
  maxpos(1), meanpos(1), minpos(1), ...
  maxpos(2), meanpos(2), minpos(2), ...
  maxpos(3), meanpos(3), minpos(3));
  output = [output outtmp];
  %-----------------%
  
  %-----------------%
  %-are all the voxels in the brain
  if size(iroi,1) ~= size(roi(i).pos,1)
    outtmp = sprintf('warning: source %1d has % 3d voxels in the brain out of % 3d\n', i, size(iroi,1), size(roi(i).pos,1));
    output = [output outtmp];
  end
  %-----------------%
  
  %-----------------%
  %-per voxel
  for v = 1:numel(isou)
    tra(cnt + (0:2),:) = filter{isou(v)};
    
    %-------%
    %-per moment
    labelnew{cnt,1} = sprintf('%s_%04d_a', roi(i).name, v);
    labelnew{cnt+1,1} = sprintf('%s_%04d_b', roi(i).name, v);
    labelnew{cnt+2,1} = sprintf('%s_%04d_c', roi(i).name, v);
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

mont.labelorg = label;
mont.labelnew = labelnew;
mont.tra = tra;
%-----------------------------------------------%
%---------------------------------------------------------%

%---------------------------------------------------------%
%-SUBFUNCTION: pcadata
%---------------------------------------------------------%
function [data] = pcadata(data, roi, fixedmom)
%PCADATA simplify data using pca on each region of interest
% keep only the first component

%-------------------------------------%
%-loop over regions of interest
trial = [];
for i1 = 1:numel(roi)
  
  %-----------------%
  %-find channels belonging to region of interest
  newname = sprintf('%s', roi(i1).name);
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

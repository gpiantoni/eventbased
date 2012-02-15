function  mont = topodipole(dip, gerp)
%TOPODIPOLE create montage based on putative dipoles 
% based on grand average, no subject specific
% Use as
%   mont = topodipole(dip, gerp)
% where dip has fields:
%  .name = 'name of dipole location'
%  .time = time of the dipole, in seconds (can be one number or [.1 .2])
%  .cond = '@(avg) avg{1}.avg' (it means: take avg from the first gerp)
% and gerp is the average ERP for the condition of interest, cfg.erpeffect
% Output mont can be used in ft_prepare_montage(data, mont)

%-----------------%
%-prepare TRA
nnew = numel(dip);
nold = numel(gerp.label);

tra = zeros(nnew, nold);

for i = 1:numel(dip)
  dipavg = ft_selectdata(gerp, 'toilim', dip(i).time + [-.5 .5] * dip(i).wndw, 'avgovertime', 'yes');
  tra(i, :) = dipavg.avg';
end
%-----------------%

%-----------------%
%-prepare MONT
mont.labelorg = gerp.label;
mont.labelnew = {dip.name}';
mont.tra = tra;
%-----------------%

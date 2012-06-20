function [leadchan] = prepare_leadchan(lead, datachan)
%PREPARE_LEADCHAN only use good channels in leadfield
%
% Part of EVENTBASED/PRIVATE

%-----------------%
%-optional: remove channels if the vol does not cover them (f.e. channels
% on the cheek or channels on the neck)
leadlabel = lead.cfg.elec.label;
%TODO: make the code below more robust
highchan = find(lead.cfg.elec.chanpos(:,3) - min(lead.cfg.vol.bnd(1).pnt(:,3)) >= 0); % index of the channels above the lowest point

leadinside = lead.inside;
if size(leadinside,1) ~= 1; leadinside = leadinside'; end
for l = leadinside
  leadchan.leadfield{l} = lead.leadfield{l}(highchan,:);
end
leadlabel = leadlabel(highchan);
%-----------------%

%-----------------%
%-select common channels
[~, ichan] = intersect(leadlabel, datachan);
ichan = sort(ichan); % ichan in a good order

for l = leadinside
  leadchan.leadfield{l} = leadchan.leadfield{l}(ichan,:);
end
%-----------------%
function [leadchan] = prepare_leadchan(lead, datachan)
%PREPARE_LEADCHAN only use good channels in leadfield
%
% Part of EVENTBASED/PRIVATE

[~, ichan] = intersect(lead.cfg.elec.label, datachan);
ichan = sort(ichan); % ichan in a good order

leadchan = lead;
leadinside = lead.inside;
if size(leadinside,1) ~= 1; leadinside = leadinside'; end
for l = leadinside
  leadchan.leadfield{l} = lead.leadfield{l}(ichan,:);
end
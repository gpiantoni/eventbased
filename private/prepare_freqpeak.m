function [param output] = prepare_freqpeak(cfg, powpeak, begtrl)
%PREPARE_LEADCHAN only use good channels in leadfield
% param has field
%   - name
%   - time
%   - wndw
%   - freq
%   - band
%   - dpss
%
% Part of EVENTBASED/PRIVATE

output = '';

%---------------------------%
%-NAME
param.name = powpeak.name;
%---------------------------%

%---------------------------%
%-TIME: time center for analysis window
param.time = powpeak.time;
%---------------------------%

%---------------------------%
%-WNDW: time window
param.wndw = powpeak.wndw;

% if there is baseline at all
if ~isempty(cfg.powsource.bline)
  
  %-----------------%
  %-check that baseline contains always data, otherwise shrink the time window
  begbline = cfg.powsource.bline - powpeak.wndw/2; % beginning of baseline
  if begbline < begtrl
    output = sprintf('%sPowpeak %s: window length was too long (% 3.2fs)\n', ...
      output, powpeak.name, powpeak.wndw);
    param.wndw = (begtrl - cfg.powsource.bline) * -2;
  end
  %-----------------%
  
end
%---------------------------%

%---------------------------%
%-FREQ: main frequency for analysis
param.freq = powpeak.freq;
%---------------------------%

%---------------------------%
%-BAND: band of interest
% check that time window is long enough for dpss smoothing
% you get 0.75 from dpss(a,b). If b < 0.75, you get one taper (but the
% last one is not used by fieldtrip). If it's between 0.75 and 1.25,
% you get two tapers, but the second is not used and you get a warning
% from fieldtrip. Still, the one taper in use is not a simple hanning
param.band = powpeak.band/2;

if param.wndw * param.band > 0.75
  
  param.dpss = true;
  output = sprintf('%sPowpeak %s: window length % 3.2fs, using dpss (tapsmofrq: +-% 2.1f Hz)\n', ...
    output, powpeak.name, param.wndw, param.band);
  
else
  param.dpss = false; % hanning, no smoothing
  output = sprintf('%sPowpeak %s: window length % 3.2fs, using hanning\n', ...
    output, powpeak.name, param.wndw);
  
end
%---------------------------%

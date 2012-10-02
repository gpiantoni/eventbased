function peaks = get_peak(info, peak, type)
%GET_PEAK get the peaks for the source analysis
% Use as:
%  peaks = get_peak(info, peak, type)
%
%-with ERP
% INFO
%  .derp: directory with ERP data
% 
% PEAK: description of the peaks
%  - manually, as struct
%    .name: 'name_of_the_peak'
%    .time: center of the time window in s
%    .wndw: length of the time window in s
%  - automatically, as comparison of interest (one cell of cfg.opt.comp)
%
%-with POW or POWCORR
% INFO
%  .dpow: directory with POW data
% 
% PEAK: description of the peaks
%  - manually, as struct
%    .name: 'name_of_the_peak'
%    .time: center of the time window in s
%    .wndw: length of the time window in s
%    .freq: center of the frequency in Hz
%    .band: width of the frequency band in Hz
%  - automatically, as comparison of interest (one cell of cfg.opt.comp)
%
%  TYPE
%  'erp' for erp_peaks, 'pow' for pow_peaks, 'powcorr' for powcorr_peaks
%
%  PEAKS
%  erp_peak with time and window info
%  pow_peak, powcorr_peak with time and window info, freq and band info.
%
% Part of EVENTBASED/PRIVATE

if isstruct(peak)
  peaks = peak;
  
else
  
  %---------------------------%
  %-read from file
  peakname = regexprep(peak{1}, '*', '');
  if numel(peak) == 2;
    peakname = [peakname '_' regexprep(peak{2}, '*', '')];
  end
  
  switch type
    case 'erp'
      load([info.derp 'erp_peak_' peakname], 'erp_peak')
      peaks = erp_peak;
      
    case 'pow'
      load([info.dpow 'pow_peak_' peakname], 'pow_peak')
      peaks = pow_peak;
      
    case 'powcorr'
      load([info.dpow 'powcorr_peak_' peakname], 'powcorr_peak')
      peaks = powcorr_peak;
      
  end
  %---------------------------%
  
end
function peaks = getpeak(cfg, type)
%GETPEAK simply get the peaks for the analysis, similar for ERP and POW
% Use as:
%  peaks = getpeak(cfg, type)
%
%  CFG
%  .powsource.peaks: how to speficy peaks to analyze, 'manual' or 'pow_peak'
%          (peaks from grandpow) or 'powcorrpeak' (peaks from grandpowcorr)
%    if 'manual'
%      .powsource.pow_peak(1).name: string ('name_of_the_time_window')
%      .powsource.pow_peak(1).time: scalar (center of the time window in s)
%      .powsource.pow_peak(1).wndw: scalar (length of the time window in s)
%      .powsource.pow_peak(1).freq = 10; % center of the frequency
%      .powsource.pow_peak(1).band = 4; % width of the frequency band
%    if 'pow_peak'
%      .pow.refcomp: string of the comparison whose peaks will be localized
%    if 'powcorrpeak'
%      .powcorr.refcomp: string of the comparison whose peaks will be localized
%
%  TYPE
%  'erp' for erp_peaks or 'pow' for pow_peaks
%
%  PEAKS
%  erp_peak with time and window info
%  pow_peak with time and window info, freq and band info.
%
% Part of EVENTBASED/PRIVATE

switch type
  case 'erp'
    %---------------------------%
    %-use predefined or ERP-peaks for areas of interest
    switch cfg.erpsource.peaks
      case 'manual'
        erp_peak = cfg.erpsource.erp_peak;
        
      case 'erp_peak'
        peakname = regexprep(cfg.erpsource.refcomp{1}, '*', '');
        if numel(cfg.erpsource.refcomp) == 2;
          peakname = [peakname '_' regexprep(cfg.erpsource.refcomp{2}, '*', '')];
        end
        load([cfg.derp 'erp_peak_' peakname], 'erp_peak')
        
    end
    
    peaks = erp_peak;
    %---------------------------%
    
  case 'pow'
    %---------------------------%
    %-use predefined or power-peaks for areas of interest
    if strcmp(cfg.powsource.peaks, 'manual')
      pow_peak = cfg.powsource.pow_peak;
      
    elseif strcmp(cfg.powsource.areas, 'pow_peak')
      peakname = regexprep(cfg.pow.refcomp, '*', '');
      load([cfg.dpow 'pow_peak_' peakname], 'pow_peak')
      
    elseif strcmp(cfg.powsource.areas, 'powcorrpeak')
      peakname = regexprep(cfg.powcorr.refcomp, '*', '');
      load([cfg.dpow 'powcorrpeak_' peakname], 'powcorrpeak')
      pow_peak = powcorrpeak;
      
    end
    
    peaks = pow_peak;
    %---------------------------%
    
end
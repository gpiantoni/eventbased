function peaks = get_peak(cfg, type)
%GETPEAK simply get the peaks for the analysis, similar for ERP and POW
% Use as:
%  peaks = get_peak(cfg, type)
%
%  CFG with type == 'erp'
%  .erpsource.areas: how to speficy peaks to analyze, 'manual' or 'erp_peak'
%          (peaks from granderp)
%    if 'manual'
%      .erpsource.erp_peak(1).name: string ('name_of_the_time_window')
%      .erpsource.erp_peak(1).time: scalar (center of the time window in s)
%      .erpsource.erp_peak(1).wndw: scalar (length of the time window in s)
%    if 'erp_peak'
%      .erp.refcomp: string of the comparison whose peaks will be localized
%
%  CFG with type == 'pow'
%  .powsource.areas: how to speficy peaks to analyze, 'manual' or 'pow_peak'
%          (peaks from grandpow) or 'powcorr_peak' (peaks from grandpowcorr)
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
    switch cfg.erpsource.areas
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
    if strcmp(cfg.powsource.areas, 'manual')
      pow_peak = cfg.powsource.pow_peak;
      
    elseif strcmp(cfg.powsource.areas, 'pow_peak') ...
        || strcmp(cfg.powsource.areas, 'powcorr_peak')
      
      peakname = regexprep(cfg.powsource.refcomp{1}, '*', '');
      if numel(cfg.powsource.refcomp) == 2;
        peakname = [peakname '_' regexprep(cfg.powsource.refcomp{2}, '*', '')];
      end
      
      if strcmp(cfg.powsource.areas, 'pow_peak')
        load([cfg.dpow 'pow_peak_' peakname], 'pow_peak')
      else
        load([cfg.dpow 'powcorr_peak_' peakname], 'powcorr_peak')
        pow_peak = powcorr_peak;
      end
    end
    
    peaks = pow_peak;
    %---------------------------%
    
end
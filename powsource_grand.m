function powsource_grand(info, opt)
%POWSOURCE_GRAND group-level analysis of POW source data
%
% INFO
%  .log: name of the file and directory to save log
%  .dpow: directory with POW data
%  .mri.template: template of the averaged MRI
%
% CFG.OPT.
%  .comp*: comparisons to test (cell within cell, e.g. {{'cond1'} {'cond2'}})
%          You can only have one cell (while in POW_GRAND you can have 2),
%          because statistics with different filters is not correct. Use POWSTAT_SUBJ to reuse filters
%  .peak*: peak of interest (see GET_PEAK)
%
% IN
%  [info.dpow 'powsource_SUBJ_COND'] 'powsource_subj_A': source data for period of interest for each subject
%  [info.dpow 'powsource_SUBJ_COND'] 'powsource_subj_B': source data for baseline for each subject
%
% OUT
%  [info.dpow 'powsource_COND'] 'powsource': source analysis for all subject
%  [info.dpow 'powsource_peak_COND'] 'powsource_peak': significant source peaks in POW
%
% FIGURES
%  gpow_peak_COND_POWPEAK: 3d plot of the source for one peak
%
% * indicates obligatory parameter
%
% Part of EVENTBASED group-analysis
% see also ERP_SUBJ, ERP_GRAND,
% ERPSOURCE_SUBJ, ERPSOURCE_GRAND, ERPSTAT_SUBJ, ERPSTAT_GRAND,
% POW_SUBJ, POW_GRAND, POW_GRP, POWCORR_SUBJ, POWCORR_GRAND,
% POWSOURCE_SUBJ, POWSOURCE_GRAND, POWSTAT_SUBJ, POWSTAT_GRAND,
% SOURCE_SUBJ, CONN_SUBJ, CONN_GRAND, CONN_STAT

%---------------------------%
%-start log
output = sprintf('%s began at %s on %s\n', ...
  mfilename,  datestr(now, 'HH:MM:SS'), datestr(now, 'dd-mmm-yy'));
tic_t = tic;
%---------------------------%

pow_peak = get_peak(info, opt.peak, 'pow');

%---------------------------%
%-load template
template = ft_read_mri(info.mri.template);
%---------------------------%

%---------------------------------------------------------%
%-statistics for main effects
%-----------------------------------------------%
%-loop over conditions
for t = 1:numel(opt.comp)
  
  %-------------------------------------%
  %-loop over peaks
  powsource_peak = [];
  powsource = [];
  for p = 1:numel(pow_peak)
    output = sprintf('%s\n%s:\n', output, pow_peak(p).name);
    
    %---------------------------%
    %-statistics for effects of interest
    if numel(opt.comp{t}) == 1
      
      %-----------------%
      %-compare against zero
      cond = opt.comp{t}{1};
      comp = regexprep(cond, '*', '');
      output = sprintf('%s\n   COMPARISON %s\n', output, cond);
      %-----------------%
      
      %-----------------%
      %-file for each cond
      [outtmp data] = load_subj(info, 'powsource', cond);
      output = [output outtmp];
      if isempty(data); continue; end
      %-----------------%
      
      %-----------------%
      %-pow or coh
      if isfield(data{1}.avg, 'coh') % coh wins
        opt.parameter = 'coh';
      else
        opt.parameter = 'pow';
      end
      %-----------------%
      
      %-----------------%
      %-grand average
      cfg = [];
      cfg.keepindividual = 'yes';
      cfg.parameter = opt.parameter;
      gpowsouPre = ft_sourcegrandaverage(cfg, data{:,1,p});
      gpowsource = ft_sourcegrandaverage(cfg, data{:,2,p});
      %-----------------%
      
      %-----------------%
      %-do stats
      [soupos powsource{p,1} outtmp] = report_source(opt, gpowsource, gpowsouPre);
      powsource_peak(p).pos = soupos;
      powsource_peak(p).center = mean(soupos,1);
      powsource_peak(p).name = pow_peak(p).name;
      output = [output outtmp];
      %-----------------%
      
    else
      output = [output sprintf('It does not make sense to compare conditions using different filters\nPlease use powstat for comparisons\n')];
      continue
      
    end
    
    %-----------------%
    %-plot source
    h = figure('vis', 'off');
    plot_volume(powsource(p, :), template, 'stat')
    
    %--------%
    pngname = sprintf('gpow_peak_%s_%s', comp, pow_peak(p).name);
    saveas(h, [info.log filesep pngname '.png'])
    close(h); drawnow
    
    [~, logfile] = fileparts(info.log);
    system(['ln ' info.log filesep pngname '.png ' info.rslt pngname '_' logfile '.png']);
    %--------%
    %-----------------%
    %---------------------------%
    
  end
  %-------------------------------------%
  
  %---------------------------%
  %-save
  save([info.dpow 'powsource_peak_' comp], 'powsource_peak')
  
  for p = 1:numel(powsource)
    powsource{p}.cfg = []; % this is huge
  end
  save([info.dpow 'powsource_' comp], 'powsource', '-v7.3')
  %---------------------------%
  
end
%---------------------------%
%---------------------------------------------------------%

%---------------------------%
%-end log
toc_t = toc(tic_t);
outtmp = sprintf('%s ended at %s on %s after %s\n\n', ...
  mfilename, datestr(now, 'HH:MM:SS'), datestr(now, 'dd-mmm-yy'), ...
  datestr( datenum(0, 0, 0, 0, 0, toc_t), 'HH:MM:SS'));
output = [output outtmp];

%-----------------%
fprintf(output)
fid = fopen([info.log '.txt'], 'a');
fwrite(fid, output);
fclose(fid);
%-----------------%
%---------------------------%
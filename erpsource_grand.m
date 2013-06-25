function erpsource_grand(info, opt)
%ERPSOURCE_GRAND group-level analysis of ERP source data
%
% INFO
%  .log: name of the file and directory to save log
%  .derp: directory with ERP data
%  .mri.template: template of the averaged MRI
%
% CFG.OPT.
%  .comp*: comparisons to test (cell within cell, e.g. {{'cond1'} {'cond2'}})
%          You can only have one cell (while in ERP_GRAND you can have 2),
%          because statistics with different filters is not correct. Use ERPSTAT_SUBJ to reuse filters
%  .peak*: peak of interest (see GET_PEAK)
%
% IN
%  [info.derp 'erpsource_SUBJ_COND'] 'erpsource_subj_A': source data for period of interest for each subject
%  [info.derp 'erpsource_SUBJ_COND'] 'erpsource_subj_B': source data for baseline for each subject
%
% OUT
%  [info.derp 'erpsource_COND'] 'erpsource': source analysis for all subject
%
% FIGURES
%  gerp_peak_COND_ERPPEAK: 3d plot of the source for one peak
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

erp_peak = get_peak(info, opt.peak, 'erp');

%---------------------------%
%-load template TODO: find elegant solution for this
load(info.mri.template)
template = mri;
% template = ft_read_mri(info.mri.template);
%---------------------------%

%---------------------------------------------------------%
%-statistics for main effects
%-----------------------------------------------%
%-loop over conditions
for t = 1:numel(opt.comp)
  
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
    [outtmp data] = load_subj(info, 'erpsource', cond);
    output = [output outtmp];
    if isempty(data); continue; end
    %-----------------%
    
    %-----------------%
    %-erp or coh
    if isfield(data{1}.avg, 'coh') % coh wins
      opt.parameter = 'coh';
    else
      opt.parameter = 'pow';
    end
    %-----------------%
    
  else
    output = [output sprintf('It does not make sense to compare conditions using different filters\nPlease use erpstat for comparisons\n')];
    continue
    
  end
  %---------------------------%
  
  %-------------------------------------%
  %-loop over peaks
  erpsource = [];
  for p = 1:numel(erp_peak)
    output = sprintf('%s\n%s:\n', output, erp_peak(p).name);
    
    %---------------------------%
    %-----------------%
    %-grand average
    cfg = [];
    cfg.keepindividual = 'yes';
    cfg.parameter = opt.parameter;
    gerpsouPre = ft_sourcegrandaverage(cfg, data{:,1,p});
    gerpsource = ft_sourcegrandaverage(cfg, data{:,2,p});
    %-----------------%
    
    %-----------------%
    %-do stats
    [erpsource{p,1} outtmp] = report_source(opt, gerpsource, gerpsouPre);
    erpsource{p}.name = erp_peak(p).name;
    output = [output outtmp];
    %-----------------%
    
    %-----------------%
    %-plot source
    h = figure('vis', 'off');
    plot_volume(erpsource(p, :), template, 'stat')
    
    %--------%
    pngname = sprintf('gerp_peak_%s_%s', comp, erp_peak(p).name);
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
  for p = 1:numel(erpsource)
    erpsource{p}.cfg = []; % this is huge
  end
  save([info.derp 'erpsource_' comp], 'erpsource', '-v7.3')
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

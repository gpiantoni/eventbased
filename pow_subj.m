function pow_subj(cfg, subj)
%POW_SUBJ create subject-specific pow
%
% CFG
%  .data: name of projects/PROJNAME/subjects/
%  .mod: name of the modality used in recordings and projects
%  .cond: name to be used in projects/PROJNAME/subjects/0001/MOD/CONDNAME/
%  .endname: includes previous steps '_seldata_gclean_preproc_redef'
%  .log: name of the file and directory with analysis log
%  .test: a cell with the condition defined by redef. This function will loop over cfg.test
%  .dpow: directory to save POW data
%
%  .pow: a structure with cfg to pass to ft_freqanalysis
%  .pow.outliers: logical (print tables with number of points above a
%  certain number of standard deviation, experimental code)
%  .pow.outliersthr: if outliers is true, this is the threshold to reject a
%  trial (if empty, no rejection. It's very hard to give an indicative
%  value for this number. First try checking the plots, then assign a threshold)
%
% Baseline correction is applied in POW_GRAND
%
% OUT
%  [cfg.dpow 'pow_001_TEST']: power analysis for single-subject
%
% Part of EVENTBASED single-subject
% see also ERP_SUBJ, ERP_GRAND, ERPSOURCE_SUBJ, ERPSOURCE_GRAND,
% POW_SUBJ, POW_GRAND, POWSOURCE_SUBJ, POWSOURCE_GRAND,
% POWCORR_SUBJ, POWCORR_GRAND,
% CONN_SUBJ, CONN_GRAND, CONN_STAT

%---------------------------%
%-start log
output = sprintf('(p%02.f) %s started at %s on %s\n', ...
  subj, mfilename,  datestr(now, 'HH:MM:SS'), datestr(now, 'dd-mmm-yy'));
tic_t = tic;
%---------------------------%

%---------------------------%
%-dir and files
ddir = sprintf('%s%04.f/%s/%s/', cfg.data, subj, cfg.mod, cfg.cond); % data
%---------------------------%

%-------------------------------------%
%-loop over conditions
for k = 1:numel(cfg.test)
  
  %-----------------%
  %-input and output for each condition
  allfile = dir([ddir cfg.test{k} cfg.endname '.mat']); % files matching a preprocessing
  
  condname = regexprep(cfg.test{k}, '*', '');
  outputfile = sprintf('pow_%02.f_%s', subj, condname);
  %-----------------%
  
  %-----------------%
  %-concatenate only if you have more datasets
  if numel(allfile) > 1
    spcell = @(name) sprintf('%s%s', ddir, name);
    allname = cellfun(spcell, {allfile.name}, 'uni', 0);
    
    cfg1 = [];
    cfg1.inputfile = allname;
    data = ft_appenddata(cfg1);
    
  elseif numel(allfile) == 1
    load([ddir allfile(1).name], 'data')
    
  else
    output = sprintf('%sCould not find any file in %s for test %s\n', ...
      output, ddir, cfg.test{k});
    continue
    
  end
  %-----------------%
  
  %-----------------%
  %-calculate power
  cfg2 = cfg.pow;
  cfg2.feedback = 'none';
  
  if cfg.pow.outliers
    cfg2.keeptrials = 'yes';
  end
  
  freq = ft_freqanalysis(cfg2, data);
  
  if isfield(cfg.pow, 'toi')
    freq.time = cfg.pow.toi;
  end
  %-----------------%
  
  %-----------------%
  %-deal with outliers
  if cfg.pow.outliers
    
    %-------%
    %-compute sd across trials
    trlsd = zeros(size(freq.powspctrm,1), 1);
    for i = 1:size(freq.powspctrm,1)
      i_powspctrm = freq.powspctrm(i, :, :, :);
      trlsd(i) = nanstd(i_powspctrm(:),1);
    end
    %-------%
    
    %-------%
    %-compute average (and reject outliers if necessary)
    cfg4 = [];
    
    if isfield(cfg.pow, 'outliersthr') && ~isempty(cfg.pow.outliersthr)
      goodtrl = find(trlsd < cfg.pow.outliersthr);
      badtrl = find(trlsd >= cfg.pow.outliersthr);
      
      output = [output sprintf('Cond %s Rejecting %5.f outliers out of %5.f trials, remaining %5.f trials\n', ...
        condname, numel(badtrl), numel(trlsd), numel(goodtrl))];
      
      
      cfg4.trials = goodtrl;
      
    end
    
    freq = ft_freqdescriptives(cfg4, freq);
    %-------%
    
    %-------%
    %-plot feedback on montage
    figure
    plot(trlsd, '.')
    
    if isfield(cfg.pow, 'outliersthr') && ~isempty(cfg.pow.outliersthr)
      hold on
      plot(badtrl, trlsd(badtrl), '+r')
    end
    
    legend('trials', 'bad trials')
    title('outliers')
    xlabel('n trials')
    ylabel('s.d.')
    
    pngfile = [cfg.log filesep 'powoutliers_' sprintf('%03.f', subj) '_' condname '.png'];
    saveas(gcf, pngfile);
    close(gcf); drawnow
    %-------%
    
  end
  %-----------------%
  
  save([cfg.dpow outputfile], 'freq')
  
end
%-------------------------------------%

%---------------------------%
%-end log
toc_t = toc(tic_t);
outtmp = sprintf('(p%02.f) %s ended at %s on %s after %s\n\n', ...
  subj, mfilename, datestr(now, 'HH:MM:SS'), datestr(now, 'dd-mmm-yy'), ...
  datestr( datenum(0, 0, 0, 0, 0, toc_t), 'HH:MM:SS'));
output = [output outtmp];

%-----------------%
fprintf(output)
fid = fopen([cfg.log '.txt'], 'a');
fwrite(fid, output);
fclose(fid);
%-----------------%
%---------------------------%

function powcorr_subj(cfg, subj)
%POWCORR_SUBJ correlate power with trialinfo
% practically identical to pow_subj but it takes log and it correlates with
% cfg.powcorr.info 

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
  allfile = dir([ddir '*' cfg.test{k} cfg.endname '.mat']); % files matching a preprocessing
  if isempty(allfile)
    continue
  end
  
  condname = regexprep(cfg.test{k}, '*', '');
  outputfile = sprintf('powcorr_%02.f_%s', subj, condname);
  %-----------------%
  
  %-----------------%
  %-concatenate only if you have more datasets
  if numel(allfile) > 1
    spcell = @(name) sprintf('%s%s', ddir, name);
    allname = cellfun(spcell, {allfile.name}, 'uni', 0);
    
    cfg1 = [];
    cfg1.inputfile = allname;
    data = ft_appenddata(cfg1);
    
  else
    load([ddir allfile(1).name], 'data')
    
  end
  %-----------------%
  
  %-----------------%
  %-calculate power
  cfg2 = cfg.powcorr;
  cfg2.feedback = 'none';
  cfg2.keeptrials = 'yes';
  data = ft_freqanalysis(cfg2, data);
  data.time = cfg2.toi;
  %-----------------%
  
  %-----------------%
  %-fix baseline
  if ~isempty(cfg.powcorr.bl.baseline)
    cfg3 = cfg.powcorr.bl;
    freq = ft_freqbaseline(cfg3, data);
    
  else
    freq = data;
    
  end
  %-----------------%
  
  %-----------------%
  %-regression at each point
  freq.powspctrm = log(freq.powspctrm);
  
  [s1 s2 s3 s4] = size(freq.powspctrm);
  
  powspctrm = nan(s2,s3,s4);
  warning('off', 'MATLAB:rankDeficientMatrix') % NaN in \
  
  for i2 = 1:s2
    for i3 = 1:s3
      for i4 = 1:s4
        regr = [ones(size(freq.trialinfo,1),1) freq.powspctrm(:,i2,i3,i4)];
        beta = regr \ freq.trialinfo(:, cfg.powcorr.info);
        powspctrm(i2,i3,i4) = beta(2);
      end
    end
  end
  warning('on', 'MATLAB:rankDeficientMatrix')
  %-----------------%
  
  %-----------------%
  %-restructure freq for multiple subjects
  freq = rmfield(freq, 'trialinfo');
  freq.dimord = freq.dimord(5:end);
  freq.powspctrm = powspctrm;
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

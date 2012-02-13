function powcorr_grand(cfg)
%POWCORR_GRAND grand power average for powcorr
% practically identical to pow_grand but not powpeaks
% pow -> powcorr; cfg.gpow -> cfg.gpowcorr

mversion = 2;
%02 12/02/07 cfg.poweffect -> cfg.powcorreffect
%01 12/02/06 practically identical to pow_grand

%---------------------------%
%-start log
output = sprintf('%s (v%02.f) started at %s on %s\n', ...
  mfilename,  mversion, datestr(now, 'HH:MM:SS'), datestr(now, 'dd-mmm-yy'));
tic_t = tic;
%---------------------------%

%---------------------------%
%-loop over conditions
gpow = [];
for k = 1:numel(cfg.test)
  
  %-----------------%
  %-file for each cond
  condname = regexprep(cfg.test{k}, '*', '');
  inputfile = sprintf('powcorr_*_%s.mat', condname);
  
  allsub = dir([cfg.dpow inputfile]);
  
  if isempty(allsub)
    outtmp = sprintf('%s does not match any file\n', condname);
    output = [output outtmp];
    continue
  end
  %-----------------%
  
  %-----------------%
  %-POW over subj
  spcell = @(name) sprintf('%s%s', cfg.dpow, name);
  allname = cellfun(spcell, {allsub.name}, 'uni', 0);
  
  cfg1 = [];
  cfg1.inputfile = allname;
  cfg1.keepindividual = 'yes';
  gfreq{k} = ft_freqgrandaverage(cfg1);
  
  cfg2 = [];
  cfg2.variance = 'yes';
  gpow{k} = ft_freqdescriptives(cfg2, gfreq{k});
  gpow{k}.tscore =  gpow{k}.powspctrm ./ gpow{k}.powspctrmsem;
  %-----------------%
  
end

%-----------------%
%-save
save([cfg.dpow cfg.proj '_grandpowcorr'], 'gpow')
%-----------------%
%---------------------------%

if ~isempty(gpow)
    
    %---------------------------%
  %-statistics for main effects
  [powcorrpeak outtmp] = reportcluster(gfreq{cfg.powcorreffect}, cfg);
  
  save([cfg.dpow cfg.proj '_powcorrpeak'], 'powcorrpeak')
  output = [output outtmp];
  %---------------------------%
  
  %---------------------------%
  %-feedback
  load(cfg.sens.layout, 'layout');
  
%   %-----------------%
%   %-singleplotTFR
%   for t = 1:numel(cfg.test)
%     
%     %--------%
%     %-loop over freq
%     for c = 1:numel(cfg.gpowcorr.chan)
%       
%       figure
%       cfg2 = [];
%       cfg2.channel = cfg.gpowcorr.chan(c).chan;
%       cfg2.zlim = [-4 4];
%       cfg2.parameter = 'tscore';
%       cfg2.layout = layout;
%       ft_singleplotTFR(cfg2, gpow{t});
%       
%       title([cfg.test{t} ' ' cfg.gpowcorr.chan(c).name])
%       
%       pngname = [cfg.log filesep cfg.proj '_grandpowcorr_TFR' num2str(k) '_' num2str(c)];
%       saveas(gcf, pngname, 'png')
%       close(gcf); drawnow
%     end
%     %--------%
%     
%   end
%   %-----------------%
  %---------------------------%
  
end

%---------------------------%
%-end log
toc_t = toc(tic_t);
outtmp = sprintf('%s (v%02.f) ended at %s on %s after %s\n\n', ...
  mfilename, mversion, datestr(now, 'HH:MM:SS'), datestr(now, 'dd-mmm-yy'), ...
  datestr( datenum(0, 0, 0, 0, 0, toc_t), 'HH:MM:SS'));
output = [output outtmp];

%-----------------%
fprintf(output)
fid = fopen([cfg.log '.txt'], 'a');
fwrite(fid, output);
fclose(fid);
%-----------------%
%---------------------------%
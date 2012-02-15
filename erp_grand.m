function erp_grand(cfg)
%ERP_GRAND grand time lock analysis

mversion = 9;
%09 12/02/08 nicer plots and link to results in cfg.rslt
%08 12/02/06 deal with cases when gerp is empty
%07 12/02/05 make topoplot as well
%06 12/02/03 renamed to erp_grand
%05 12/01/12 include reportcluster to get erppeak
%04 11/11/21 also compute keepindividual = 'yes'
%03 11/11/20 fix wild card in DIR
%02 11/09/27 cfg.erp.cond -> cfg.test
%01 11/07/21 created

%---------------------------%
%-start log
output = sprintf('%s (v%02.f) started at %s on %s\n', ...
  mfilename,  mversion, datestr(now, 'HH:MM:SS'), datestr(now, 'dd-mmm-yy'));
tic_t = tic;
%---------------------------%

%---------------------------%
%-loop over conditions
gerp = [];
for k = 1:numel(cfg.test)
  
  %-----------------%
  %-file for each cond
  condname = regexprep(cfg.test{k}, '*', '');
  inputfile = sprintf('erp_*_%s.mat', condname);
  
  allsub = dir([cfg.derp inputfile]);
  
  if isempty(allsub)
    outtmp = sprintf('%s does not match any file\n', condname);
    output = [output outtmp];
    continue
  end
  %-----------------%
  
  %-----------------%
  %-erp over subj
  spcell = @(name) sprintf('%s%s', cfg.derp, name);
  allname = cellfun(spcell, {allsub.name}, 'uni', 0);
  
  cfg1 = [];
  cfg1.inputfile = allname;
  gerp{k} = ft_timelockgrandaverage(cfg1);
  cfg1.keepindividual = 'yes';
  gerpall{k} = ft_timelockgrandaverage(cfg1);
  %-----------------%
  
end

%-----------------%
%-save
save([cfg.derp cfg.proj '_granderp'], 'gerp')
%-----------------%
%---------------------------%

%-----------------------------------------------%
%-stats and plots
load(cfg.sens.layout, 'layout');

if ~isempty(gerp)
  
  %---------------------------%
  %-statistics for main effects
  [erppeak outtmp] = reportcluster(gerpall{cfg.erpeffect}, cfg);
  
  save([cfg.derp cfg.proj '_erppeak'], 'erppeak')
  output = [output outtmp];
  %---------------------------%
  
  %---------------------------%
  %-singleplotER (multiple conditions at once)
  for c = 1:numel(cfg.gerp.chan)
    
    %--------%
    %-figure
    h = figure;
    set(h, 'Renderer', 'painters')
    %--------%
    
    %--------%
    %-plot
    cfg3 = [];
    cfg3.channel = cfg.gerp.chan(c).chan;
    cfg3.baseline = cfg.gerp.bline;
    cfg3.ylim = 'maxabs';
    ft_singleplotER(cfg3, gerp{:});
    
    legend(cfg.test)
    
    title(cfg.gerp.chan(c).name)
    %--------%
    
    %--------%
    %-save and link
    pngname = sprintf('gerp_erp_c%02.f', c);
    saveas(gcf, [cfg.log filesep pngname '.png'])
    close(gcf); drawnow
    
    [~, logfile] = fileparts(cfg.log);
    system(['ln ' cfg.log filesep pngname '.png ' cfg.rslt pngname '_' logfile '.png']);
    %--------%
    
  end
  %---------------------------%
  
  %---------------------------%
  %-topoplotTFR (loop over tests)
  for t = 1:numel(cfg.test)
    
    %--------%
    %-figure
    h = figure;
    set(h, 'Renderer', 'painters')
    %--------%
    
    %--------%
    %-plot
    colorlim = max(max(abs(gerp{cfg.erpeffect}.avg)));
    timelim = gerp{cfg.erpeffect}.time([1 end]);
    
    cfg4 = [];
    cfg4.xlim = timelim(1):.1:timelim(2); % one plot every 100 ms
    cfg4.zlim = [-1 1] * colorlim;
    cfg4.layout = layout;
    cfg4.style = 'straight';
    cfg4.marker = 'off';
    cfg4.comment = 'xlim';
    cfg4.commentpos = 'title';

    ft_topoplotER(cfg4, gerp{t});
    %--------%
    
    %--------%
    %-save and link
    condname = regexprep(cfg.test{t}, '*', '');
    pngname = sprintf('gerp_topo_%s', condname);
    saveas(gcf, [cfg.log filesep pngname '.png'])
    close(gcf); drawnow
    
    [~, logfile] = fileparts(cfg.log);
    system(['ln ' cfg.log filesep pngname '.png ' cfg.rslt pngname '_' logfile '.png']);
    %--------%
    
  end
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
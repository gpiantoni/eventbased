function erpsource_subj(cfg, subj)
%LCMV on interesting parts, defined by cfg.erpsource.source

mversion = 7;
%07 12/02/02 renamed to  erpsource_subj
%06 12/01/16 calculate only erpeffect (otherwise the file is too big)
%05 12/01/13 calculates nai based on baseline
%04 12/01/13 lcmv power depends on length of covariance window, now baseline and source have the same baseline duration
%03 12/01/12 use template or subject-specific vol
%02 12/01/12 save souPre as well, but don't calculate NAI
%01 12/01/11 created from gosdpowsource

%-----------------%
%-input
if nargin == 1
  subj = cfg.subj;
end
%-----------------%

%---------------------------%
%-start log
output = sprintf('(p%02.f) %s (v%02.f) started at %s on %s\n', ...
  subj, mfilename,  mversion, datestr(now, 'HH:MM:SS'), datestr(now, 'dd-mmm-yy'));
tic_t = tic;
%---------------------------%

%---------------------------%
%-dir and files
ddir = sprintf('%s%04.f/%s/%s/', cfg.data, subj, cfg.mod, cfg.cond); % data

%-----------------%
%-head shape
if strcmp(cfg.voltype, 'template')
  load(cfg.leadfile, 'vol', 'lead', 'sens')
  
else
  mod = 'smri';
  cond = 't1';
  mdir = sprintf('%s%04.f/%s/%s/', cfg.data, subj, mod, cond); % mridata dir
  mfile = sprintf('%s_%s_%04.f_%s_%s', cfg.proj, cfg.rec, subj, mod, cond); % mridata
  
  load([mdir mfile '_elec.mat'], 'elec')
  sens = elec;
  load([mdir mfile '_vol_' cfg.voltype '.mat'], 'vol')
  load([mdir mfile '_lead_' cfg.voltype '.mat'], 'lead')
  
end
%-----------------%

%-----------------%
%-use predefined or ERP-peaks for areas of interest
if strcmp(cfg.erpsource.areas, 'manual')
  erppeak = cfg.erpsource.erppeak;
  
elseif strcmp(cfg.erpsource.areas, 'erppeak')
  load([cfg.derp cfg.proj '_erppeak'], 'erppeak')
  
end
%-----------------%
%---------------------------%

%-------------------------------------%
%-loop over conditions
for e = 1:numel(cfg.erpeffect)
  k = cfg.erpeffect(e);
  
  %-----------------%
  %-input and output for each condition
  allfile = dir([ddir '*' cfg.test{k} cfg.endname '.mat']); % files matching a preprocessing
  if isempty(allfile)
    continue
  end
  
  condname = regexprep(cfg.test{k}, '*', '');
  outputfile = sprintf('erpsource_%02.f_%s', subj, condname);
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
  %-find bad channels (bad channel if bad in at least one dataset)
  if ~iscell(data.cfg.previous) % not appenddata
    badchan = ft_findcfg(data.cfg, 'badchannel');
    
  else
    
    badchan = {};
    for i = 1:numel(data.cfg.previous)
      badtemp = ft_findcfg(data.cfg.previous{i}, 'badchannel');
      badchan = union(badchan, badtemp);
    end
    
  end
  %-----------------%
  
  %-----------------%
  %-remove bad channels from leadfield
  datachan = ft_channelselection([{'all'}; cellfun(@(x) ['-' x], badchan, 'uni', false)], data.label);
  
  [~, ichan] = intersect(lead.cfg.elec.label, datachan);
  ichan = sort(ichan); % ichan in a good order
  
  leadchan = lead;
  leadinside = lead.inside;
  if size(leadinside,1) ~= 1; leadinside = leadinside'; end
  for l = leadinside
    leadchan.leadfield{l} = lead.leadfield{l}(ichan,:);
  end
  %-----------------%
  
  for f = 1:numel(erppeak)
    
    %---------------------------%
    %-baseline
    %-----------------%
    %-covariance window
    cfg2 = cfg.erpsource.erp;
    cfg2.covariance = 'yes';
    cfg2.covariancewindow = cfg.erpsource.bline  + erppeak(f).wndw * [-.5 .5];
    cfg2.feedback = 'none';
    cfg2.channel = datachan;
    
    avgPre = ft_timelockanalysis(cfg2, data);
    %-----------------%
    
    %-----------------%
    %-source analysis
    cfg3              = [];
    
    cfg3.method   = 'lcmv';
    cfg3.lcmv.lambda = cfg.erpsource.lambda;
    cfg3.lcmv.powmethod = cfg.erpsource.powmethod;
    cfg3.lcmv.feedback = 'none';
    
    cfg3.vol      = vol;
    cfg3.grid     = leadchan;
    cfg3.elec     = sens;
    cfg3.feedback = 'none';
    
    souPre{f} = ft_sourceanalysis(cfg3, avgPre);
    %-----------------%
    %---------------------------%
    
    %---------------------------%
    %-main effect
    %-----------------%
    %-covariance
    cfg2.covariancewindow = erppeak(f).time + erppeak(f).wndw * [-.5 .5];
    avgPost = ft_timelockanalysis(cfg2, data);
    %-----------------%
    
    %-----------------%
    %-source
    cfg3.lcmv.keepfilter = 'yes'; % it's used by gosdconn
    cfg3.lcmv.realfilter = 'yes';
    cfg3.lcmv.keepmom    = 'no';
    source{f} = ft_sourceanalysis(cfg3, avgPost);
    
    source{f}.avg.nai = log(source{f}.avg.pow ./ souPre{f}.avg.pow);
    %-----------------%
    %---------------------------%
    
  end
  
  %-----------------%
  %-save source
  save([cfg.derp outputfile], 'source', 'souPre', '-v7.3')
  %-----------------%
  
end
%-------------------------------------%

%---------------------------%
%-end log
toc_t = toc(tic_t);
outtmp = sprintf('(p%02.f) %s (v%02.f) ended at %s on %s after %s\n\n', ...
  subj, mfilename, mversion, datestr(now, 'HH:MM:SS'), datestr(now, 'dd-mmm-yy'), ...
  datestr( datenum(0, 0, 0, 0, 0, toc_t), 'HH:MM:SS'));
output = [output outtmp];

%-----------------%
fprintf(output)
fid = fopen([cfg.log '.txt'], 'a');
fwrite(fid, output);
fclose(fid);
%-----------------%
%---------------------------%  
  


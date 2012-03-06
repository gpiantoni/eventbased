function erpsource_subj(cfg, subj)
%ERPSOURCE_SUBJ: identify sources from erp peaks using LCMV
%
% CFG
%  .data: name of projects/PROJNAME/subjects/
%  .mod: name of the modality used in recordings and projects
%  .cond: name to be used in projects/PROJNAME/subjects/0001/MOD/CONDNAME/
%  .endname: includes previous steps '_seldata_gclean_preproc_redef'
%  .test: a cell with the condition defined by redef. This function will loop over cfg.test
%  .derp: directory to save ERP data
%
%  .erp: a structure with cfg to pass to ft_timelockanalysis
%
%  .vol.type: 'template' or subject-specific ('dipoli' or 'openmeeg')
%    if template, specify template .vol.template (should contain vol, lead, sens)
%    if not template, specify 
%      .vol.mod: 'smri'
%      .vol.cond: 't1'
%      .proj: because the project name is part of the MRI name
%
%  .erpsource.areas: how to speficy peaks to analyze, 'manual' or 'erppeak' (peaks from granderp)
%    if 'manual'
%      .erpsource.erppeak(1).name = 'name_of_the_time_ window';
%      .erpsource.erppeak(1).time = 0.10; % center of the time window
%      .erpsource.erppeak(1).wndw = 0.05; % length of the time window
%    if 'erppeak', it reads the significant peaks calculated by erp_grand
%
%  .erpsource.erp: a structure with cfg to pass to ft_timelockanalysis (it's better if identical to cfg.erp)
%  .erpsource.bline: one number in s, the center of the covariance window of the baseline (the window length depends on erppeak)
%
%  .erpsource.lambda: regularization parameter of beamformer ('10%')
%  .erpsource.powmethod: power method of beamformer ('trace' or 'lambda1')
%
% OUT
%  [cfg.derp 'erpsource_001_TEST']: source data for period of interest and baseline for each subject
%
% Part of EVENTBASED single-subject
% see also ERP_SUBJ, ERP_GRAND, ERPSOURCE_SUBJ, ERPSOURCE_GRAND

%---------------------------%
%-start log
output = sprintf('(p%02.f) %s started at %s on %s\n', ...
  subj, mfilename,  datestr(now, 'HH:MM:SS'), datestr(now, 'dd-mmm-yy'));
tic_t = tic;
%---------------------------%

%---------------------------%
%-dir and files
ddir = sprintf('%s%04.f/%s/%s/', cfg.data, subj, cfg.mod, cfg.cond); % data

%-----------------%
%-head shape
if strcmp(cfg.vol.type, 'template')
  load(cfg.vol.template, 'vol', 'lead', 'sens')
  
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
  load([cfg.derp cfg.cond '_erppeak'], 'erppeak')
  
end
%-----------------%
%---------------------------%

%-------------------------------------%
%-loop over conditions
for e = 1:numel(cfg.erpeffect)
  k = cfg.erpeffect(e);
  
  %-----------------%
  %-input and output for each condition
  allfile = dir([ddir cfg.test{k} cfg.endname '.mat']); % files matching a preprocessing
  
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
    
  elseif numel(allfile) == 1
    load([ddir allfile(1).name], 'data')
    
  else
    output = sprintf('%sCould not find any file in %s for test %s\n', ...
      output, ddir, cfg.test{k});
    continue
    
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
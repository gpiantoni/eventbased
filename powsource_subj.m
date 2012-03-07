function powsource_subj(cfg, subj)
%ERPSOURCE_SUBJ: identify sources from pow peaks using DICS
%
% CFG
%  .data: name of projects/PROJNAME/subjects/
%  .mod: name of the modality used in recordings and projects
%  .cond: name to be used in projects/PROJNAME/subjects/0001/MOD/CONDNAME/
%  .endname: includes previous steps '_seldata_gclean_preproc_redef'
%  .test: a cell with the condition defined by redef. This function will loop over cfg.test
%  .dpow: directory to save ERP data
%
%  .pow: a structure with cfg to pass to ft_freqanalysis
%
%  .vol.type: 'template' or subject-specific ('dipoli' or 'openmeeg')
%    if template, specify template .vol.template (should contain vol, lead, sens)
%    if not template, specify 
%      .vol.mod: 'smri'
%      .vol.cond: 't1'
%      .proj: because the project name is part of the MRI name
%
%  .powsource.areas: how to speficy peaks to analyze, 'manual' or 'powpeak' (peaks from grandpow)
%    if 'manual'
%      .powsource.powpeak(1).name = 'name_of_the_time_ window';
%      .powsource.powpeak(1).time = 0.10; % center of the time window
%      .powsource.powpeak(1).wndw = 0.05; % length of the time window
%      .powsource.powpeak(1).freq = 10; % center of the frequency
%      .powsource.powpeak(1).band = 4; % width of the frequency band
%    if 'powpeak', it reads the significant peaks calculated by pow_grand
%
%  .powsource.pow: a structure with cfg to pass to ft_timelockanalysis (it's better if identical to cfg.pow)
%  .powsource.bline: one number in s, the center of the covariance window of the baseline (the window length depends on powpeak)
%
%  .powsource.lambda: regularization parameter of beamformer ('10%')
%  .powsource.powmethod: power method of beamformer ('trace' or 'lambda1')
%
% OUT
%  [cfg.dpow 'powsource_001_TEST']: source data for period of interest and baseline for each subject
%
% Part of EVENTBASED single-subject
% see also ERP_SUBJ, ERP_GRAND, ERPSOURCE_SUBJ, ERPSOURCE_GRAND, 
% POW_SUBJ, POW_GRAND, POWSOURCE_SUBJ, POWSOURCE_GRAND, 
% POWCORR_SUBJ, POWCORR_SUBJ,
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

%-----------------%
%-head shape
if strcmp(cfg.voltype, 'template')
  load(cfg.leadfile, 'vol', 'lead', 'sens')
  
else
  mdir = sprintf('%s%04.f/%s/%s/', cfg.data, subj, cfg.vol.mod, cfg.vol.cond); % mridata dir
  mfile = sprintf('%s_%04.f_%s_%s', cfg.rec, subj, cfg.vol.mod, cfg.vol.cond); % mridata
  
  load([mdir mfile '_elec.mat'], 'elec')
  sens = elec;
  load([mdir mfile '_vol_' cfg.vol.type '.mat'], 'vol')
  load([mdir mfile '_lead_' cfg.vol.type '.mat'], 'lead')
  
end
%-----------------%

%-----------------%
%-use predefined or power-peaks for areas of interest
if strcmp(cfg.powsource.areas, 'manual')
  powpeak = cfg.powsource.powpeak;
  
elseif strcmp(cfg.powsource.areas, 'powpeak')
  load([cfg.dpow cfg.proj '_powpeak'], 'powpeak')
  
elseif strcmp(cfg.powsource.areas, 'powcorrpeak')
  load([cfg.dpow cfg.proj '_powcorrpeak'], 'powcorrpeak')
  powpeak = powcorrpeak;
  
end
%-----------------%
%---------------------------%

%-------------------------------------%
%-loop over conditions
for e = 1:numel(cfg.poweffect)
  k = cfg.poweffect(e);
  
  %-----------------%
  %-input and output for each condition
  allfile = dir([ddir cfg.test{k} cfg.endname '.mat']); % files matching a preprocessing
  
  condname = regexprep(cfg.test{k}, '*', '');
  outputfile = sprintf('powsource_%02.f_%s', subj, condname);
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
  
  for f = 1:numel(powpeak)
    
    fprintf('\n   ->->-> Running % 2.f powpeak (%s) <-<-<-\n', f, powpeak(f).name);
    
    %-----------------%
    %-more precise source reconstruction
    %-------%
    %-check that baseline contains always data, otherwise shrink the time
    % window
    timelim = powpeak(f).wndw;
    begbline = cfg.powsource.bline - timelim/2; % beginning of baseline
    begtrl = data.time{1}(1); % beginning of trial
    if begbline < begtrl
      output = sprintf('%sPowpeak % 2.f (%s), window length was too long (% 3.2fs)\n', ...
        output, f, powpeak(f).name, timelim);
      timelim = (begtrl - cfg.powsource.bline) * -2;
    end
    %-------%
    
    %-------%
    %-check that time window is long enough for dpss smoothing
    % you get 0.75 from dpss(a,b). If b < 0.75, you get one taper (but the
    % last one is not used by fieldtrip). If it's between 0.75 and 1.25,
    % you get two tapers, but the second is not used and you get a warning
    % from fieldtrip. Still, the one taper in use is not a simple hanning
    tapsmofrq = powpeak(f).band/2;
    if timelim * tapsmofrq > 0.75 
      
      isdpss = true;
      output = sprintf('%sPowpeak % 2.f (%s), window length % 3.2fs, using dpss (tapsmofrq: +-% 2.1f Hz)\n', ...
        output, f, powpeak(f).name, timelim, tapsmofrq);
      
    else
      isdpss = false; % hanning, no smoothing
      output = sprintf('%sPowpeak % 2.f (%s), window length % 3.2fs, using hanning\n', ...
        output, f, powpeak(f).name, timelim);
      
    end
    %-------%
    
    cfg1 = [];
    cfg1.method = 'mtmconvol';
    cfg1.output = 'fourier';
    cfg1.foi = powpeak(f).freq;
    
    if isdpss
      cfg1.taper = 'dpss';
      cfg1.tapsmofrq = powpeak(f).band/2;
    else
      cfg1.taper = 'hanning';
    end
    
    cfg1.t_ftimwin = timelim;
    cfg1.toi = [cfg.powsource.bline powpeak(f).time];
    cfg1.feedback = 'none';
    cfg1.channel = datachan;
    freq = ft_freqanalysis(cfg1, data);
    %-----------------%
    
    %-----------------%
    %-source analysis
    cfg1 = [];
    cfg1.latency = cfg.powsource.bline;
    cfg1.frequency = powpeak(f).freq;
    
    cfg1.method = 'dics';
    cfg1.dics.feedback = 'none';
    cfg1.dics.lambda = cfg.powsource.lambda;
    cfg1.dics.powmethod = cfg.powsource.powmethod;
    
    cfg1.vol          = vol;
    cfg1.grid         = leadchan;
    cfg1.elec         = sens;
    
    souPre{f}       = ft_sourceanalysis(cfg1, freq);

    %-keep filters only for source analysis
    cfg1.(cfg1.method).keepfilter   = 'yes';
    cfg1.(cfg1.method).realfilter   = 'yes';

    cfg1.latency    = powpeak(f).time;
    source{f}       = ft_sourceanalysis(cfg1, freq);
    %-----------------%
    
  end

  %-----------------%
  %-save source
  save([cfg.dpow outputfile], 'source', 'souPre', '-v7.3')
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

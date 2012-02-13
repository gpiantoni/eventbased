function powsource_subj(cfg, subj)
%POWSOURCE_SUBJ DICS on interesting parts, defined by powpeaks

mversion = 14;
%14 12/03/04 calculate max length of timelim based on baseline and beginning of trial
%13 12/02/03 keep track of bad channels, use only good channels
%12 12/02/02 renamed to powsource_subj
%11 12/01/18 fixed small bug - -> + (begtrl = cfg.def.trialdef.prestim + cfg.powsource.bline;)
%10 12/01/13 baseline can change
%09 12/01/13 add cfg.pow.t_ftimwin to powpeak(f).wndw (which can be 0)
%08 12/01/12 check that time window is long enough for dpss smoothing, otherwise hanning
%07 12/01/12 use template or subject-specific vol
%06 12/01/12 tfr using improved methods (dpss, better time window)
%05 12/01/12 save both souPre and source
%04 12/01/11 cfg.powsource.source -> cfg.powsource.powpeak
%03 12/01/10 use predefined or power-peaks for areas of interest
%02 12/01/10 save real filters
%01 11/12/06 created from gosdpow

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
for e = 1:numel(cfg.erpeffect)
  k = cfg.erpeffect(e);
  
  %-----------------%
  %-input and output for each condition
  allfile = dir([ddir '*' cfg.test{k} cfg.endname '.mat']); % files matching a preprocessing
  if isempty(allfile)
    continue
  end
  
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
  
  for f = 1:numel(powpeak)
    
    %-----------------%
    %-more precise source reconstruction
    %-------%
    %-check that baseline contains always data, otherwise shrink the time
    % window
    % This time window is extremely tricky to work with: powpeak is
    % selected based on the significant time points of the TFR. However,
    % each of those points represents an FFT around that time window, with
    % length defined by cfg.pow.t_ftimwin. powpeak(f).wndw reports only the
    % significant timepoints (it can 0, if only one time point is
    % significant), so we add cfg.pow.t_ftimwin here (half in the beginning
    % and half in the end, but it's centered around powpeak(f).time anyway)
    ifreq = nearest(cfg.pow.foi, powpeak(f).freq); % length of cfg.pow.t_ftimwin depends on the frequency
    timelim = powpeak(f).wndw + cfg.pow.t_ftimwin(ifreq); 
    
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
    tapsmofrq = powpeak(f).band/2;
    if timelim > 1/tapsmofrq 
      isdpss = true;
      output = sprintf('%sPowpeak % 2.f (%s), window length % 3.2fs, using dpss (tapsmofrq: +-% 2.1f Hz)\n', ...
        output, f, powpeak(f).name, timelim, tapsmofrq);
      
    else % when ==, then dpss uses only one taper, hanning (but use hanning, otherwise we get one warning per trial, useless)
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

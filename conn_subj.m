function conn_subj(cfg, subj)
%SUBJ CONN connectivity on single-subject data

mversion = 15;
%15 12/02/02 renamed to conn_subj
%14 12/01/18 do transpose of cca, to keep consistent with fieldtrip
%13 12/01/15 use cca for time-domain granger
%12 12/01/12 toi is implemented by subfunctions
%11 12/01/11 toi looks at changes over time
%10 12/01/11 added source2mont, with possibility for fixed mom
%09 11/12/06 give feedback on topodipole (TODO: normalize)
%08 11/12/02 added topodipole
%07 11/12/01 do not use cfg.pow but simple mtmfft
%06 11/11/19 redefine trials and apply montage
%05 11/10/04 allows different types of connectivity methods
%04 11/09/27 cfg.conn.cond -> cfg.test
%03 11/09/12 2nd argument for subj (and cfg.subj -> subj)
%02 11/05/20 datdir -> ddir
%01 11/05/16 created

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
% %-----------------%
% %-add toolbox
% if strcmp(cfg.conn.toolbox, 'bsmart')
%   addpath /usr/local/toolbox/bsmart/
% else
%   addpath /usr/local/toolbox/biosig/t200_FileAccess/
%   addpath /usr/local/toolbox/biosig/tsa/
% end
% %-----------------%

ddir = sprintf('%s%04.f/%s/%s/', cfg.data, subj, cfg.mod, cfg.cond); % data

%-----------------%
%-prepare montage
if strcmp(cfg.conn.areas, 'channel')
  
  if ischar(cfg.seldata.channel)
    error('You need to specify all the channels in cfg.prepr.channel');
  else
    mont = prepare_montage(cfg.conn.chan, cfg.seldata.channel');
  end
  
elseif strcmp(cfg.conn.areas, 'dipole')
  
  %-------%
  %-read data
  if ~exist([cfg.derp cfg.proj '_granderp.mat'], 'file')
    error([cfg.derp cfg.proj '_granderp.mat does not exist'])
  end
  
  load([cfg.derp cfg.proj '_granderp'], 'gerp')
  %-------%
  
  mont = topodipole(cfg.conn.dip, gerp{cfg.erpeffect});
  
elseif strcmp(cfg.conn.areas, 'erppeak') || strcmp(cfg.conn.areas, 'powpeak')
  
  if strcmp(cfg.conn.areas, 'erppeak')
    
    %-------%
    %-load sources of erp
    condname = regexprep(cfg.test{cfg.erpeffect}, '*', '');
    inputfile = sprintf('erpsource_%02.f_%s', subj, condname);
    load([cfg.derp inputfile], 'source')
    
    load([cfg.derp cfg.proj '_soupeak'], 'soupeak')
    %-------%
    
  elseif strcmp(cfg.conn.areas, 'powpeak')
    
    %-------%
    %-load source of pow
    condname = regexprep(cfg.test{cfg.poweffect}, '*', '');
    inputfile = sprintf('powsource_%02.f_%s', subj, condname);
    load([cfg.dpow inputfile], 'source')
    
    load([cfg.dpow cfg.proj '_soupeak'], 'soupeak')
    %-------%
    
  end
  
  %-------%
  %-calculate montage
  [mont outtmp] = source2mont(source, soupeak);
  output = [output outtmp];
  %-------%
  
end
%-----------------%

%-----------------%
%-plot feedback on montage
figure
imagesc(mont.tra, [-1 1] * max(abs(mont.tra(:))))
colorbar

pngfile = [cfg.log filesep cfg.proj '_gosdconn_mont' sprintf('%03.f', subj) '.png'];
saveas(gcf, pngfile);
close(gcf); drawnow
%-----------------%
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
  outputfile = sprintf('%s_%s_%02.f_%s', cfg.proj, cfg.conn.method, subj, condname);
  %-----------------%
  
  %-----------------%
  %-concatenate
  spcell = @(name) sprintf('%s%s', ddir, name);
  allname = cellfun(spcell, {allfile.name}, 'uni', 0);
  
  cfg1 = [];
  cfg1.inputfile = allname;
  data = ft_appenddata(cfg1);
  %-----------------%
  
  %-----------------%
  %-apply montage
  data = ft_apply_montage(data, mont, 'feedback', 'none');
  
  if strcmp(cfg.conn.areas, 'erppeak') || strcmp(cfg.conn.areas, 'powpeak')
    data = pcadata(data, soupeak, cfg.conn.fixedmom);
  end
  %-----------------%
  
  switch cfg.conn.type
    case 'ft'
      %---------------------------%
      %-FIELDTRIP calculate model
      %-----------------%
      %-mvar
      if strcmpi(cfg.conn.mvar, 'yes')
        cfg2 = [];
        cfg2.order = cfg.conn.order;
        cfg2.toolbox = cfg.conn.toolbox;
        cfg2.feedback = 'none';
        
        cfg2.toi = cfg.conn.toi;
        cfg2.t_ftimwin = cfg.conn.t_ftimwin;
        
        data = ft_mvaranalysis(cfg2, data); % ft_mvaranalysis can do it on multiple time points, but freqanalysis does not handle it anymore
      end
      %-----------------%
      
      %-----------------%
      %-freq on mvar
      if strcmpi(cfg.conn.freq, 'yes')
        
        if strcmpi(cfg.conn.mvar, 'yes')
          
          %--------%
          %-use special freq analysis for mvar data
          cfg3 = [];
          cfg3.method    = 'mvar';
          
          data = ft_freqanalysis_mvar(cfg3, data);
          %--------%
          
        else
          
          %--------%
          %-no mvar
          cfg3 = [];
          cfg3.method = 'mtmconvol';
          cfg3.taper = 'hanning';
          cfg3.foi  = cfg.pow.foi;
          cfg3.output = 'fourier';
          cfg3.feedback = 'none';
          cfg3.toi = cfg.conn.toi;
          cfg3.t_ftimwin = cfg.conn.t_ftimwin .* ones(numel(cfg3.foi));
          
          data = ft_freqanalysis(cfg3, data);
          %--------%
          
        end
        
      end
      %-----------------%
      
      %-----------------%
      %-connectivity
      cfg4         = [];
      cfg4.method  = cfg.conn.method;
      cfg4.channelcmb = {'all', 'all'};
      cfg4.outputfile = [cfg.dcon outputfile];
      ft_connectivityanalysis(cfg4, data);
      %-----------------%
      %---------------------------%
      
    case 'cca'
      
      addpath(genpath('/data1/toolbox/gcca'))
      
      %---------------------------%
      %-CCA calculate model
      nchan = numel(data.label);
      gcmat = NaN(nchan, nchan, 1, numel(cfg.conn.toi));
      
      for t = 1:numel(cfg.conn.toi)
        
        %-----------------%
        %-into cca format
        dat = ft_selectdata(data, 'toilim', cfg.conn.toi(t)+[-.5 .5]*cfg.conn.t_ftimwin);
        X = [dat.trial{:}]; % after ft_redefinetrial
        Nr = numel(dat.trial);
        Nl = numel(dat.time{1});
        %-----------------%
        
        %-----------------%
        %- granger for each time window
        X = cca_detrend_mtrial(X, Nr, Nl);
        X = cca_rm_temporalmean_mtrial(X, Nr, Nl);
        
        ret = cca_granger_regress_mtrial(X, Nr, Nl, cfg.conn.order);
        gcmat(:,:,1,t) = transpose(ret.gc);
        %-----------------%
        
      end
      
      %-----------------%
      %-convert into fieldtrip format
      stat = [];
      stat.dimord = 'chan_chan_freq_time';
      stat.label = data.label;
      stat.gc = gcmat;
      stat.freq = Inf; % only one frequency
      stat.time = cfg.conn.toi;
      stat.cfg.previous = data.cfg;
      stat.cfg.method = 'cca';
      stat.cfg.order = cfg.conn.order;
      %-----------------%
      
      save([cfg.dcon outputfile], 'stat')
      %---------------------------%
  end
  
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


%---------------------------------------------------------%
%-SUBFUNCTION: pcadata
%---------------------------------------------------------%
function [data] = pcadata(data, soupeak, fixedmom)
%PCADATA simplify data using pca on each region of interest
% keep only the first component

%-------------------------------------%
%-loop over regions of interest
trial = [];
for i1 = 1:numel(soupeak)
  
  %-----------------%
  %-find channels belonging to region of interest
  newname = sprintf('%s', soupeak(i1).name);
  label{i1,1} = newname;
  ivox = ~cellfun(@isempty, strfind(data.label, newname));
  %-----------------%
  
  if strcmp(fixedmom, 'no')
    %-----------------%
    %-the moment is different in each trial
    for t = 1:numel(data.trial)
      [~, pcatrl] = princomp(data.trial{t}(ivox,:)');
      trial{t}(i1,:) = pcatrl(:,1)';
    end
    %-----------------%
    
  else
    
    %-----------------%
    %-get coefficient for all the trials, then apply it to the data
    alltrl = [data.trial{:}];
    [coef] = princomp(alltrl(ivox,:)');
    
    for t = 1:numel(data.trial)
      [pcatrl] = data.trial{t}(ivox,:)' * coef;
      trial{t}(i1,:) = pcatrl(:,1)';
    end
    %-----------------%
  end
  
end
%-------------------------------------%

data.trial = trial;
data.label = label;
%---------------------------------------------------------%

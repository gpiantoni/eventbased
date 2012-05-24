function conn_subj(cfg, subj)
%CONN_SUBJ connectivity on single-subject data
%
% CFG
%  .data: path of /data1/projects/PROJ/subjects/
%  .rec: REC in /data1/projects/PROJ/recordings/REC/
%  .nick: NICK in /data1/projects/PROJ/subjects/0001/MOD/NICK/
%  .mod: modality, MOD in /data1/projects/PROJ/subjects/0001/MOD/NICK/
%  .endname: includes preprocessing steps (e.g. '_seldata_gclean_redef')
%
%  .log: name of the file and directory to save log
%  .dconn: directory for connectivity data
%  .conn.cond: cell with conditions (e.g. {'*cond1' '*cond2'})'
%
%-ROI parameters
%  .conn.areas: 'channel' or 'dipole' or 'erppeak' or 'powpeak'
%    if 'erppeak' (use beamformer to construct virtual electrode):
%      .derp: directory with ERP data
%      .erpsource.refcond: string of the condition used for source location
%      .conn.fixedmom: logical (use the same moment for source or change it every time)
%    if 'powpeak' (use beamformer to construct virtual electrode):
%      .dpow: directory with POW data
%      .powsource.refcond: string of the condition used for source location
%      .conn.fixedmom: logical (use the same moment for source or change it every time)
%
%-Connectivity parameters
%  .conn.toi: vector with time points to run connectivity on
%  .conn.t_ftimwin: scalar with duration of time window
%  .conn.type: 'cca' or 'ft'
%    if 'cca' (use Anil Seth's implementation):
%      .conn.method: 'gc' (only this method is available, it's time-domain based) 
%      .conn.order: scalar indicating model order or 'aic' or 'bic'
%                   If 'aic' or 'bic', it estimates the model order and
%                   writes it down estimated model order in csv in log dir.
%
%    if 'ft' (use Fieldtrip implementation):
%      .conn.method: can be symmetric or asymmetric (a string of one of the below)
%SYM:  'coh', 'csd', 'plv', 'powcorr', 'amplcorr', 'ppc', 'psi'
%ASYM: 'granger', 'dtf', 'pdc'
%      .conn.freq: logical, needs to be TRUE (in future, Fieldtrip might implement time-domain connectivity measures)
%      .conn.mvar: logical, estimate coefficients or not
%        if TRUE:
%          .conn.toolbox: 'biosig' or 'bsmart'
%          .conn.order: scalar indicating model order
%        if FALSE:
%          .conn.foi: frequency of interest (best if identical to .pow.foi)
%
% IN:
%  data in /PROJ/subjects/SUBJ/MOD/NICK/
%  if .conn.areas == 'erppeak'
%     [cfg.derp 'erpsource_SUBJ_COND']: source data for period of interest for each subject
%     [cfg.derp 'NICK_COND_soupeak']: significant source peaks in the ERP
%  if .conn.areas == 'powpeak'
%     [cfg.dpow 'powsource_SUBJ_COND']: source data for period of interest for each subject
%     [cfg.dpow 'NICK_COND_soupeak']: significant source peaks in the POW
%
% OUT
%  [cfg.dcon 'conn_CONNMETHOD_SUBJ_COND']: connectivty analysis for each subject
%  If .conn.type is 'cca', it also returns the model check. The file has
%    subj, condition, toi, ADF, KPSS, residuals, consistency
%  where 
%    ADF, KPSS check the covariance stationary (should be 0, can give different results)
%    residuals tells you if the residuals, after fitting GC, are not white
%    consistency gives you how much variance the model explains (better if > 80)
%
% Part of EVENTBASED single-subject
% see also ERP_SUBJ, ERP_GRAND, 
% ERPSOURCE_SUBJ, ERPSOURCE_GRAND, ERPSTAT_SUBJ, ERPSTAT_GRAND,
% POW_SUBJ, POW_GRAND, POWCORR_SUBJ, POWCORR_GRAND,
% POWSOURCE_SUBJ, POWSOURCE_GRAND, POWSTAT_SUBJ, POWSTAT_GRAND,
% CONN_SUBJ, CONN_GRAND, CONN_STAT


% 
%  .conn.areas: 'channel' or 'dipole' or 'erppeak' or 'powpeak'
%    if 'channel' (take average over groups of channels)
%      .conn.chan: a struct with: 
%                               .name: 'name of group elec'
%                               .chan: {'elec1' 'elec2' 'elec3'};
%      .seldata.label: names of the channels to match with cfg.conn.chan
%      (not very robust, might require more testing)
%   
%    if 'dipole' (use weights from erp at specific time window):
%      .derp: directory with ERP data
%      .erpeffect: vector of condition whose ERP is used here
%      .conn.dip: a struct with:
%                              .name: 'name of dipole'
%                              .time: one or two scalars, time of the ERP of interest
%   
    


%---------------------------%
%-start log
output = sprintf('%s (%04d) began at %s on %s\n', ...
  mfilename, subj, datestr(now, 'HH:MM:SS'), datestr(now, 'dd-mmm-yy'));
tic_t = tic;
%---------------------------%

%---------------------------%
%-dir and files

%-----------------%
%-prepare montage
if strcmp(cfg.conn.areas, 'channel')
  % TODO
  %   if ischar(cfg.seldata.label)
  %     error('You need to specify all the channels in cfg.seldata.label');
  %   else
  %     mont = prepare_montage(cfg.conn.chan, cfg.seldata.label');
  %   end
  
elseif strcmp(cfg.conn.areas, 'dipole')
  % TODO
  %   %-------%
  %   %-read data
  %   if ~exist([cfg.derp cfg.cond '_granderp.mat'], 'file')
  %     error([cfg.derp cfg.cond '_granderp.mat does not exist'])
  %   end
  %
  %   load([cfg.derp cfg.cond '_granderp'], 'gerp')
  %   %-------%
  %
  %   mont = topodipole(cfg.conn.dip, gerp{cfg.erpeffect(1)});
  
elseif strcmp(cfg.conn.areas, 'erppeak') || strcmp(cfg.conn.areas, 'powpeak')
  
  if strcmp(cfg.conn.areas, 'erppeak')
    
    %-------%
    %-load source of pow
    condname = regexprep(cfg.erpsource.refcond, '*', ''); % DOC: cfg.erpsource.refcond
    inputfile = sprintf('erpsource_%04d_%s', subj, condname);

    load([cfg.derp inputfile], 'source')
    load([cfg.derp cfg.nick '_' condname '_soupeak'], 'soupeak')
    %-------%
    
  elseif strcmp(cfg.conn.areas, 'powpeak')
    
    %-------%
    %-load source of pow
    condname = regexprep(cfg.powsource.refcond, '*', ''); % DOC: cfg.powsource.refcond
    inputfile = sprintf('powsource_%04d_%s', subj, condname);

    load([cfg.dpow inputfile], 'source') 
    load([cfg.dpow cfg.nick '_' condname '_soupeak'], 'soupeak') % DOC: in
    %-------%
    
  end
  
  %-------%
  %-calculate montage
  [mont outtmp] = source2mont(source, soupeak);
  output = [output outtmp];
  %-------%
  
end
%-----------------%
%---------------------------%

%-------------------------------------%
%-loop over conditions
for k = 1:numel(cfg.conn.cond)
  cond     = cfg.conn.cond{k};
  condname = regexprep(cond, '*', '');
  
  %---------------------------%
  %-read data
  [data] = load_data(cfg, subj, cond);
  if isempty(data)
    output = sprintf('%sCould not find any file for condition %s\n', ...
      output, cond);
    continue
  end
  
  outputfile = sprintf('conn_%s_%04d_%s', cfg.conn.method, subj, condname);
  %---------------------------%
  
  %-----------------%
  %-apply montage
  data = ft_apply_montage(data, mont, 'feedback', 'none');
  
  if strcmp(cfg.conn.areas, 'erppeak') || strcmp(cfg.conn.areas, 'powpeak')
    data = pcadata(data, soupeak, cfg.conn.fixedmom);
  end
  
  data = ft_checkdata(data, 'hassampleinfo', 'yes'); % recreate sampleinfo which got lost (necessary for selfromraw)
  %-----------------%
  
  switch cfg.conn.type
    case 'ft'
      %---------------------------%
      %-FIELDTRIP calculate model
      %-----------------%
      %-mvar
      if cfg.conn.mvar
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
      if cfg.conn.freq
        
        if cfg.conn.mvar
          
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
          cfg3.foi  = cfg.conn.foi;
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
        X = cca_rm_ensemblemean(X, Nr, Nl);
        
        %-------%
        %-model order
        if ~ischar(cfg.conn.order)
          order = cfg.conn.order;
          
        else % estimate model order
          [bic, aic] = cca_find_model_order_mtrial(X,Nr, Nl, 2, 20);
          
          if strcmpi(cfg.conn.order, 'bic')
            order = bic;
          elseif strcmpi(cfg.conn.order, 'aic')
            order = aic;
          end
          
          %-output
          output = sprintf('%s At time % 6.2f, model order estimated with %s is% 3d\n', ...
            output, cfg.conn.toi(t), cfg.conn.order, order);
          
          csvfile = sprintf('%s/cca_modelorder_%s.csv', cfg.log, cfg.conn.order);
          csvtext = sprintf('''%04d'',''%s'',%1.2f,%d\n', ...
            subj, condname, cfg.conn.toi(t), order);
          
          fid = fopen(csvfile, 'a+');
          fwrite(fid, csvtext);
          fclose(fid);
          
        end
        %-------%
        
        ret = cca_granger_regress_mtrial(X, Nr, Nl, order);
        gcmat(:,:,1,t) = transpose(ret.gc);
        %-----------------%
        
        %-----------------%
        %-check the model
        %-------%
        %-covariance ADF
        % is covariance non-stationary? should be 0
        adf = cca_check_cov_stat_mtrial(X, Nr, Nl, 10);
        adfr = numel(find(adf)) / numel(adf);
        %-------%
        
        %-------%
        %-covariance KPSS
        % is covariance non-stationary? should be 0
        [kh,kpss] = cca_kpss_mtrial(X, Nr, Nl);
        khr = numel(find(kh==0)) / numel(kh);
        %-------%
        
        %-------%
        %-white residuals
        % ratio of significant residuals, should be 0
        nvar = size(X,1);
        dwthresh = 0.05/nvar;    % critical threshold, Bonferroni corrected
        dw = ret.waut < dwthresh;
        dwr = numel(find(dw)) / numel(dw);
        %-------%
        
        %-------%
        %-model consistency
        % check model consistency, ie. proportion of correlation structure of the
        % data accounted for by the MVAR model
        % should be more than 80
        %-------%
        
        %-------%
        %-write to file
        % subj, condition, toi, ADF, KPSS, white residuals, consistency
        csvfile = sprintf('%s/cca_modelcheck.csv', cfg.log);
        csvtext = sprintf('''%04d'',''%s'',%1.2f,%1.2f,%1.2f,%1.2f,%1.2f\n', ...
          subj, condname, cfg.conn.toi(t), adfr, khr, dwr, ret.cons);
        
        fid = fopen(csvfile, 'a+');
        fwrite(fid, csvtext);
        fclose(fid);
        %-------%
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
outtmp = sprintf('%s (%04d) ended at %s on %s after %s\n\n', ...
  mfilename, subj, datestr(now, 'HH:MM:SS'), datestr(now, 'dd-mmm-yy'), ...
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
  
  if ~fixedmom
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

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
%  .conn.cond: cell with conditions (e.g. {'*cond1' '*cond2'})
%  .conn.source: read virtual electrode data (logical)
%
%-Connectivity parameters
%  .conn.type: 'cca' or 'ft'
%    if 'cca' (use Anil Seth's implementation):
%      .conn.method: 'gc' (only this method is available, it's time-domain based)
%      .conn.order: scalar indicating model order or 'aic' or 'bic'
%                   If 'aic' or 'bic', it estimates the model order and
%                   writes it down estimated model order in csv in log dir.
%      .conn.toi: vector with time points to run connectivity on
%      .conn.t_ftimwin: scalar with duration of time window
%
%    if 'ft' (use Fieldtrip implementation):
%      .conn.method: can be symmetric or asymmetric (a string of one of the below)
%SYM:  'coh', 'csd', 'plv', 'powcorr', 'amplcorr', 'ppc', 'psi'
%ASYM: 'granger', 'dtf', 'pdc'
%      .conn.mvar: logical, estimate coefficients or not
%        if TRUE:
%          .conn.toolbox: 'biosig' or 'bsmart'
%          .conn.order: scalar indicating model order
%          .conn.toi: vector with time points to run connectivity on
%          .conn.t_ftimwin: scalar with duration of time window
%        if FALSE:
%          .conn.foi: frequency of interest (best if identical to .pow.foi)
%      .conn.freq: 'mtmconvol' or 'mtmfft'
%        in either case,
%          .conn.avgoverfreq: average power spectrum
%          .conn.planar for MEG if you want to run planar
%          (but does it make sense to do planar on fourier data?)
%        if 'mtmconvol':
%          .conn.foi: frequency of interest
%          .conn.toi: vector with time points to run connectivity on
%          .conn.t_ftimwin: scalar with duration of time window (same length as .conn.foi)
%        if 'mtmfft':
%          .conn.foilim: two values for the frequency of interest
%
% IN:
%  data in /PROJ/subjects/SUBJ/MOD/NICK/
%  if .conn.areas == 'erppeak' or ('dip' and .conn.beamformer == 'erp')
%     [cfg.derp 'erpsource_SUBJ_COND']: source data for period of interest for each subject
%     [cfg.derp 'NICK_COND_soupeak']: significant source peaks in the ERP
%  if .conn.areas == 'powpeak'  or ('dip' and .conn.beamformer == 'pow')
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
% SOURCE_SUBJ, CONN_SUBJ, CONN_GRAND, CONN_STAT

%---------------------------%
%-start log
output = sprintf('%s (%04d) began at %s on %s\n', ...
  mfilename, subj, datestr(now, 'HH:MM:SS'), datestr(now, 'dd-mmm-yy'));
tic_t = tic;
%---------------------------%

%-------------------------------------%
%-loop over conditions
for k = 1:numel(cfg.conn.cond)
  cond     = cfg.conn.cond{k};
  condname = regexprep(cond, '*', '');
  
  %---------------------------%
  %-read data
  if ~isfield(cfg.conn, 'source') || ~cfg.conn.source
    [data] = load_data(cfg, subj, cond);
  else
    [data] = load_source(cfg, subj, cond);
  end
  if isempty(data)
    output = sprintf('%sCould not find any file for condition %s\n', ...
      output, cond);
    continue
  end
  
  outputfile = sprintf('conn_%04d_%s', subj, condname);
  %---------------------------%
  
  switch cfg.conn.type
    
    case 'ft'
      
      %---------------------------%
      %-Fieltrip
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
      %---------------------------%
      
      %---------------------------%
      %-fieldtrip way or use fieldtrip function on mvar of space-state
      if cfg.conn.mvar
        
        %-----------------%
        %-freq on mvar
        %--------%
        %-use special freq analysis for mvar data
        tmpcfg = [];
        tmpcfg.method    = 'mvar';
        
        data = ft_freqanalysis(tmpcfg, data);
        %--------%
        %-----------------%
        
      else
        
        %-----------------%
        %-freq on raw data
        
        %-------%
        %-planar
        if isfield(data, 'grad') && cfg.conn.planar
          
          tmpcfg = [];
          tmpcfg.grad = data.grad;
          tmpcfg.method = 'distance';
          tmpcfg.neighbourdist = cfg.sens.dist;
          nbor = ft_prepare_neighbours(tmpcfg);
          
          tmpcfg = [];
          tmpcfg.neighbours = nbor;
          data = ft_megplanar(tmpcfg, data);
          
        end
        %-------%
        
        %-------%
        %-no mvar
        tmpcfg = [];
        
        tmpcfg.taper = 'hanning';
        tmpcfg.feedback = 'none';
        tmpcfg.keeptrials = 'yes';
        
        if strcmp(cfg.conn.freq, 'mtmconvol')
          tmpcfg.method = 'mtmconvol';
          tmpcfg.toi = cfg.conn.toi;
          tmpcfg.foi = cfg.conn.foi;
          tmpcfg.t_ftimwin = cfg.conn.t_ftimwin .* ones(numel(tmpcfg.foi));
          
        elseif strcmp(cfg.conn.freq, 'mtmfft')
          tmpcfg.method = 'mtmfft';
          tmpcfg.foilim = cfg.conn.foilim;
          
        end
        
        switch cfg.conn.method
          case {'granger'}
            tmpcfg.output = 'fourier';
          case {'powcorr'}
            tmpcfg.output = 'pow';
        end
        
        data = ft_freqanalysis(tmpcfg, data);
        %-------%
        
        %-------%
        %-planar
        if isfield(data, 'grad') && cfg.conn.planar
          tmpcfg = [];
          data = ft_combineplanar(tmpcfg, data);
        end
        %-------%
        
        %-------%
        %-average over frequency
        if isfield(cfg.conn, 'avgoverfreq') && cfg.conn.avgoverfreq
          data = ft_selectdata(data, 'avgoverfreq', 'yes');
        end
        %-------%
        %-----------------%
        
      end
      
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

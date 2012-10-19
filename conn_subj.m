function conn_subj(info, opt, subj)
%CONN_SUBJ connectivity on single-subject data
%
% INFO
%  .dcon: directory with CONN data
%  .log: name of the file and directory to save log
% 
% CFG.OPT
%  .source: read virtual electrode data (logical)
%  .cond*: cell with conditions (e.g. {'*cond1' '*cond2'})
%
%-Connectivity parameters
%  .conn.type*: 'cca' or 'ft'
%    if 'cca' (use Anil Seth's implementation):
%      .conn.method*: 'gc' (only this method is available, it's time-domain based)
%      .conn.order*: scalar indicating model order or 'aic' or 'bic'
%                   If 'aic' or 'bic', it estimates the model order and
%                   writes it down estimated model order in csv in log dir.
%      .conn.toi*: vector with time points to run connectivity on
%      .conn.t_ftimwin*: scalar with duration of time window
%
%    if 'ft' (use Fieldtrip implementation):
%      .conn.method: 'coh', 'csd', 'plv', 'powcorr', 'amplcorr', 'ppc',
%      'psi' (symmetric) or 'granger', 'dtf', 'pdc' (asymmetric)
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
% * indicates obligatory parameter
%
% IN
%  LOAD_DATA: data in /PROJ/subjects/SUBJ/MOD/NICK/
% OR if cfg.opt.source
%  LOAD_SOURCE: source in info.dsou after SOURCE_SUBJ
%
% OUT
%  [info.dcon 'conn_CONNMETHOD_SUBJ_COND'] 'conn_s' connectivty analysis for each subject
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
for k = 1:numel(opt.cond)
  cond = opt.cond{k};
  condname = regexprep(cond, '*', '');
  
  %---------------------------%
  %-read data
  if ~isfield(opt, 'source') || ~opt.source
    [data] = load_data(info, subj, cond);
  else
    [data] = load_source(info, subj, cond);
  end
  if isempty(data)
    output = sprintf('%sCould not find any file for condition %s\n', ...
      output, cond);
    continue
  end
  
  outputfile = sprintf('conn_%04d_%s', subj, condname);
  %---------------------------%
  
  switch opt.conn.type
    
    case 'ft'

      if opt.conn.mvar
        
        %-----------------%
        %-mvar
        cfg = [];
        cfg.order = opt.conn.order;
        cfg.toolbox = opt.conn.toolbox;
        cfg.feedback = 'none';
        
        cfg.toi = opt.conn.toi;
        cfg.t_ftimwin = opt.conn.t_ftimwin;
        
        data = ft_mvaranalysis(cfg, data); % ft_mvaranalysis can do it on multiple time points, but freqanalysis does not handle it anymore
        %-----------------%
        
        %-----------------%
        %-freq on mvar
        cfg = [];
        cfg.method    = 'mvar';
        
        data = ft_freqanalysis(cfg, data);
        %-----------------%
        
      else
        
        %-----------------%
        %-freq on raw data
        
        %-------%
        %-planar
        if isfield(data, 'grad') && opt.conn.planar
          
          cfg = [];
          cfg.grad = data.grad;
          cfg.method = 'distance';
          cfg.neighbourdist = info.sens.dist;
          nbor = ft_prepare_neighbours(cfg);
          
          cfg = [];
          cfg.neighbours = nbor;
          data = ft_megplanar(cfg, data);
          
        end
        %-------%
        
        %-------%
        %-no mvar
        cfg = [];
        
        cfg.taper = 'hamming'; % 'hanning'
        cfg.feedback = 'none';
        cfg.keeptrials = 'yes';
        
        if strcmp(opt.conn.freq, 'mtmconvol') % multiple time window, but slow
          cfg.method = 'mtmconvol';
          cfg.toi = opt.conn.toi;
          cfg.foi = opt.conn.foi;
          cfg.t_ftimwin = opt.conn.t_ftimwin .* ones(size(cfg.foi));
          
        elseif strcmp(opt.conn.freq, 'mtmfft') % one time window, faster
          cfg.method = 'mtmfft';
          cfg.foilim = [min(opt.conn.foilim(:,1)) max(opt.conn.foilim(:,2))]; % doc
          
        end
        
        switch opt.conn.method
          case {'granger'}
            cfg.output = 'fourier';
          case {'powcorr'}
            cfg.output = 'pow';
        end
        
        data = ft_freqanalysis(cfg, data);
        %-------%
        
        %-------%
        %-planar
        if isfield(data, 'grad') && opt.conn.planar
          cfg = [];
          data = ft_combineplanar(cfg, data);
        end
        %-------%
        
        %-------%
        %-average over frequency
        % it doesn't average for mvar
        if isfield(opt.conn, 'avgoverfreq') && opt.conn.avgoverfreq

          if isfield(opt.conn, 'foilim') && size(opt.conn.foilim,1) > 1
            clear dat
            for i = 1:size(opt.conn.foilim,1)
              dat(i) = ft_selectdata(data, 'avgoverfreq', 'yes', 'foilim', opt.conn.foilim(i,:));
            end
            
            data = dat(1);
            data.freq = cat(2, dat.freq);
            data.powspctrm = cat(3, dat.powspctrm);
            
          else
            data = ft_selectdata(data, 'avgoverfreq', 'yes');
          end
          
        end
        %-------%
        %-----------------%
        
      end
      
      %-----------------%
      %-connectivity
      cfg = [];
      cfg.method  = opt.conn.method;
      cfg.channelcmb = {'all', 'all'};
      conn_s = ft_connectivityanalysis(cfg, data);
      save([info.dcon outputfile], 'conn_s')
      %-----------------%
      %---------------------------%
      
    case 'cca'
      
      addpath(genpath('/data1/toolbox/gcca'))
      
      %---------------------------%
      %-CCA calculate model
      nchan = numel(data.label);
      gcmat = NaN(nchan, nchan, 1, numel(opt.conn.toi));
      
      for t = 1:numel(opt.conn.toi)
        
        %-----------------%
        %-into cca format
        dat = ft_selectdata(data, 'toilim', opt.conn.toi(t)+[-.5 .5]*opt.conn.t_ftimwin);
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
        if ~ischar(opt.conn.order)
          order = opt.conn.order;
          
        else % estimate model order
          [bic, aic] = cca_find_model_order_mtrial(X,Nr, Nl, 2, 20);
          
          if strcmpi(opt.conn.order, 'bic')
            order = bic;
          elseif strcmpi(opt.conn.order, 'aic')
            order = aic;
          end
          
          %-output
          output = sprintf('%s At time % 6.2f, model order estimated with %s is% 3d\n', ...
            output, opt.conn.toi(t), opt.conn.order, order);
          
          csvfile = sprintf('%s/cca_modelorder_%s.csv', info.log, opt.conn.order);
          csvtext = sprintf('''%04d'',''%s'',%1.2f,%d\n', ...
            subj, condname, opt.conn.toi(t), order);
          
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
        csvfile = sprintf('%s/cca_modelcheck.csv', info.log);
        csvtext = sprintf('''%04d'',''%s'',%1.2f,%1.2f,%1.2f,%1.2f,%1.2f\n', ...
          subj, condname, opt.conn.toi(t), adfr, khr, dwr, ret.cons);
        
        fid = fopen(csvfile, 'a+');
        fwrite(fid, csvtext);
        fclose(fid);
        %-------%
        %-----------------%
        
      end
      
      %-----------------%
      %-convert into fieldtrip format
      conn_s = [];
      conn_s.dimord = 'chan_chan_freq_time';
      conn_s.label = data.label;
      conn_s.gc = gcmat;
      conn_s.freq = Inf; % only one frequency
      conn_s.time = opt.conn.toi;
      conn_s.cfg.previous = data.cfg;
      conn_s.cfg.method = 'cca';
      conn_s.cfg.order = opt.conn.order;
      %-----------------%
      
      save([info.dcon outputfile], 'conn_s')
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
fid = fopen([info.log '.txt'], 'a');
fwrite(fid, output);
fclose(fid);
%-----------------%
%---------------------------%

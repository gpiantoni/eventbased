function [dataout output] = conn_statespace(cfg, data)
%CONN_STATESPACE calculate MVAR model at source space
% using state-space models

addpath(fullfile(fileparts(mfilename('fullpath')), 'ssm_em'))

%-------------------------------------%
%-get the options
%-----------------%
%-data
cfg.channel      = ft_getopt(cfg, 'channel',  'all');
cfg.resamplefs   = ft_getopt(cfg, 'resamplefs',  []);

cfg.scaling      = ft_getopt(cfg, 'scaling', []);
cfg.scaling.data = ft_getopt(cfg.scaling, 'data', []);
cfg.scaling.lf   = ft_getopt(cfg.scaling, 'lf', []);
cfg.bases        = ft_getopt(cfg, 'bases', 1); % TODO
%-----------------%

%-----------------%
%-initiation parameters
cfg.order        = ft_getopt(cfg, 'order', 2);
cfg.A_range      = ft_getopt(cfg, 'A_range', [-0.5 0.5]);
cfg.Q_range      = ft_getopt(cfg, 'Q_range', [0.1 0.9]);
%-----------------%

%-----------------%
%-EM parameters
cfg.maxiter      = ft_getopt(cfg, 'maxiter', 20);
cfg.tol          = ft_getopt(cfg, 'tol', 1e-5);
%-----------------%

%-----------------%
%-remember random state
stream = RandStream.getGlobalStream;
output = sprintf('random stream: Type %s Seed %d\n', ...
  stream.Type, stream.Seed);
reset(stream); % go back to beginning of the seed (values are different because fexec initializes a new stream)
%-----------------%
%-------------------------------------%

%-------------------------------------%
%-prepare forward model
%-----------------%
%-channels
cfg.channel = ft_channelselection(cfg.channel, data.label);
data = ft_selectdata(data, 'channel', cfg.channel);
%-----------------%

%-----------------%
% cortical patches
tmpcfg = [];
tmpcfg.channel = cfg.channel;
tmpcfg.vol = cfg.vol;
tmpcfg.elec = cfg.elec; % this can be electrodes or gradiometers
tmpcfg.grid.pos = cat(1,cfg.roi.pos);
tmpcfg.grid.mom = cat(2,cfg.roi.mom);
grid = ft_prepare_leadfield(tmpcfg);
%-----------------%

%-----------------%
%-TODO: reduce with svd
lf = [grid.leadfield{:}];
%-----------------%

%-----------------%
% how many dipoles in each roi
roi = arrayfun( @(x) size(x.pos,1), cfg.roi);
nroi = numel(roi);
%-----------------%
%-------------------------------------%

%-------------------------------------%
%-----------------%
%-prepare data
if ~isempty(cfg.resamplefs)
  tmpcfg = [];
  tmpcfg.resamplefs = cfg.resamplefs;
  tmpcfg.detrend = 'yes';
  data = ft_resampledata(tmpcfg, data);
end
%-----------------%

%-----------------%
%-scaling to prevent numerical inconsistency
cat_data = [data.trial{:}];
output = [output sprintf('data units: % 10.4f, lead units: % 10.4f\n', std(cat_data(:)), std(lf(:)))];
if ~isempty(cfg.scaling)
  
  if ~isempty(cfg.scaling.lf)
    output = [output sprintf('applying scaling to lead: % 2d (magnitude)\n', log10(cfg.scaling.lf))];
    lf = lf * cfg.scaling.lf;
  end
  
  if ~isempty(cfg.scaling.data)
    output = [output sprintf('applying scaling to data: % 2d (magnitude)\n', log10(cfg.scaling.data))];
    for i = 1:numel(data.trial)
      data.trial{i} = data.trial{i} * cfg.scaling.data;
    end
  end
  cat_data = [data.trial{:}];
  output = [output sprintf('data units: % 10.4f, lead units: % 10.4f\n', std(cat_data(:)), std(lf(:)))];
end

clear cat_data
%-----------------%
%-------------------------------------%

%-------------------------------------%
%-compute EM
%-----------------%
%-covariance
tmpcfg = [];
tmpcfg.covariance = 'yes';
tmpcfg.covariancewindow = 'all';
tmpcfg.removemean = 'no'; % Cheung's algorithm does not remove the mean, but only rounding differences
tmpcfg.feedback = 'none';
timelock = ft_timelockanalysis(tmpcfg, data);
output = [output sprintf('number of channels % 3d, rank of covariance % 3d\n', ...
  size(timelock.cov,1), rank(timelock.cov))];
%-----------------%

%-----------------%
%-estimate activity in the sources, to initialize the values
xhat = xhat_lcmv(data.trial, lf, timelock.cov); 
%-----------------%

%-----------------%
%-initialize values
tmpcfg = [];
tmpcfg.roi = roi;
tmpcfg.order = cfg.order;
tmpcfg.A_range = cfg.A_range;
tmpcfg.Q_range = cfg.Q_range;
init = ssm_em_init(tmpcfg, xhat, lf, timelock.cov);
%-----------------%

%-----------------%
%-learn the Kalman filter
tmpcfg = [];
tmpcfg.init = init;
tmpcfg.C = lf;
tmpcfg.tol = cfg.tol;
tmpcfg.maxiter = cfg.maxiter;
tmpcfg.order = cfg.order;
tmpcfg.roi = roi;
[EM, LL] = ssm_em(tmpcfg, data.trial);
%-----------------%

coeffs = reshape(EM.A, [nroi nroi cfg.order]);
noisecov = EM.Q;
%-------------------------------------%

%-------------------------------------%
%-deal with the output
dataout.dimord = 'chan_chan_lag';
dataout.label = {cfg.roi.label};
dataout.dof = cfg.order;
dataout.fsampleorig = data.fsample;
dataout.coeffs = coeffs;
dataout.noisecov = noisecov;
dataout.cfg = cfg;
dataout.cfg.vol = []; % don't keep vol, it's huge
dataout.cfg.EM = EM;
dataout.cfg.LL = LL;
output = [output sprintf('Log-likelihood: %f after % 3d steps\n', max(LL), numel(find(LL)))]; % TODO: remove zeros/find
%-------------------------------------%
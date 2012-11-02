function init = ssm_em_init(cfg, xhat, C, Cf)
%SSM_EM_INIT generate initial conditions for state-space model
%
% x(t) = A * x(t-1) + w(t), with w ~ N(0, Q)
% y(t) = C * lambda * x(t) + v(t), with v ~ N(0, R)
%
% Use as:
%    init = ssm_em_init(cfg, xhat, C, Cf)
%
% CFG
%   .roi: number of voxels in each ROI
%   .order: model order
%   .Q_range: range for Q (random values will be in the interval range)
%   .R_range: range for R (random values will be in the interval range)
% 
% X_HAT: reconstructed activity of the sources, using MNE or beamforming
%        (a cell with nroi X ntimepoints)
%
% C: forward model (nchan X nroi)
% 
% CF: covariance matrix of the channels (nchan X nchan)
%
% INIT
%   .A: autoregressive model ( (nroi * order) X nroi )
%   .x0: initial state ( (nroi * order) X 1 )
%   .P0: covariance of sources ( (nroi * order) X (nroi * order) X 1 )
%   .Q: covariance of AR model error ( nroi X nroi ),qrs
%   .R: covariance of measuramente error ( nchan X nchan )
% 
% See also SSM_EM, SSM_EM_INIT, SSM_EM_INIT_AR, SSM_EM_KALMAN_FILTER,
% SSM_EM_KALMAN_SMOOTH, SSM_EM_AR

%-------------------------------------%
%-estimate AR model
nroi = numel(cfg.roi);
%-----------------%
%-in cortical patch, take the main principal component
xhat_c = [xhat{:}]; % all cortical patches, but concatenated

phi = zeros(sum(cfg.roi), nroi);
xhat_r = zeros(nroi, size(xhat_c,2)); % xhat for each roi
for i = 1:nroi
  %-index of the sources belonging to one roi
  roi_beg = sum(cfg.roi(1:(i-1))) + 1;
  roi_end = sum(cfg.roi(1:i));
  
  [p1, p2] = princomp(xhat_c(roi_beg:roi_end,:)');
  phi(roi_beg:roi_end,i) = p1(:,1);
  xhat_r(i, :) = p2(:,1)';
  
end

n_trl = numel(xhat);
n_smp = size(xhat{1},2); % make it more flexible, with different trial lengths

xhat_r = mat2cell(xhat_r, nroi, repmat(n_smp, n_trl, 1)); % redefine it in trials
%-----------------%

[AR_hat, Q_hat, P_hat] = ssm_em_init_ar(xhat_r, cfg.order, nroi);
%-------------------------------------%

%-------------------------------------%
%-x(t) = A * x(t-1) + w(t), with w ~ N(0, Q)
%-get initial random values

%-----------------%
%-slightly faster implementation than Cheung
% values will be different because different numbers of calls to rand.
A = zeros(nroi * cfg.order, nroi);
for i = 1:nroi
  A(i:nroi:end, i) = rand(cfg.order, 1) * diff(cfg.A_range) + cfg.A_range(1);
end
%-----------------%

Q = diag(rand(nroi,1) * diff(cfg.Q_range) + cfg.Q_range(1));

x0 = 0.1 * randn(nroi * cfg.order, 1); 
P0 = trace(Q_hat) / nroi * eye(nroi * cfg.order); 
%-------------------------------------%

%-------------------------------------%
%-y(t) = C * lambda * x(t) + v(t), with v ~ N(0, R)
R = Cf - C * phi * P_hat(1:nroi, 1:nroi) * phi' * C'; % covariance which is not explained by source covariance

%-------%
%-no easy interpretation when R is not positive semidefinite, so get rid of
% negative eigenvalues
[V, D] = eig(R);
D(D < 0) = 1e-4; % very small number
R = V * D * V';
R = diag(diag(R));
%-------%
%-------------------------------------%

%-------------------------------------%
%-output
init.phi = phi;
init.A = A;
init.x0 = x0;
init.P0 = P0;
init.Q = Q;
init.R = R;

cfg.init = init;
%-------------------------------------%

function [em, LL] = ssm_em(cfg, y)
%SSM_EM estimate state-space model
%
% x(t) = A * x(t-1) + w(t), with w ~ N(0, Q)
% y(t) = C * lambda * x(t) + v(t), with v ~ N(0, R)
%
% Use as:
%    [em, ll] = ssm_em(cfg, y)
%
% CFG
%   .init
%     .A: autoregressive model ( (nroi * order) X nroi )
%     .x0: initial state ( (nroi * order) X 1 )
%     .P0: covariance of sources ( (nroi * order) X (nroi * order) X 1 )
%     .Q: covariance of AR model error ( nroi X nroi )
%     .R: covariance of measuramente error ( nchan X nchan )
%   .order: model order
%   .maxiter: maximum number of iterations
%   .tol: tolerance
% 
% Y: original data (a cell with nchan X ntimepoints)
%
%
% 
% See also SSM_EM, SSM_EM_INIT, SSM_EM_INIT_AR, SSM_EM_KALMAN_FILTER,
% SSM_EM_KALMAN_SMOOTH, SSM_EM_AR

%p: MVAR model order
%cfg.init.A: initial A matrix with dim = mp x m
%cfg.init.Q: initial Q matrix 
%cfg.init.R: initial R matrix 
%cfg.init.x0: initial x_t_init vector 

%-------------------------------------%
%-input
%-----------------%
%-scalars
order = cfg.order;
order_nsource = size(cfg.init.x0,1);
nroi = order_nsource / order;

ntrl = numel(y);
nchan = size(y{1},1);
nrest = nroi * (order - 1);
%-----------------%

%-----------------%
%-x(t) = A * x(t-1) + w(t), with w ~ N(0, Q)
A  = cfg.init.A'; % TODO: transpose???

B = zeros(order_nsource, nroi);
B(1:nroi,1:nroi) = eye(nroi);

phi = cfg.init.phi;
x0 = cfg.init.x0;
P0 = cfg.init.P0;
Q  = cfg.init.Q;
%-----------------%

%-----------------%
%-y(t) = C * lambda * x(t) + v(t), with v ~ N(0, R)
R  = cfg.init.R;

C_zero = [cfg.C * phi zeros(nchan, nrest)];
%-----------------%

warning off
%-------------------------------------%

%-------------------------------------%
%-iteration
unittol = floor(-1*log10(cfg.tol));

LL = zeros(1, cfg.maxiter);

for i = 1:cfg.maxiter

  %---------------------------%
  %-initialize values for each iteration
  A = [A; eye(nrest) zeros(nrest, nroi)];
  
  x_smooth = cell(1, ntrl);
  X0 = zeros(order_nsource,1);
  E_auto_1 = zeros(order_nsource,order_nsource); % 1:end-1
  E_auto_2 = zeros(order_nsource,order_nsource); % 2:end
  E_cross = zeros(order_nsource,order_nsource);
  E_auto0 = zeros(order_nsource,order_nsource);
  LL(i) = 0;
  %---------------------------%
  
  %---------------------------%
  %-loop over trials
  for e = 1:ntrl
    
    %-----------------%
    %-Kalman filter
    [x_smooth{e}, e_auto, e_cross, LL_e] = ssm_em_kalman_filter(A, B, C_zero, Q, R, x0, P0, y{e});
    %-----------------%
    
    %-----------------%
    %-single-epoch info
    X0 = X0 + x_smooth{e}(:,1);
    E_auto_1 = E_auto_1 + sum(e_auto(:,:,1:end-1),3);
    E_auto_2 = E_auto_2 + sum(e_auto(:,:,2:end),3);
    E_cross = E_cross + sum(e_cross,3);
    E_auto0 = E_auto0 + e_auto(:,:,1);
    LL(i) = LL(i) + LL_e;
    %-----------------%
    
  end
  
  fprintf('% 5g LL: % 15.2f', i, LL(i))
  %---------------------------%
  
  %---------------------------%
  %-check convergence
  if i > 1
    
    LL_d = (LL(i) - LL(i-1)) / abs(LL(i-1));
    fprintf(['      D = % 10.' num2str(unittol) 'f'], LL_d);
    
    if LL_d < cfg.tol ...% TODO: If negative, it should use the A,Q of the previous iteration
        || isinf(LL_d) % if inf, it means that it matches perfectlY (it can happen in simulations)
      fprintf(' CONVERGED \n')
      break
    end
    
  end
  
  fprintf('\n')
  %---------------------------%
  
  %---------------------------%
  %-estimate AR model
  [A, Q, R, phi] = ssm_em_ar(x_smooth, E_auto_1, E_auto_2, E_cross, y, cfg.roi, cfg.C);
  C_zero = [cfg.C * phi zeros(nchan, nrest)];
  
  x0 = X0 / ntrl;
  P0 = E_auto0 / ntrl;
  %---------------------------%
  
end
if i == cfg.maxiter
  fprintf('Reached Max Iterations\n')
end

LL = LL(1:find(LL, 1, 'last'));
%-------------------------------------%

%-------------------------------------%
%-output
em.A = A(1:nroi, 1: (nroi * order)); % remove eye
em.Q = Q;
em.R = R;
em.x0 = x0;
em.x0 = P0;
em.phi = phi;
warning on
%-------------------------------------%
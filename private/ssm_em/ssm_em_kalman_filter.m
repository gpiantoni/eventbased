function  [x_smooth, e_auto, e_cross, LL] = ssm_em_kalman_filter(A, B, C_zero, Q, R, x0, P0, y) 
%INCLUDE FILTERING AND SMOOTHING
% C_zero, because it projects only the activity from the actual sources,
% not from the previous time points

%This function implements the Kalman filter based on the following state-space model
%x[n] = Ax[n-1] + w[n], with w ~ N(0, Q)
%y[n] = Cx[n] + v[n], with v ~ N(0, R) 

%
%Input
%-----
%A_t: State-transition model
%B_t: Control-input model 
%C_t: Observation model
%Q_t: State noise covariance matrix
%R_t: Observation noise covariance matrix
%x_t_init: initial state vector
%P_t_init: covariance of state noise
%y_t_true: observed data
%m_times_k: m*k

%
%Output
%------
%x_t_predicted: predicted state vector
%x_t_updated: updated state vector
%y_t: observation vector
%P_t_predicted: Prediction state error auto covariance matrix
%P_t_updated: Updated state error auto covariance matrix
%P_t_cross_cov: Updated state error cross covariance matrix
%inv_P_t_predicted_err: Inverse of covariance of prediction error sequence (optional)

% x: (modelorder X nsource) X nsmp
% x_hat: (modelorder X nsource) X nsmp
% P: (modelorder X nsource) X (modelorder X nsource) X nsmp

%-----------------%
%-collect values
modelorder_nsource = size(x0,1);
nsmp = size(y,2);
nchan = size(y,1);

LL = 0;
%-----------------%

%-------------------------------------%
%-Kalman Filter
%-----------------%
%-source model
x_p = zeros(modelorder_nsource, nsmp); % priori x_k|k-1
x_k = zeros(modelorder_nsource, nsmp+1); % posteriori x_k|k
P_p = zeros(modelorder_nsource, modelorder_nsource, nsmp); % P_k|k-1
P_k = zeros(modelorder_nsource, modelorder_nsource, nsmp+1); % P_k|k
P_cross = zeros(modelorder_nsource, modelorder_nsource, nsmp); % P_k|k
res = zeros(nchan, nsmp);
%-----------------%

%-----------------%
%-starting points
x_k(:,1) = x0;
P_k(:,:,1) = P0;
%-----------------%

for i = 1:nsmp
  
  %-----------------%
  %-Predict
  x_p(:,i) = A * x_k(:,i);
  P_p(:,:,i) = A * P_k(:,:,i) * A' + B * Q * B';
  %-----------------%
  
  %-----------------%
  %-Update
  res(:,i) = y(:,i) - C_zero * x_p(:,i);
  
  S = C_zero * P_p(:,:,i) * C_zero' + R; % error covariance
  % invS = inv(S); %TODO: make it faster
  K = P_p(:,:,i) * C_zero' / S; % Kalman gain
  
  x_k(:,i+1) = x_p(:,i) + K * res(:,i);
  I_KC = (eye(modelorder_nsource) - K * C_zero);
  P_k(:,:,i+1) = I_KC * P_p(:,:,i);
  P_cross(:,:,i) = I_KC * A * P_k(:,:,i+1);
  %-----------------%
  
  %-----------------%
  %-Log-Likelihood
  LL = LL - 0.5*(log(det(S)) + res(:,i)' / S * res(:,i));
  %-----------------%
  
end
%-------------------------------------%

%-------------------------------------%
%-Kalman Smoother
x_smooth = zeros(modelorder_nsource, nsmp+1); 
e_auto = zeros(modelorder_nsource, modelorder_nsource, nsmp+1);
e_cross = zeros(modelorder_nsource, modelorder_nsource, nsmp);

x_smooth(:,end) = x_k(:,end);
e_auto(:,:,end) = P_k(:,:,end);

for i = (nsmp+1):-1:2
  S = P_k(:,:,i-1) * A' / P_p(:,:,i-1);
  x_smooth(:,i-1) = x_k(:,i-1) + S * (x_smooth(:,i) - x_p(:,i-1));
  e_auto(:,:,i-1) = P_k(:,:,i-1) + S * (e_auto(:,:,i) - P_p(:,:,i-1)) * S';
  e_cross(:,:,i-1) = P_cross(:,:,i-1) + (e_auto(:,:,i) - P_k(:,:,i)) / P_k(:,:,i) * P_cross(:,:,i-1);
end
%-------------------------------------%

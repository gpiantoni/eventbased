function  [x_p, x_k, P_p, P_k, P_cross, LL] = ssm_em_kalman_filter(A, B, C_zero, Q, R, x0, P0, y) 
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


% x: (modelorder X n_source) X n_timepoint
% x_hat: (modelorder X n_source) X n_timepoint
% P: (modelorder X n_source) X (modelorder X n_source) X n_timepoint

modelorder_nsource = size(x0,1);
n_source = size(A,2);
modelorder = modelorder_nsource / n_source;
n_timepoint = size(y,2);
n_chan = size(y,1);

LL = 0;

%-source model
x_p = zeros(modelorder_nsource, n_timepoint); % priori x_k|k-1
x_k = zeros(modelorder_nsource, n_timepoint+1); % posteriori x_k|k
P_p = zeros(modelorder_nsource, modelorder_nsource, n_timepoint); % P_k|k-1
P_k = zeros(modelorder_nsource, modelorder_nsource, n_timepoint+1); % P_k|k
P_cross = zeros(modelorder_nsource, modelorder_nsource, n_timepoint); % P_k|k

C = C_zero(:,1:modelorder_nsource);

res = zeros(n_chan, n_timepoint);

x_k(:,1) = x0;
P_k(:,:,1) = P0;

for i = 1:n_timepoint
  
  % Predict
  x_p(:,i) = A * x_k(:,i);
  P_p(:,:,i) = A * P_k(:,:,i) * A' + B * Q * B';
  
  % Update
  res(:,i) = y(:,i) - C_zero * x_p(:,i); 
  
  S = C_zero * P_p(:,:,i) * C_zero' + R; % error covariance
  % invS = inv(S); %TODO: make it faster
  K = P_p(:,:,i) * C_zero' / S; % Kalman gain
  
  x_k(:,i+1) = x_p(:,i) + K * res(:,i);
  I_KC = (eye(modelorder_nsource) - K * C_zero);
  P_k(:,:,i+1) = I_KC * P_p(:,:,i);
  P_cross(:,:,i) = I_KC * A * P_k(:,:,i+1);
  
  % Log-Likelihood
  LL = LL - 0.5*(log(det(S)) + res(:,i)' / S * res(:,i));
        
end

function [x_smooth, e_auto, e_cross] = ssm_em_kalman_smooth(A, x_p, x_k, P_p, P_k, P_cross)

n_source = size(x_p, 1);
n_point = size(x_p, 2);

x_smooth = zeros(n_source, n_point+1); 
e_auto = zeros(n_source, n_source, n_point+1);
e_cross = zeros(n_source, n_source, n_point);

x_smooth(:,end) = x_k(:,end);
e_auto(:,:,end) = P_k(:,:,end);

for i = (n_point+1):-1:2
  S = P_k(:,:,i-1) * A' / P_p(:,:,i-1);
  x_smooth(:,i-1) = x_k(:,i-1) + S * (x_smooth(:,i) - x_p(:,i-1));
  e_auto(:,:,i-1) = P_k(:,:,i-1) + S * (e_auto(:,:,i) - P_p(:,:,i-1)) * S';
  e_cross(:,:,i-1) = P_cross(:,:,i-1) + (e_auto(:,:,i) - P_k(:,:,i)) / P_k(:,:,i) * P_cross(:,:,i-1);
end

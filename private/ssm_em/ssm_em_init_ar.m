function [A, Q, R] = ssm_em_init_ar(x0, order, nroi)
%SSM_EM_INIT_AR Yule-Walker on the "estimated" roi activity

%-----------------%
%-initialize
ntrl = numel(x0);
cnt = 0;

R_aa = 0;
R_ab = 0;
R_bb = 0;
%-----------------%

%-----------------%
%-estimate covariance for lags
for e = 1:ntrl
  
  modelorder_nroi = order * nroi;
  x_a = zeros(modelorder_nroi, size(x0{e},2) - modelorder_nroi);
  x_b = zeros(modelorder_nroi, size(x0{e},2) - modelorder_nroi);
  
  for i = 1:order
    x_a((i-1) * nroi + 1: i * nroi,:) = x0{e}(1:nroi, 2 * order + 1 - i + 1:end - i + 1);
    x_b((i-1) * nroi + 1: i * nroi,:) = x0{e}(1:nroi, 2 * order + 1 - i    :end - i);
  end
  
  R_aa = R_aa + x_a * x_a';
  R_ab = R_ab + x_a * x_b';
  R_bb = R_bb + x_b * x_b';
  cnt = cnt + size(x_a,2);
end
%-----------------%

%-----------------%
%-output values
A = R_ab(1:nroi,:) / R_bb;
Q = 1 / cnt * (R_aa(1:nroi,1:nroi) - A * R_ab(1:nroi,:)');

R = R_aa / cnt;
%-----------------%
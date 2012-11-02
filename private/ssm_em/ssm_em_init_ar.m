function [w, noise_cov, Rxx_n] = ssm_em_init_ar(x0, order, nroi)


Rxx_n = 0;
Rxx_rhs = 0;
Rxx_n_minus_1 = 0;

n_trl = numel(x0);

aa = 0;

for e = 1:n_trl
  
  x_n_vec = [];
  x_n_minus_1_vec = [];
  
  for i = 1:order
    x_n_vec((i-1) * nroi + 1: i * nroi,:) = x0{e}(1:nroi, 2 * order + 1 - i + 1:end - i + 1);
    x_n_minus_1_vec((i-1) * nroi + 1: i * nroi,:) = x0{e}(1:nroi, 2 * order + 1 - i:end - i);
  end
  
  Rxx_n = Rxx_n + x_n_vec * x_n_vec';
  Rxx_rhs = Rxx_rhs + x_n_vec * x_n_minus_1_vec';
  Rxx_n_minus_1 = Rxx_n_minus_1 + x_n_minus_1_vec * x_n_minus_1_vec';
  aa = aa + size(x_n_vec,2);
end

w = Rxx_rhs(1:nroi,:) / Rxx_n_minus_1;
noise_cov = 1 / aa * (Rxx_n(1:nroi,1:nroi) - (Rxx_rhs(1:nroi,:)/Rxx_n_minus_1) * Rxx_rhs(1:nroi,:)');

Rxx_n = Rxx_n / aa;

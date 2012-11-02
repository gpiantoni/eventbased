function [A_out, Q_out, R_out, phi] = ssm_em_ar(x_smooth, E_auto_1, E_auto_2, E_cross, y, roi, C_h)
% one roi contains multiple sources

nsource = sum(roi);
nroi = numel(roi);

modelorder_nroi = size(x_smooth{1},1);
ntrl = numel(x_smooth);
nchan = size(y{1},1);
modelorder = modelorder_nroi / nroi;

A = zeros(modelorder_nroi, modelorder_nroi);
B = zeros(modelorder_nroi, modelorder_nroi);
C = zeros(modelorder_nroi, modelorder_nroi);
R = zeros(nchan, nchan);
yx = zeros(nchan, nroi);
yy = zeros(nchan, nchan);

cnt = 0;
for e = 1:ntrl
  
  %Forming Rxx[k] = E{x[n]x[n-k]'}.
  for i = modelorder+1 : size(x_smooth{e},2)-1
    A = A + x_smooth{e}(:,i)   * x_smooth{e}(:,i)';
    B = B + x_smooth{e}(:,i+1) * x_smooth{e}(:,i)';
    C = C + x_smooth{e}(:,i+1) * x_smooth{e}(:,i+1)';
    
    if any(roi ~= 1) % at least one roi is not a dipole
      yx = yx + y{e}(:,i) * x_smooth{e}(1:nroi,i+1)';
      yy = yy + y{e}(:,i) * y{e}(:,i)';
      
    else
      
      R = R + (y{e}(:,i) - C_h * x_smooth{e}(1:nroi, i+1)) * (y{e}(:,i) - C_h * x_smooth{e}(1:nroi, i+1))';
    end
    
    cnt = cnt + 1;
    
  end
end

A = A + E_auto_1;
B = B + E_cross;
C = C + E_auto_2;

A_out = B(1:nroi,:) / A; % TODO: transpose?
Q_out = (C(1:nroi,1:nroi) - B(1:nroi,:) / A' * B(1:nroi,:)') / cnt;

if any(roi ~= 1) % TODO: use same implementation for patch and dipole
  
  phi = zeros(nsource, nroi);
  invS = inv( (yy - yx / C(1:nroi, 1:nroi) * yx') / cnt);
  C_interest = C(1:nroi, 1:nroi); % why this choice?
  
  A = [];
  b = [];

  for i = 1:nroi
    
    %-index of the sources belonging to one roi
    roi_beg = sum(roi(1:(i-1))) + 1;
    roi_end = sum(roi(1:i));

    psi_i = C_h(:,roi_beg:roi_end);
    temp = inv(psi_i' * invS * psi_i) * psi_i' * invS;
        
    A_row_wise = [];
    for j = 1:nroi
      
      if j == i
        B = eye(roi(i));
        C1 = 1;
        
      else
        %-index of the sources belonging to one roi
        roi_beg = sum(roi(1:(j-1))) + 1;
        roi_end = sum(roi(1:j));
        
        B = temp * C_h(:,roi_beg:roi_end);
        C1 = C_interest(j,i) / C_interest(i,i);
        
      end
      
      A_row_wise = [A_row_wise kron(C1',B)];

    end
    
    A = [A; A_row_wise];
    b_mtx = temp * yx(:,i) / C_interest(i,i);
    b = [b; b_mtx(:)];
  end
  phi_new_vec = A\b;

  for i = 1:nroi
    %-index of the sources belonging to one roi
    roi_beg = sum(roi(1:(i-1))) + 1;
    roi_end = sum(roi(1:i));
    
    source_in_roi = phi_new_vec(roi_beg:roi_end,:);
    
    %Normalization on each column
    phi(roi_beg:roi_end,i) =  source_in_roi / norm(source_in_roi);
    
  end
  
  C_h = C_h * phi;
  R_out = (yy - C_h * yx' - yx * C_h' + C_h * C(1:nroi,1:nroi) * C_h') / cnt;
  
else
  phi = eye(nroi);
  R_out = (R + C_h * E_auto_2(1:nroi,1:nroi) * C_h') / cnt;
  
end
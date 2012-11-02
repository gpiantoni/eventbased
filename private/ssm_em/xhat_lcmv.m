function x = xhat_lcmv(x, lf, cov)
%XHAT_LCMV: compute xhat using a straightforward beamforming

%-----------------%
%-calculate weights
invRz = pinv(cov); % pinv instead of inv for rank-deficient cases
weight = inv(lf' * invRz * lf) * lf' * invRz;
%-----------------%

%-----------------%
%-Apply weight for each epoch
for i = 1:numel(x)
    x{i} = weight * x{i};
end
%-----------------%
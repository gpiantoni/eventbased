function [w,noise_cov,Rxx_n,Rxx_rhs,Rxx_n_minus_1]=multichannel_YW_v3(z,p,m)
%function [w,noise_cov,x_n_vec,x_n_minus_1_vec]=multichannel_YW_v3(z,p,m)
%Implementation matches with the M-step of the EM algorithm

%Input
%-----
%z: input data (mxTxJ)
%p: model order
%m: dim of multivariate

%Output
%------
%w: AR model parameters

%Form Rij[0],Rij[1], ..., Rij[-(p-1)]
%Note: Rij[k]' = Rij[-k]

Rxx_n = 0;
Rxx_rhs = 0;
Rxx_n_minus_1 = 0;

if (iscell(z))
    cell_struct = 1;
else
    cell_struct = 0;
end

if (cell_struct)
    N_epch = size(z,2);
else
     N_epch = size(z,3);
end

aa = 0;

for kk = 1:N_epch
    %Forming Rxx[k] = E{x[n]x[n-k]'}. 
    x_n_vec = [];
    x_n_minus_1_vec = [];
    for ii = 1:p
        if (cell_struct)
            temp = z{kk};
            x_n_vec((ii-1)*m+1:ii*m,:) = temp(1:m,2*p+1-ii+1:end-ii+1);
            x_n_minus_1_vec((ii-1)*m+1:ii*m,:) = temp(1:m,2*p+1-ii:end-ii);
        else
            % Transient samples z[1] ... z[p] are not used.
            x_n_vec((ii-1)*m+1:ii*m,:) = z(1:m,2*p+1-ii+1:end-ii+1,kk);
            x_n_minus_1_vec((ii-1)*m+1:ii*m,:) = z(1:m,2*p+1-ii:end-ii,kk);
        end
    end
    
    Rxx_n = Rxx_n + x_n_vec*x_n_vec';
    Rxx_rhs = Rxx_rhs + x_n_vec*x_n_minus_1_vec';
    Rxx_n_minus_1 = Rxx_n_minus_1 + x_n_minus_1_vec*x_n_minus_1_vec';
    aa = aa + size(x_n_vec,2);
end

w = Rxx_rhs(1:m,:)/Rxx_n_minus_1;
noise_cov = 1/aa*(Rxx_n(1:m,1:m)-(Rxx_rhs(1:m,:)/Rxx_n_minus_1)*Rxx_rhs(1:m,:)');

Rxx_rhs = Rxx_rhs/aa;
Rxx_n = Rxx_n/aa;
Rxx_n_minus_1 = Rxx_n_minus_1/aa;

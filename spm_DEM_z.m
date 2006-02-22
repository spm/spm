function [z,C,K] = spm_DEM_z(M,N)
% creates hierarchical innovations for generating data
% FORMAT [z] = spm_DEM_z(M,N)
% M    - model structure
% N    - length of data sequence
%
% z{i} - innovations for level i (N.B. z{end} corresponds to causes)
% C    - component covariance
% K    - autocovariance
%
%--------------------------------------------------------------------------
% If the precision is zero, unit variance innovations are assumed
%
%                    vec(z{i}) ~ N(0,kron(K,C))
%__________________________________________________________________________
% Copyright (C) 2005 Wellcome Department of Imaging Neuroscience

% Karl Friston
% $Id$
 
% set model structure
%--------------------------------------------------------------------------
M     = spm_DEM_M_set(M);
 
% temporal convolution matrix
%--------------------------------------------------------------------------
s     = M(1).E.s + eps;
dt    = M(1).E.dt;
t     = ([1:N] - 1)*dt;
K     = toeplitz(exp(-t.^2/(2*s^2)));
 
 
% create innovations z{i} [assume unit variance if precision is 0]
%--------------------------------------------------------------------------
for i = 1:length(M)
 
    P     = M(i).V;
    for j = 1:length(M(i).Q)
        P = P + M(i).Q{j}*exp(M(i).h(j));
    end
    if ~normest(P,1)
        C{i} = speye(M(i).l,M(i).l);
    else
        C{i} = inv(P);
    end
    z{i}  = spm_sqrtm(C{i})*randn(M(i).l,N)*K;
    
end
 
% create implicit temporal covariance matrix
%--------------------------------------------------------------------------
if nargout > 2
    K     = K*K';
end
 



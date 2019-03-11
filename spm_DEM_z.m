function [z,w] = spm_DEM_z(M,N)
% creates hierarchical innovations for generating data
% FORMAT [z w] = spm_DEM_z(M,N)
% M    - model structure
% N    - length of data sequence
%
% z{i} - innovations for level i (N.B. z{end} corresponds to causes)
% w{i} - innovations for level i (state noise)
%
% If there is no fixed or hyper parameterized precision, then unit noise is
% created. It is assumed that this will be later modulated by state
% dependent terms, specified by M.ph and M.pg in spm_DEM_int
%__________________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% Karl Friston
% $Id: spm_DEM_z.m 7540 2019-03-11 10:44:51Z karl $

% temporal convolution matrix (with unit variance)
%--------------------------------------------------------------------------
s  = M(1).E.s + exp(-16);
dt = M(1).E.dt;
t  = ((1:N) - 1)*dt;
K  = toeplitz(exp(-t.^2/(2*s^2)));
K  = K*diag(1./sqrt(diag(K*K')));

% create innovations z{i} and w{i}
%--------------------------------------------------------------------------
for i = 1:length(M)
    
    % precision of causes
    %======================================================================
    P     = M(i).V;
    
    % plus prior expectations
    %----------------------------------------------------------------------
    try
        for j = 1:length(M(i).Q)
            P = P + M(i).Q{j}*exp(M(i).hE(j));
        end
    end
    
    % create causes: assume i.i.d. if precision is zero
    %----------------------------------------------------------------------
    if norm(P,1) == 0;
        z{i}  = randn(M(i).l,N)*K;
    elseif norm(P,1) >= exp(16)
        z{i}  = sparse(M(i).l,N);
    else
        z{i}  = spm_sqrtm(inv(P))*randn(M(i).l,N)*K;
    end
    
    % precision of states
    %======================================================================
    P     = M(i).W;
    
    % plus prior expectations
    %----------------------------------------------------------------------
    try
        for j = 1:length(M(i).R)
            P = P + M(i).R{j}*exp(M(i).gE(j));
        end
    end
    
    % create states: assume i.i.d. if precision (P) is zero
    %----------------------------------------------------------------------
    if ~isempty(P)
        if norm(P,1) == 0;
            w{i} = randn(M(i).n,N)*K*dt; 
        elseif norm(P,1) >= exp(16)
            w{i} = sparse(M(i).n,N);
        else
            w{i} = spm_sqrtm(inv(P))*randn(M(i).n,N)*K*dt;
        end
    else
        w{i} = sparse(0,0);
    end
    
end



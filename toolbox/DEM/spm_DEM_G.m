function [z,w] = spm_DEM_G(M)
% Create generalised innovations for generative process
% FORMAT [z,w] = spm_DEM_G(M)
% M    - model structure
%
% z{i} - innovations for level i (N.B. z{end} corresponds to causes)
% w{i} - innovations for level i (state noise)
%
% If there is no fixed or hyper parameterized precision, then unit noise is
% created. It is assumed that this will be later modulated by state
% dependent terms, specified by M.ph and M.pg in spm_DEM_int
%__________________________________________________________________________

% Karl Friston
% Copyright (C) 2005-2022 Wellcome Centre for Human Neuroimaging


% create innovations z{i} and w{i}
%--------------------------------------------------------------------------
for i = 1:length(M)
    
    % precision of causes
    %======================================================================
    P     = M(i).V;
    
    % plus prior expectations
    %----------------------------------------------------------------------
    if isfield(M,'Q')
        for j = 1:length(M(i).Q)
            P = P + M(i).Q{j}*exp(M(i).hE(j));
        end
    end
    
    % create causes: assume i.i.d. if precision is zero
    %----------------------------------------------------------------------
    if norm(P,1) == 0
        z{i}  = randn(M(i).l,1);
    else
        z{i}  = spm_sqrtm(inv(P))*randn(M(i).l,1);
    end
    
    % precision of states
    %======================================================================
    P     = M(i).W;
    
    % plus prior expectations
    %----------------------------------------------------------------------
    if isfield(M,'R')
        for j = 1:length(M(i).R)
            P = P + M(i).R{j}*exp(M(i).gE(j));
        end
    end
    
    % create states: assume i.i.d. if precision (P) is zero
    %----------------------------------------------------------------------
    if ~isempty(P)
        if norm(P,1) == 0
            w{i} = randn(M(i).n,1); 
        else
            w{i} = spm_sqrtm(inv(P))*randn(M(i).n,1);
        end
    else
        w{i} = sparse(0,0);
    end
    
end

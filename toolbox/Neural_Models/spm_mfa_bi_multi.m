function [M0,M1,L,M2] = spm_mfa_bi_multi(S,C)
% Bilinear form for multiple Gibb's ensembles
% FORMAT [M0,M1,L,M2] = spm_mfa_bi_multi(S,C)
%--------------------------------------------------------------------------
% S(s)   - MFA system specification structure for s ensembles
% Required fields:
%     S(i).M:  [1x1 struct]    - dynamic model structure
%     S(i).J0: [n x n double]  - Jacobian
%     S(i).J1: {1 x M.m cell}  - dJ0/du
%     S(i).L : [l x n double]  - d<y>/dp
%     S(i).u:  [n x m double]  - probability modes
%     S(i).v:  [m x n double]  - v*u = 1
%     S(i).X:  [n x d double]  - evaluation points of state [d]-space
%     S(i).x:  {1 x d cell}    - range of state [d]-space
%     S(i).p0: [n x 1 sparse]  - expansion point
%
% C{s x s cell} - coupling cell = dP/d<y> (change in parameters of S(i).M with
%                                   mean outputs of S(j).M - [p x l double])
%
% M0  [ns + 1 x ns + 1double]  - 1st order Bilinear matrix dq/dt;
% M1  {M.m x s}                - 2nd order Bilinear matrix dM0/du
% M2  {s x s}                  - 2nd order Bilinear matrix dM0/dC
% L   [1s x ns + 1 double]     - output matrix          <y> = L*q;
%
% Transformed probability states:  q = [1; v*(p(X) - p0)];
%__________________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging
 
% Karl Friston
% $Id: spm_mfa_bi_multi.m 4936 2012-09-18 19:47:55Z karl $
 
 
% indices
%--------------------------------------------------------------------------
s     = length(S);                      % number of ensembles
m     = length(S(1).J1);                % number of inputs
I     = {(1:size(S(1).v,1)) + 1};       % ensemble indices
for i = 2:s
    I{i} = (1:size(S(i).v,1)) + I{i - 1}(end);
end
q     = I{s}(end);
M0    = sparse(q,q);
M1    = {};
L     = sparse(0,q);
 
% compute M0 = df/dx. M1 = d(df/dx)/du amd L = dy/dx
%--------------------------------------------------------------------------
for i = 1:s
 
    % Ist order bilinear operator M0
    %======================================================================
    M0(I{i},1)       = S(i).v*S(i).J0*S(i).p0;
    M0(I{i},I{i})    = S(i).v*S(i).J0*S(i).u;
 
    % 2nd order bilinear operators M1{i} - dM0/du
    %======================================================================
    for j = 1:m
        M            = sparse(q,q);
        M(I{i},1)    = S(i).v*S(i).J1{j}*S(i).p0;
        M(I{i},I{i}) = S(i).v*S(i).J1{j}*S(i).u;
        M1{m,i}      = M;
    end
 
    % output matrix <y>
    %======================================================================
    j         = (1:size(S(i).L,1)) + size(L,1);
    L(j,1)    = S(i).L*S(i).p0;
    L(j,I{i}) = S(i).L*S(i).u;
 
 
end
 
% return unless coupling matrices are required
%--------------------------------------------------------------------------
if nargout < 4
    return
end
 
% 2nd order bilinear operators M2{i} - dM0/dC
%==========================================================================
M2    = cell(s,s);
dP    = 1e-1;
for i = 1:s
 
    % target ensemble df/dP = dJ/dP*p0
    %----------------------------------------------------------------------
    M     = S(i).M;
    p     = length(M.pE);
    for k = 1:p
        M.pE      = S(i).M.pE;
        M.pE(k)   = M.pE(k) + dP;
        dJdP      = (spm_mfa(M,S(i).x) - S(i).J0)/dP;
        dfdP(:,k) = S(i).v*dJdP*S(i).p0;
    end
 
    for j = 1:s
 
        % source ensemble dP/dx = dP/d<y>*d<y>/dx = C*L
        %------------------------------------------------------------------
        dfdx      = sparse(q,q);
        if ~isempty(C{i,j})
            dPdx            = C{i,j}*S(j).L*S(j).u;
            dfdx(I{i},I{j}) = dfdP*dPdx;
        end
 
        % M2 = df/dx = df/dP*dP/dx
        %------------------------------------------------------------------
        M2{i,j} = dfdx;
 
    end
end

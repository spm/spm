function [M0,M1,L1] = spm_mfa_bi(M,P)
% Bilinear form as a function of coupling parameters
% FORMAT [M0,M1,L1] = spm_mfa_bi(M,P)
%--------------------------------------------------------------------------
% M  - MFA network specification structure
% Required fields:
%   N.bi  - 'spm_mfa_bi';
%   N.M0  - 1st order bilinear operator;
%   N.M1  - dM1/dPu;
%   N.M2  - dM0/dPc;
%   N.L   - d<y>/dq;
% P     - input and coupling parameters P = [Pu, Pc]
%
% M0 [n x n double] - 1st order Lie matrix 
% M1 {1 x m cell}   - 2nd order Lie matrix 
% L1 {l x n cell}   - output matrix (1st order)    <y> = L1*q 
%
% Transformed states     q = [1; v*(p(X) - p0)];
%
%                dq/dt = M0*q + u(1)*M1{1}*q + ...;
%
%__________________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging
 
% Karl Friston
% $Id: spm_mfa_bi.m 5615 2013-08-15 14:37:24Z spm $

% reshape coupling parameters
%--------------------------------------------------------------------------
[m s] = size(M.M1);             % number of inputs
[d s] = size(M.M2);             % number of regions

p     = m*s;
Pu    = reshape(P(1:p),m,s);
Pc    = reshape(P((p + 1):end),s,s);

% compute new M0
%--------------------------------------------------------------------------
M0    = M.M0;
for i = 1:s
    for j = 1:s
        M0 = M0 + M.M2{i,j}*Pc(i,j);
    end
end

% compute new M1
%--------------------------------------------------------------------------
for i = 1:m
    M1{i} = sparse(1,1);
    for j = 1:s
        M1{i} = M1{i} + M.M1{i,j}*Pu(i,j);
    end
end

% set output L
%--------------------------------------------------------------------------
L1   = M.L;

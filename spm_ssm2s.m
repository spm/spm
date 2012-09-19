function [S] = spm_ssm2s(P,M,U)
% Converts state-space (M) representation to eigenspectrum
% FORMAT [S] = spm_ssm2s(P,M,U)
%
% P    - model parameters
% M    - model (with flow M.f and expansion point M.x and M.u)
% U    - exogenous inputs
%
% S    - (SPM) eigenspectrum
% Hz   - vector of frequencies (Hz)
%
% csd  - cross spectral density
%__________________________________________________________________________
% Copyright (C) 2012 Wellcome Trust Centre for Neuroimaging
 
% Karl Friston
% $Id: spm_ssm2s.m 4937 2012-09-19 10:40:26Z guillaume $
 

% Jacobian and eigenspecturm
%--------------------------------------------------------------------------
dfdx   = spm_diff(M.f,M.x,M.u,P,M,1);
[u,s]  = eig(full(dfdx));
s      = diag(s);

% principal eigenmodes (highest frequency)
%--------------------------------------------------------------------------
[q,i]  = sort(imag(s),'descend');
i      = i(1:3);
S      = s(i);
s(i)   = [];


% principal eigenmodes (most unstable)
%--------------------------------------------------------------------------
[u,i]  = sort(real(s),'descend');
S      = [S; u(1)];

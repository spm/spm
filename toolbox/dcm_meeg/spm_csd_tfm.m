function [y,w,t,x] = spm_csd_tfm(P,M,U)
% Time frequency response of a neural mass model
% FORMAT [y,w,t,x] = spm_csd_tfm(P,M,U)
%
% P - parameters
% M - neural mass model structure
% U - time-dependent input
%
% y - {y(t,w,nc,nc}} - cross-spectral density for nc channels {trials}
%                    - for w frequencies over time t in M.Hz
% w - frequencies
% t – peristimulus time
% x - expectation or hidden (neuronal) states
%
%__________________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% Karl Friston
% $Id: spm_csd_tfm.m 4714 2012-04-10 13:30:44Z karl $


% check input
%--------------------------------------------------------------------------
if nargin < 3
    U.dt = 1;
    U.u  = sparse(1,M.m);
end

% peristimulus time
%--------------------------------------------------------------------------
t    = (1:size(U.u,1))*U.dt;

% solve for expected hidden states (by removing observer function)
%==========================================================================
S    = M;
S.g  = {};
x    = spm_int_L(P,S,U);

% cycle over time – expanding around expected states and input
%--------------------------------------------------------------------------
for i = 1:length(t)
    
    % update hidden neuronal states and time-dependent parameters
    %----------------------------------------------------------------------
    M.x   = spm_unvec(x(i,:),M.x);
    [u,Q] = spm_fx_tfm_P(U.u(i,:),P);
    
    % compute complex cross spectral density
    %----------------------------------------------------------------------
    [g,w] = spm_csd_mtf(Q,M,U);
    
    % and place in response
    %----------------------------------------------------------------------
    for c = 1:length(g)
        y{c}(i,:,:,:) = g{c};
    end
    
end

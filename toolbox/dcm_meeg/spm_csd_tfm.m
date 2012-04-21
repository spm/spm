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
% $Id: spm_csd_tfm.m 4721 2012-04-21 08:53:29Z karl $


% check input - one trial (no between-tria effects)
%--------------------------------------------------------------------------
if nargin < 3
    U.dt = 0.004;
    U.u  = sparse(1,M.m);
    U.X  = sparse(1,0);
end

% peristimulus time
%--------------------------------------------------------------------------
t   = (1:size(U.u,1))*U.dt;

% within-trial (exogenous) inputs
%==========================================================================
if ~isfield(U,'u')
    
    % check input u = f(t,P,M)
    %----------------------------------------------------------------------
    try, fu  = M.fu; catch,  fu  = 'spm_erp_u'; end
    try, ns  = M.ns; catch,  ns  = 128;         end
    try, dt  = U.dt; catch,  dt  = 0.004;       end
    
    % peri-stimulus time inputs
    %----------------------------------------------------------------------
    U.u = feval(fu,(1:ns)*dt,P,M);
    
end

% between-trial (experimental) inputs
%==========================================================================
try
    X = U.X;
    if ~size(X,1)
        X = sparse(1,0);
    end
catch
    X = sparse(1,0);
end

% disperson (FWHM) of time-frequency responses
%==========================================================================
st    = 4*1e-3/U.dt;
sw    = 4;

% cycle over trials
%--------------------------------------------------------------------------
for  c = 1:size(X,1)
    
    % baseline parameters
    %----------------------------------------------------------------------
    Q  = P;
    
    % trial-specific effects
    %----------------------------------------------------------------------
    for i = 1:size(X,2)
        
        % extrinsic (forward and backwards) connections
        %------------------------------------------------------------------
        for j = 1:length(A)
            Q.A{j} = Q.A{j} + X(c,i)*P.B{i};
        end
        
        % intrinsic connections
        %------------------------------------------------------------------
        try
            Q.H(:,1) = Q.H(:,1) + X(c,i)*diag(P.B{i});
        catch
            Q.G(:,1) = Q.G(:,1) + X(c,i)*diag(P.B{i});
        end
        
    end
    
    % solve for expected hidden states (by removing observer function)
    %======================================================================
    S     = rmfield(M,'g');
    x     = spm_int_J(Q,S,U);
    
    % reconfigure model for endogenous inputs
    %----------------------------------------------------------------------
    X     = Q.X;
    Y     = Q.Y;
    Q     = rmfield(Q,{'X','Y'});
    M.f   = M.h;
    
    
    
    % cycle over time – expanding around expected states and input
    %======================================================================
    for i = 1:length(t)
        
        % update states-dependent parameters
        %------------------------------------------------------------------
        dQ    = X*U.u(i,:)' + Y*x(i,:)';
        R     = spm_unvec(spm_vec(Q) + dQ,Q);
        
        
        % update expansion point
        %------------------------------------------------------------------
        M.x   = spm_unvec(x(i,:),M.x);
        M.u   = 1/128 + exp(R.C)*U.u(i,:)';
        
        
        % compute complex cross spectral density
        %------------------------------------------------------------------
        [g,w] = spm_csd_mtf(R,M);
        
        
        % and place in response
        %------------------------------------------------------------------
        y{c}(i,:,:,:) = g{:};

    end
    
    % model dispersion in time and frequency
    %------------------------------------------------------------------
    for j = size(y{c},3);
        for k = size(y{c},4);
            y{c}(:,:,j,k) = spm_conv(y{c}(:,:,j,k),st,sw);
        end
    end
    
end







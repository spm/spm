function [Y,w,t,x,G,S,E] = spm_csd_int(P,M,U)
% Time frequency response of a neural mass model
% FORMAT [Y,w,t,x,G,S,E] = spm_csd_int(P,M,U)
%
% P - parameters
% M - neural mass model structure
% U - time-dependent input
%
% Y - {Y(t,w,nc,nc}} - cross-spectral density for nc channels {trials}
%                    - for w frequencies over time t in M.Hz
% w - frequencies
% t – peristimulus time
% x - expectation of hidden (neuronal) states
% G - {G(t,w,nc,nc}} - cross spectrum density before dispersion
% S - {S(t,w,nc,nu}} - transfer functions
% E - {E(t,nc}}      - event-related average (sensor space)
%__________________________________________________________________________
%
% This integration routine evaluates the responses of a neural mass model
% to exogenous input – in terms of neuronal states. These are then used as
% expansion point to generate complex cross spectral responses due to
% random neuronal fluctuations. The ensuing spectral (induced) response is
% then convolved (in time) with a window that corresponds to the window of
% a standard wavelet transform. In other words, this routine  generates
% predictions of data features based upon a wavelet transforms
% characterisation of induced responses.
%__________________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% Karl Friston
% $Id: spm_csd_int.m 4814 2012-07-30 19:56:05Z karl $


% check input - one trial (no between-tria effects)
%--------------------------------------------------------------------------
if nargin < 3
    U.dt = 0.004;
    U.u  = sparse(1,M.m);
    U.X  = sparse(1,0);
end


% check function format
%--------------------------------------------------------------------------
f   = fcnchk(M.f);

% check input function  u = f(t,P,M)
%--------------------------------------------------------------------------
try, fu  = M.fu;    catch, fu  = 'spm_erp_u'; end
try, ns  = M.ns;    catch, ns  = 128;         end
try, Rft = M.Rft;   catch, Rft = 4;           end
try, dt  = U.dt;    catch, dt  = 0.004;       end


% within-trial (exogenous) inputs
%==========================================================================
if ~isfield(U,'u') 
    u = feval(fu,(1:ns)*dt,P,M)';
else
    u = U.u';
end

% peristimulus time
%--------------------------------------------------------------------------
t   = (1:size(u,2))*U.dt;

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
nu    = length(P.A{1});

% accelerate bilinear reduction by assuming inputs enter linearly
%--------------------------------------------------------------------------
[dfdxu{1:nu}] = deal(sparse(M.n,M.n));


% cycle over trials
%--------------------------------------------------------------------------
for c = 1:size(X,1)
    
    % baseline parameters
    %----------------------------------------------------------------------
    Q = P;
    
    % trial-specific effects
    %----------------------------------------------------------------------
    for i = 1:size(X,2)
        
        % extrinsic (forward and backwards) connections
        %------------------------------------------------------------------
        for j = 1:length(P.A)
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
    
    % initialise steady state
    %----------------------------------------------------------------------
    try
        x = spm_vec(spm_dcm_neural_x(Q,M));
    catch
        x = M.x;
    end
    
    % cycle over time – expanding around expected states and input
    %======================================================================
    if isfield(Q,'X')
        Q = rmfield(Q,'X');
        Q = rmfield(Q,'Y');
    end
    for i = 1:length(t)
        
        % hidden states
        %------------------------------------------------------------------
        if i == 1
            x(:,i) = spm_vec(M.x);
        else
            x(:,i) = x(:,i - 1);
        end
        
        
        % state-dependent parmeters
        %==================================================================
        
        % update state-dependent parameters (first order)
        %------------------------------------------------------------------
        dQ  = 0;
        if isfield(P,'X')
            dQ  = P.X*u(:,i) + P.Y*x(:,i);
        end
        
        
        % update state-dependent parameters (second order)
        %------------------------------------------------------------------
        % dQ = dQ + spm_csd_cch(Q,M); dQ
        R   = spm_unvec(spm_vec(Q) + dQ,Q);
        
        
        % flow dx(t)/dt and Jacobian df/dx
        %------------------------------------------------------------------
        if nargout(f) == 3
            [fx dfdx D] = f(x(:,i),u(:,i),R,M);
            
        elseif nargout(f) == 2
            [fx dfdx]   = f(x(:,i),u(:,i),R,M);
            D           = 1;
            
        else
            fx          = f(x(:,i),u(:,i),P,M);
            dfdx        = spm_cat(spm_diff(f,x(:,i),u(:,i),R,M,1));
            D           = 1;
        end
        

        % update expansion point (hidden states and endogenous input)
        %------------------------------------------------------------------
        M.x     = spm_unvec(x(:,i),M.x);
        M.u     = sparse(nu,1) + exp(R.C)*u(:,i);
        
        M.dfdxu = dfdxu;
        M.dfdx  = dfdx;
        M.dfdu  = spm_diff(f,M.x,M.u,R,M,2);
        M.f0    = fx; 
        M.D     = D;      
        
        % compute complex cross spectral density
        %==================================================================
        [g,w,s] = spm_csd_mtf(R,M);
        
        M       = rmfield(M,'u');
        M.x     = spm_unvec(x(:,1),M.x);
        
        
        % and place in response
        %------------------------------------------------------------------
        G{c}(i,:,:,:) = g{1};
        S{c}(i,:,:,:) = s{1};
        
        % update states: dx = (expm(dt*J) - I)*inv(J)*fx
        %------------------------------------------------------------------
        x(:,i)  = x(:,i) + spm_dx(D*dfdx,D*fx,dt);
        
        % and response
        %------------------------------------------------------------------
        e(:,i)  = feval(M.g,x(:,i),u(:,i),R,M);
        
    end
    
    % model dispersion associated with wavelet transforms
    %------------------------------------------------------------------
    Y{c}  = spm_morlet_conv(G{c},w*dt,Rft);
    E{c}  = e';
    
end







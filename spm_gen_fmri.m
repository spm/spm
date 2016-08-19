function [y,lfp,csd,w] = spm_gen_fmri(P,M,U)
% Generates a prediction of (multimodal) source activity
% FORMAT [y,lfp,csd,w] = spm_gen_fmri(P,M,U)
%
% P - parameters
% M - neural-mass model structure
% U - trial-effects
%   U.u  - inputs
%   U.dt - (micro) time bins for within-trial effects
%
% y    - BOLD predictions (for every TR)
% lfp  - voltages and conductances (for every micotime bin)
% csd  - spectral density (for every TR)
% w    - frequencies
%
% This integration scheme returns a prediction of neuronal responses to
% experimental inputs, in terms of BOLD responses and, if requested, local
% field potentials and spectral density responses in each region or source.
%
% The scheme uses a canonical microcircuit nneuron mass model of each
% region to evaluate the new fixed point of neuronal activity every time
% the input changes. This is evaluated in microtime (usually a 16th of the
% TR). These neuronal states are then used to compute the pre-synaptic
% activity of (extrinsic and intrinsic) afferents to each subpopulation to
% furnish a neurovascular signal. The ensuing haemodynamic response is then
% estimated by integrating a haemodynamic model. Neurovascular coupling
% depends upon the mixtures of pre-synaptic activity driving haemodynamic
% model. The associated weights are free parameters.
%__________________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% Karl Friston
% $Id: spm_gen_fmri.m 6857 2016-08-19 15:17:06Z karl $


% persistent variables to speed integration schemes
%==========================================================================
persistent pP pf qf rf

try
    INT = any(spm_vec(pP) - spm_vec(P.N));
catch
    INT = 1;
end
pP      = P.N;

% integrate neural system - when imput changes
%==========================================================================
if INT || nargout > 1
    
    TOL = exp(-6);                     % tolerance for fixed point
    x   = spm_vec(M.x.x);
    u   = U.u(1,:) - 1;
    D   = 1;
    
    % neuronal model
    %----------------------------------------------------------------------
    N.f = spm_funcheck(M.fn);
    N.x = M.x.x;
    N.m = numel(u);
    f   = @M.f;
    
    % lead field model
    %----------------------------------------------------------------------
    N.g            = @spm_gx_erp;
    N.l            = M.l;
    N.dipfit.model = 'TFM';
    N.dipfit.type  = 'LFP';
    N.dipfit.Ns    = M.l;
    N.dipfit.Nc    = M.l;
    

    
    for t = 1:size(U.u,1)
        
        % recompute connectivity parameters if input has changed
        %------------------------------------------------------------------
        if any(u - U.u(t,:));
            
            % input
            %--------------------------------------------------------------
            u  = U.u(t,:);
            
            % recompute input dependent parameters
            %--------------------------------------------------------------
            Q  = P.N;
            for i = 1:N.m
                
                % extrinsic connections
                %----------------------------------------------------------
                B{1}   = Q.B{1}(:,:,i);
                B{2}   = Q.B{2}(:,:,i);
                G{1}   = diag(B{1});
                G{2}   = diag(B{2});
                Q.A{1} = Q.A{1} + u(i)*(B{1} - diag(G{1}));
                Q.A{3} = Q.A{3} + u(i)*(B{2} - diag(G{2}));
                
                % intrinsic connections
                %----------------------------------------------------------
                Q.G(:,2) = Q.G(:,2) + u(i)*G{1};
                Q.G(:,4) = Q.G(:,4) + u(i)*G{2};
                
            end
            
            % Jacobian df/dx, checking for neuronal delay operator
            %--------------------------------------------------------------
            if nargout(N.f) >= 3
                [fx,dfdx,D] = N.f(x,u,Q,N);
            elseif nargout(f) == 2
                [fx,dfdx]   = N.f(x,u,Q,N);
            else
                dfdx = spm_cat(spm_diff(N.f,x,u,Q,N,1));
            end
            dfdx  = D*dfdx;
                       
            % condition unstable eigenmodes
            %--------------------------------------------------------------
            [v,s] = eig(full(dfdx),'nobalance');
            s     = diag(s);
            s     = 1j*imag(s) + real(s) - exp(real(s));
            
            % dx  = (expm(dt*df/dx) - I)*inv(df/dx)*f(x,u)
            %--------------------------------------------------------------
            s     = (exp(U.dt*s) - 1)./s;
            J     = real(v*diag(s)*pinv(v));
            
            % reset update
            %--------------------------------------------------------------
            dx = 1;
            
            % cross spectral density if required
            %--------------------------------------------------------------
            if nargout > 1
                Q  = spm_L_priors(N.dipfit,Q);
                Q  = spm_ssr_priors(Q);
            end
            if nargout > 2;
                [g,w] = spm_csd_mtf(Q,N);
            end
            
        end
        
        % if neuronal activity is not at its fixed point
        %------------------------------------------------------------------
        if any(abs(dx) > TOL)
            dx      = J*N.f(x,u,Q,N);
            x       = x + dx;
            [p,q,r] = N.f(x,u,Q,N,'activity');
        end
        
        % presynaptic neuronal activity at this time point
        %------------------------------------------------------------------
        pf(:,:,t) = p;             % intrinsic (inhibitory) input
        qf(:,:,t) = q;             % intrinsic (excitatory) input
        rf(:,:,t) = r;             % extrinsic (excitatory) input
        
        % electrophysiological predictions
        %------------------------------------------------------------------
        if nargout > 1
            lfp(:,t) = spm_gx_erp(x,u,Q,N);
        end
        
        % induced responses
        %------------------------------------------------------------------
        if nargout > 2
            for i = 1:M.l
                psd(t,:,i) = real(g{1}(:,i,i));
            end
        end
    end
end

% spectral responses every TR (with speckle smoothing)
%--------------------------------------------------------------------------
if nargout > 2
    j     = fix(linspace(1,t,M.ns));
    for i = 1:M.l
        for t = 1:M.ns
            csd(t,:,i) = spm_morlet_conv(psd(j(t),:,i),w,U.dt,16);
        end
    end
end


% neurovascular input in terms of intrinsic and extrinsic afferents
%==========================================================================
for i = 1:size(pf,1)
    a(:,i) = squeeze(pf(i,:,:))'*P.J(:,1) + ...
             squeeze(qf(i,:,:))'*P.J(:,2) + ...
             squeeze(rf(i,:,:))'*P.J(:,3);
end
U.u   = a;

% integrate dynamics and generate prediction
%==========================================================================

% haemodynamic model
%--------------------------------------------------------------------------
H.f   = @spm_fx_hdm;
H.g   = @spm_gx_hdm;
H.x   = M.x.h;
H.m   = M.m;
H.l   = M.l;
H.ns  = M.ns;

% solve for haemodynamic responses
%--------------------------------------------------------------------------
y        = spm_int(P.H,H,U);









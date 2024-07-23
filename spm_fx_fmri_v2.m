function [f,dfdx,D,dfdu] = spm_fx_fmri_v2(x,u,P,M)
% State equation for a dynamic [bilinear/nonlinear/Balloon] model of fMRI
% responses
%
% FORMAT [f,dfdx,D,dfdu] = spm_fx_fmri(x,u,P,M)
% x      - state vector
%   x(:,1) - excitatory neuronal activity            ue
%   x(:,2) - vascular signal                          s
%   x(:,3) - rCBF (inflow)                         ln(f)
%   x(:,4) - venous volume                         ln(v)
%   x(:,5) - deoyxHb                               ln(q)
%  [x(:,6) - inhibitory neuronal activity            ui
%
% f      - dx/dt
% dfdx   - df/dx
% dfdu   - df/du
% D      - delays
% 
% Updates relative to spm_fx_fmri.m:
% 1. Clarification of variable names


% Options
%--------------------------------------------------------------------------
if nargin < 4, M = struct([]); end

% Is this a one-state DCM (5 hidden states per region of which one is 
% neural)? Otherwise it's a two-state DCM.
one_state = (size(x,2) == 5);

% Unpack states
%--------------------------------------------------------------------------
ue  = x(:,1);
s   = x(:,2);
cbf = x(:,3);
v   = x(:,4);
q   = x(:,5);
if ~one_state
    ui = x(:,6);
end

% exponentiation of haemodynamic states
cbf = exp(cbf);
v   = exp(v);
q   = exp(q);

% Neural model
% -------------------------------------------------------------------------
if one_state 
        
    % if P.A encodes the eigenvalues of the (average) connectivity matrix
    if isvector(P.A) && numel(P.A) > 1
        
        % excitatory connections
        EE = spm_dcm_fmri_mode_gen(P.A,M.modes);
        
        % input dependent modulation
        for i = 1:size(P.B,3)
            EE = EE + u(i)*P.B(:,:,i);
        end
        
        % and nonlinear (state) terms
        for i = 1:size(P.D,3)
            EE = EE + ue(i)*P.D(:,:,i);
        end
        
    else
        
        % input dependent modulation
        for i = 1:size(P.B,3)
            P.A(:,:,1) = P.A(:,:,1) + u(i)*P.B(:,:,i);
        end
        
        % and nonlinear (state) terms
        for i = 1:size(P.D,3)
            P.A(:,:,1) = P.A(:,:,1) + ue(i)*P.D(:,:,i);
        end
        
        % combine forward and backward connections if necessary
        if size(P.A,3) > 1
            P.A  = exp(P.A(:,:,1)) - exp(P.A(:,:,2));
        end
        
        % one neuronal state per region: diag(A) is a log self-inhibition
        SE = diag(P.A);
        EE = P.A - diag(exp(SE)/2 + SE);     
    end
    
    % flow
    f(:,1) = EE*ue + P.C*u(:);
    
else
    % otherwise two neuronal states per region
    
    % input dependent modulation
    for i = 1:size(P.B,3)
        P.A(:,:,1) = P.A(:,:,1) + u(i)*P.B(:,:,i);
    end
    
    % and nonlinear (state) terms
    for i = 1:size(P.D,3)
        P.A(:,:,1) = P.A(:,:,1) + ue(i)*P.D(:,:,i);
    end
    
    % extrinsic (two neuronal states): enforce positivity
    n     = size(P.A,1);            % number of regions
    EE    = exp(P.A(:,:,1))/8;
    IE    = diag(diag(EE));         % intrinsic inhibitory to excitatory
    EE    = EE - IE;                % extrinsic excitatory to excitatory
    EI    = eye(n,n);               % intrinsic excitatory to inhibitory
    SE    = eye(n,n)/2;             % intrinsic self-inhibition (excitatory)
    SI    = eye(n,n);               % intrinsic self-inhibition (inhibitory)
    
    % excitatory proportion
    if size(P.A,3) > 1
        phi = spm_phi(P.A(:,:,2)*2);
        EI  = EI + EE.*(1 - phi);
        EE  = EE.*phi - SE;
    else
        EE  = EE - SE;
    end
    
    % motion - excitatory and inhibitory: f = dx/dt
    f(:,1) = EE*ue - IE*ui + P.C*u(:);
    f(:,6) = EI*ue - SI*ui;
    
end

% Hemodynamic motion
%==========================================================================

% hemodynamic parameters
%--------------------------------------------------------------------------
%   DECAY - signal decay (Hz)                              d(ds/dt)/ds)
%   FEEDBACK - autoregulation (Hz)                         d(ds/dt)/df)
%   TRANSIT - transit time (secs)                          (t0)
%   ALPHA - exponent for Fout(v)                           (alpha)
%   E0 - resting oxygen extraction                         (E0)
%   EPSILON - ratio of intra- to extra-vascular components (epsilon)
%          of the gradient echo signal
%--------------------------------------------------------------------------

% Parameters with constant values
DECAY    = 0.64;
FEEDBACK = 0.32; 
TRANSIT  = 2.00;
ALPHA    = 0.32;
E0       = 0.4; % Should this be 0.04?

% signal decay
sd = DECAY*exp(P.decay);

% transit time
tt = TRANSIT*exp(P.transit);

% Fout = f(v) - outflow
fv = v.^(1/ALPHA);

% e = f(f) - oxygen extraction
ff = (1 - (1 - E0).^(1./cbf))/E0;

% implement differential state equation f = dx/dt (haemodynamic)
f(:,2)   = ue - sd.*s - FEEDBACK*(cbf - 1); % vascular signal
f(:,3)   = s./cbf;                          % rCBF
f(:,4)   = (cbf - fv)./(tt.*v);             % volume
f(:,5)   = (ff.*cbf - fv.*q./v)./(tt.*q);   % dHb
f        = f(:);

if nargout < 2, return, end

% Calculate derivatives
% -------------------------------------------------------------------------

% Neuronal Jacobian
[n,m] = size(x);
if one_state
    % one neuronal state per region
    dfdx{1,1} = EE;
    for i = 1:size(P.D,3)
        D  = P.D(:,:,i) + diag((diag(EE) - 1).*diag(P.D(:,:,i)));
        dfdx{1,1}(:,i) = dfdx{1,1}(:,i) + D*ue;
    end  
else
    % two neuronal states: NB nonlinear (D) effects not implemented)
    dfdx{1,1} = EE;
    dfdx{1,6} = - IE;
    dfdx{6,1} = EI;
    dfdx{6,6} = - SI;
end

% input
dfdu{1,1} = P.C;
for i = 1:size(P.B,3)
    B  = P.B(:,:,i) + diag((diag(EE) - 1).*diag(P.B(:,:,i)));
    dfdu{1,1}(:,i) = dfdu{1,1}(:,i) + B*ue;
end
dfdu{2,1} = sparse(n*(m - 1),length(u(:)));

% Hemodynamic Jacobian
dfdx{2,1} = speye(n,n);
dfdx{2,2} = speye(n,n)*(-sd);
dfdx{2,3} = diag(-FEEDBACK*cbf);
dfdx{3,2} = diag(1./cbf);
dfdx{3,3} = diag(-s./cbf);
dfdx{4,3} = diag( cbf./(tt.*v));
dfdx{4,4} = diag(-v.^(1/ALPHA - 1)./(tt*ALPHA) - (1./v.*(cbf - v.^(1/ALPHA)))./tt);
dfdx{5,3} = diag((cbf + log(1 - E0).*(1 - E0).^(1./cbf) - cbf.*(1 - E0).^(1./cbf))./(tt.*q*E0));
dfdx{5,4} = diag((v.^(1/ALPHA - 1)*(ALPHA - 1))./(tt*ALPHA));
dfdx{5,5} = diag((cbf./q).*((1 - E0).^(1./cbf) - 1)./(tt*E0));

% concatenate
dfdx      = spm_cat(dfdx);
dfdu      = spm_cat(dfdu);
D         = 1;

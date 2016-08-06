function [f,dfdx,D,dfdu] = spm_fx_fmri(x,u,P,M)
% State equation for a dynamic [bilinear/nonlinear/Balloon] model of fMRI
% responses
% FORMAT [f,dfdx,D,dfdu] = spm_fx_fmri(x,u,P,M)
% x      - state vector
%   x(:,1) - excitatory neuronal activity            ue
%   x(:,2) - vascular signal                          s
%   x(:,3) - rCBF                                  ln(f)
%   x(:,4) - venous volume                         ln(v)
%   x(:,5) - deoyxHb                               ln(q)
%  [x(:,6) - inhibitory neuronal activity             ui
%
% f      - dx/dt
% dfdx   - df/dx
% dfdu   - df/du
% D      - delays
%
%__________________________________________________________________________
%
% References for hemodynamic & neuronal state equations:
% 1. Buxton RB, Wong EC & Frank LR. Dynamics of blood flow and oxygenation
%    changes during brain activation: The Balloon model. MRM 39:855-864,
%    1998.
% 2. Friston KJ, Mechelli A, Turner R, Price CJ. Nonlinear responses in
%    fMRI: the Balloon model, Volterra kernels, and other hemodynamics.
%    Neuroimage 12:466-477, 2000.
% 3. Stephan KE, Kasper L, Harrison LM, Daunizeau J, den Ouden HE,
%    Breakspear M, Friston KJ. Nonlinear dynamic causal models for fMRI.
%    Neuroimage 42:649-662, 2008.
% 4. Marreiros AC, Kiebel SJ, Friston KJ. Dynamic causal modelling for
%    fMRI: a two-state model.
%    Neuroimage. 2008 Jan 1;39(1):269-78.
%__________________________________________________________________________
% Copyright (C) 2002-2014 Wellcome Trust Centre for Neuroimaging

% Karl Friston & Klaas Enno Stephan
% $Id: spm_fx_fmri.m 6855 2016-08-06 10:06:35Z karl $

% options
%--------------------------------------------------------------------------
if nargin < 4, M = struct([]); end
if isfield(M,'symmetry'), symmetry = M.symmetry; else symmetry = 0; end


% Neuronal motion
%==========================================================================
P.A   = full(P.A);                       %    linear parameters
P.B   = full(P.B);                       % bi-linear parameters
P.C   = P.C/16;                          % exogenous parameters
P.D   = full(P.D);                       % nonlinear parameters


% implement differential state equation y = dx/dt (neuronal)
%--------------------------------------------------------------------------
f    = x;

% if there are five hidden states per region, only one is neuronal
%==========================================================================
if size(x,2) == 5
    
    
    % if P.A encodes the eigenvalues of the (average) connectivity matrix
    %======================================================================
    if isvector(P.A)
        
        % excitatory connections
        %------------------------------------------------------------------
        EE = spm_dcm_fmri_mode_gen(P.A,M.modes);
        
        % input dependent modulation
        %------------------------------------------------------------------
        for i = 1:size(P.B,3)
            EE = EE + u(i)*P.B(:,:,i);
        end
        
        % and nonlinear (state) terms
        %------------------------------------------------------------------
        for i = 1:size(P.D,3)
            EE = EE + x(i,1)*P.D(:,:,i);
        end
        
    else % otherwise average connections are encoded explicitly
        %==================================================================
        
        % input dependent modulation
        %------------------------------------------------------------------
        for i = 1:size(P.B,3)
            P.A(:,:,1) = P.A(:,:,1) + u(i)*P.B(:,:,i);
        end
        
        % and nonlinear (state) terms
        %------------------------------------------------------------------
        for i = 1:size(P.D,3)
            P.A(:,:,1) = P.A(:,:,1) + x(i,1)*P.D(:,:,i);
        end
        
        % combine forward and backward connections if necessary
        %------------------------------------------------------------------
        if size(P.A,3) > 1
            P.A  = exp(P.A(:,:,1)) - exp(P.A(:,:,2));
        end
        
        % one neuronal state per region: diag(A) is a log self-inhibition
        %------------------------------------------------------------------
        SE     = diag(P.A);
        EE     = P.A - diag(exp(SE)/2 + SE);
        
        % symmetry constraints for demonstration purposes
        %------------------------------------------------------------------
        if symmetry, EE = (EE + EE')/2; end
        
    end
    
    % flow
    %----------------------------------------------------------------------
    f(:,1) = EE*x(:,1) + P.C*u(:);
    
    
else
    
    % otherwise two neuronal states per region
    %======================================================================
    
    % input dependent modulation
    %----------------------------------------------------------------------
    for i = 1:size(P.B,3)
        P.A(:,:,1) = P.A(:,:,1) + u(i)*P.B(:,:,i);
    end
    
    % and nonlinear (state) terms
    %----------------------------------------------------------------------
    for i = 1:size(P.D,3)
        P.A(:,:,1) = P.A(:,:,1) + x(i,1)*P.D(:,:,i);
    end
    
    % extrinsic (two neuronal states): enforce positivity
    %----------------------------------------------------------------------
    n     = size(P.A,1);            % number of regions
    EE    = exp(P.A(:,:,1))/8;
    IE    = diag(diag(EE));         % intrinsic inhibitory to excitatory
    EE    = EE - IE;                % extrinsic excitatory to excitatory
    EI    = eye(n,n);               % intrinsic excitatory to inhibitory
    SE    = eye(n,n)/2;             % intrinsic self-inhibition (excitatory)
    SI    = eye(n,n);               % intrinsic self-inhibition (inhibitory)
    
    % excitatory proportion
    %----------------------------------------------------------------------
    if size(P.A,3) > 1
        phi = spm_phi(P.A(:,:,2)*2);
        EI  = EI + EE.*(1 - phi);
        EE  = EE.*phi - SE;
    else
        EE  = EE - SE;
    end
    
    % motion - excitatory and inhibitory: f = dx/dt
    %----------------------------------------------------------------------
    f(:,1) = EE*x(:,1) - IE*x(:,6) + P.C*u(:);
    f(:,6) = EI*x(:,1) - SI*x(:,6);
    
end

% Hemodynamic motion
%==========================================================================

% hemodynamic parameters
%--------------------------------------------------------------------------
%   H(1) - signal decay                                   d(ds/dt)/ds)
%   H(2) - autoregulation                                 d(ds/dt)/df)
%   H(3) - transit time                                   (t0)
%   H(4) - exponent for Fout(v)                           (alpha)
%   H(5) - resting oxygen extraction                      (E0)
%   H(6) - ratio of intra- to extra-vascular components   (epsilon)
%          of the gradient echo signal
%--------------------------------------------------------------------------
H        = [0.64 0.32 2.00 0.32 0.4];

% exponentiation of hemodynamic state variables
%--------------------------------------------------------------------------
x(:,3:5) = exp(x(:,3:5));

% signal decay
%--------------------------------------------------------------------------
sd       = H(1)*exp(P.decay);

% transit time
%--------------------------------------------------------------------------
tt       = H(3)*exp(P.transit);

% Fout = f(v) - outflow
%--------------------------------------------------------------------------
fv       = x(:,4).^(1/H(4));

% e = f(f) - oxygen extraction
%--------------------------------------------------------------------------
ff       = (1 - (1 - H(5)).^(1./x(:,3)))/H(5);


% implement differential state equation f = dx/dt (hemodynamic)
%--------------------------------------------------------------------------
f(:,2)   = x(:,1) - sd.*x(:,2) - H(2)*(x(:,3) - 1);
f(:,3)   = x(:,2)./x(:,3);
f(:,4)   = (x(:,3) - fv)./(tt.*x(:,4));
f(:,5)   = (ff.*x(:,3) - fv.*x(:,5)./x(:,4))./(tt.*x(:,5));
f        = f(:);


if nargout < 2, return, end


% Neuronal Jacobian
%==========================================================================
[n,m] = size(x);
if m == 5
    
    % one neuronal state per region
    %----------------------------------------------------------------------
    dfdx{1,1} = EE;
    for i = 1:size(P.D,3)
        D  = P.D(:,:,i) + diag((diag(EE) - 1).*diag(P.D(:,:,i)));
        dfdx{1,1}(:,i) = dfdx{1,1}(:,i) + D*x(:,1);
    end
    
else
    
    % two neuronal states: NB nonlinear (D) effects not implemented)
    %----------------------------------------------------------------------
    dfdx{1,1} = EE;
    dfdx{1,6} = - IE;
    dfdx{6,1} = EI;
    dfdx{6,6} = - SI;
    
end

% input
%==========================================================================
dfdu{1,1} = P.C;
for i = 1:size(P.B,3)
    B  = P.B(:,:,i) + diag((diag(EE) - 1).*diag(P.B(:,:,i)));
    dfdu{1,1}(:,i) = dfdu{1,1}(:,i) + B*x(:,1);
end
dfdu{2,1} = sparse(n*(m - 1),length(u(:)));


% Hemodynamic Jacobian
%==========================================================================
dfdx{2,1} = speye(n,n);
dfdx{2,2} = diag(-sd);
dfdx{2,3} = diag(-H(2)*x(:,3));
dfdx{3,2} = diag( 1./x(:,3));
dfdx{3,3} = diag(-x(:,2)./x(:,3));
dfdx{4,3} = diag( x(:,3)./(tt.*x(:,4)));
dfdx{4,4} = diag(-x(:,4).^(1/H(4) - 1)./(tt*H(4)) - (1./x(:,4).*(x(:,3) - x(:,4).^(1/H(4))))./tt);
dfdx{5,3} = diag((x(:,3) + log(1 - H(5)).*(1 - H(5)).^(1./x(:,3)) - x(:,3).*(1 - H(5)).^(1./x(:,3)))./(tt.*x(:,5)*H(5)));
dfdx{5,4} = diag((x(:,4).^(1/H(4) - 1)*(H(4) - 1))./(tt*H(4)));
dfdx{5,5} = diag((x(:,3)./x(:,5)).*((1 - H(5)).^(1./x(:,3)) - 1)./(tt*H(5)));


% concatenate
%--------------------------------------------------------------------------
dfdx      = spm_cat(dfdx);
dfdu      = spm_cat(dfdu);
D         = 1;

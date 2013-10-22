function [f,dfdx,D,dfdu] = spm_fx_fmri(x,u,P,M)
% state equation for a dynamic [bilinear/nonlinear/Balloon] model of fMRI
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
%___________________________________________________________________________
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

% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% Karl Friston & Klaas Enno Stephan
% $Id: spm_fx_fmri.m 5708 2013-10-22 09:20:59Z karl $


% Neuronal motion
%==========================================================================
P.A   = full(P.A);                       %    linear parameters
P.B   = full(P.B);                       % bi-linear parameters
P.C   = P.C/16;                          % exogenous parameters
P.D   = full(P.D);                       % nonlinear parameters

% excitatory connections
%--------------------------------------------------------------------------
for i = 1:size(P.B,3)
    P.A(:,:,1) = P.A(:,:,1) + u(i)*P.B(:,:,i);
end

% and nonlinear (state) terms
%--------------------------------------------------------------------------
for i = 1:size(P.D,3)
    P.A(:,:,1) = P.A(:,:,1) + x(i,1)*P.D(:,:,i);
end

% implement differential state equation y = dx/dt (neuronal)
%--------------------------------------------------------------------------
f    = x;
if size(x,2) == 5
    
    % combine forward and backward connections if necessary
    %--------------------------------------------------------------------------
    if size(P.A,3) > 1
        P.A  = exp(P.A(:,:,1)) - exp(P.A(:,:,2));
    end
    
    % one neuronal state per region: diag(A) is a log self-inhibition
    %----------------------------------------------------------------------
    SI     = diag(P.A);
    P.A    = P.A - diag(exp(SI)/2 + SI);
    
    % flow
    %----------------------------------------------------------------------
    f(:,1) = P.A*x(:,1) + P.C*u(:);
    
else
    
    % extrinsic (two neuronal states): enforce positivity
    %----------------------------------------------------------------------
    P.A   = exp(P.A(:,:,1))/8;
    DA    = diag(diag(P.A));        % intrinsic connectivity
    EE    = P.A - DA;               % excitatory to excitatory
    n     = length(P.A);            % number of regions
    IE    = eye(n,n);               % inhibitory to excitatory
    EI    = eye(n,n);               % excitatory to inhibitory
    SE    = eye(n,n);               % self-inhibition (excitatory)
    SI    = eye(n,n);               % self-inhibition (inhibitory)
    
    
    % <<< intrinsic connectivity >>>
    %----------------------------------------------------------------------
    IE    = DA;                     % self-inhibition (inhibitory)
    
    % <<< switch excitatory to excitatory -> excitatory to inhibitory >>>
    %----------------------------------------------------------------------
    in    = {};
    for i = 1:length(in)
        EI(in{i}(1),in{i}(2)) = EE(in{i}(1),in{i}(2));
        EE(in{i}(1),in{i}(2)) = 0;
    end
    
    % motion - excitatory and inhibitory: f = dx/dt
    %----------------------------------------------------------------------
    f(:,1) = (EE - SE)*x(:,1) - IE*x(:,6) + P.C*u(:);
    f(:,6) = EI*x(:,1) - SI*2*x(:,6);
    
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
H        = [0.64 0.32 2.00 0.32 0.32];

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
[n m] = size(x);
if m == 5
    
    % one neuronal state per region
    %----------------------------------------------------------------------
    dfdx{1,1} = P.A;
    for i = 1:size(P.D,3)
        D  = P.D(:,:,i) + diag((diag(P.A) - 1).*diag(P.D(:,:,i)));
        dfdx{1,1}(:,i) = dfdx{1,1}(:,i) + D*x(:,1);
    end
    
else
    
    % two neuronal states: NB nonlinear (D) effects not implemented)
    %----------------------------------------------------------------------
    dfdx{1,1} = EE - SE;
    dfdx{1,6} = - IE;
    dfdx{6,1} = EI;
    dfdx{6,6} = - SI*2;
    
end

% input
%==========================================================================
dfdu{1,1} = P.C;
for i = 1:size(P.B,3)
    B  = P.B(:,:,i) + diag((diag(P.A) - 1).*diag(P.B(:,:,i)));
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








function [f] = spm_fx_fnirs(x,u,P,M)
% State equation for a dynamic model of fNIRS responses
% FORMAT [f] = spm_fx_nirs(x,u,P,M)
%
% x      - state vector
%--------------------------------------------------------------------------
%   x(:,1) - excitatory neuronal activity            ue
%   x(:,2) - vasodilatory signal                          s
%   x(:,3) - rCBF                                  ln(f)
%   x(:,4) - venous volume                         ln(v)
%   x(:,5) - deoxyHb                               ln(q)
%   x(:,6) - totalHb                               ln(p)
%   [x(:,7)] - inhibitory neuronal activity            ui
%--------------------------------------------------------------------------
% u     experimental inputs 
% P     prior of latent variables 
% M    model structure
%
% f      - dx/dt
%___________________________________________________________________________
%
% References for hemodynamic & neuronal state equations:
% 1. Friston KJ, Mechelli A, Turner R, Price CJ. Nonlinear responses in
%    fMRI: the Balloon model, Volterra kernels, and other hemodynamics.
%    Neuroimage 12:466-477, 2000.
% 2. Stephan KE, Kasper L, Harrison LM, Daunizeau J, den Ouden HE,
%    Breakspear M, Friston KJ. Nonlinear dynamic causal models for fMRI.
%    Neuroimage 42:649-662, 2008.
% 3. Marreiros AC, Kiebel SJ, Friston KJ. Dynamic causal modelling for
%    fMRI: a two-state model.
%    Neuroimage. 39(1):269-78, 2008. 
% 4. Buxton RB, Uludag, K, Dubowitz, DJ, Liu, TT. Modeling the hemodynamic
%    response to brain activation. Neuroimage. 2004, 23: 220-233. 
% 5. X Cui and S Bray and A Reiss. Functional near infrared spectroscopy (NIRS)
%    signal improvement based on negative correlation between oxygenated and
%    deoxygenated hemoglobin dynamics.
%    Neuroimage 49:3039-3046, 2010.
%
% This script is based on spm_fx_fmri.m written by 
% Karl Friston & Klaas Enno Stephan.
%__________________________________________________________________________
% Copyright (C) 2015 Wellcome Trust Centre for Neuroimaging

% Will Penny & Sungho Tak
% $Id: spm_fx_fnirs.m 6942 2016-11-21 13:17:44Z guillaume $

% Neuronal motion
%==========================================================================
P.A   = full(P.A);                       %    linear parameters
P.B   = full(P.B);                       % bi-linear parameters
P.C   = P.C/16;                          % exogenous parameters

% implement differential state equation y = dx/dt (neuronal)
%--------------------------------------------------------------------------
f    = x;

% input dependent modulation
%--------------------------------------------------------------------------
for i = 1:size(P.B,3)
    P.A(:,:,1) = P.A(:,:,1) + u(i)*P.B(:,:,i);
end

if size(x, 2) == 6 % one neuronal state per region
    %======================================================================
    
    % one neuronal state per region: diag(A) is a log self-inhibition
    %----------------------------------------------------------------------
    SI     = diag(P.A);
    P.A    = P.A - diag(exp(SI)/2 + SI);
    
    % flow
    %----------------------------------------------------------------------
    f(:,1) = P.A*x(:,1) + P.C*u(:);
    
elseif size(x,2) == 7 % two neuronal states per region
    %======================================================================
    
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
    f(:,1) = EE*x(:,1) - IE*x(:, 7) + P.C*u(:);
    f(:,6 + 1) = EI*x(:,1) - SI*x(:, 7);
end

% Hemodynamic motion
%==========================================================================

% hemodynamic parameters
%--------------------------------------------------------------------------
%   H(1) - signal decay                                   d(ds/dt)/ds)
%   H(2) - autoregulation                                 d(ds/dt)/df)
%   H(3) - transit time                                      (t0)
%   H(4) - exponent for Fout(v)                           (alpha)
%   H(5) - resting oxygen extraction                      (E0)
%   H(6) - viscoelastic time constant              (tv)
%--------------------------------------------------------------------------
H = [0.64 0.32 2.00 0.32 0.32 2.00];

% exponentiation of hemodynamic state variables
%--------------------------------------------------------------------------
x(:,3:6) = exp(x(:,3:6));

% signal decay
sd = H(1)*exp(P.decay);

% autoregulatory feedback
af = H(2).*exp(P.afback);

% transit time
tt = H(3).*exp(P.transit);
% viscoelastic time constant 
tv = H(6).*exp(P.tv); 

% outflow f(v) (steady state)
fv_s = x(:,4).^(1/H(4));

% implement differential state equation f = dx/dt (hemodynamic)
%--------------------------------------------------------------------------
f(:,2)   = x(:,1) - sd.*x(:,2) - af.*(x(:,3) - 1);
f(:,3)   = x(:,2)./x(:,3);
f(:,4)   = (x(:,3) - fv_s)./((tt+tv).*x(:,4));

% outflow (dynamic state) 
fv_d = fv_s + tv.*x(:,4).*f(:,4);

% oxygen extraction fraction 
ff = (1 - (1 - H(5)).^(1./x(:,3)))/H(5);

f(:,5)   = (x(:,3).*ff - fv_d.*x(:,5)./x(:,4))./(tt.*x(:,5));
f(:,6)   = (x(:,3) - fv_d) ./ (tt.*x(:,4));

f        = f(:);







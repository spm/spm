function [y] = spm_fx_fmri(x,u,P,M)
% state equation for a dynamic [bilinear/nonlinear/Balloon] model of fMRI
% responses
% FORMAT [y] = spm_fx_fmri(x,u,P,M)
% x      - state vector
%   x(:,1) - excitatory neuronal activity             ue
%   x(:,2) - vascular signal                          s
%   x(:,3) - rCBF                                  ln(f)
%   x(:,4) - venous volume                         ln(v)
%   x(:,5) - deoyxHb                               ln(q)
%  [x(:,6) - inhibitory neuronal activity             ui]
%
% y      - dx/dt
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
% $Id: spm_fx_fmri.m 5457 2013-04-30 14:13:20Z karl $


% Neuronal motion
%==========================================================================
P.B   = full(P.B);                       % bi-linear parameters
P.C   = P.C/16;                          % exogenous parameters
P.D   = full(P.D);                       % nonlinear parameters

% excitatory connections
%--------------------------------------------------------------------------
for i = 1:size(P.B,3)
    P.A = P.A + u(i)*P.B(:,:,i);
end

% and nonlinear (state) terms
%--------------------------------------------------------------------------
for i = 1:size(P.D,3)
    P.A = P.A + x(i,1)*P.D(:,:,i);
end

% implement differential state equation y = dx/dt (neuronal)
%--------------------------------------------------------------------------
y    = x;
if size(x,2) == 5  
    
    % one neuronal state per region: diag(A) is a log self-inhibition
    %----------------------------------------------------------------------
    SI     = diag(P.A);
    P.A    = P.A - diag(exp(SI)/2 + SI);
    y(:,1) = P.A*x(:,1) + P.C*u(:);

else

    % extrinsic (two neuronal states)
    %----------------------------------------------------------------------
    n     = length(P.A);            % number of regions
    A     = exp(P.A)/8;             % enforce positivity
    DA    = diag(diag(A));          % intrinsic connectivity
    EE    = A - DA;                 % excitatory to excitatory
    IE    = eye(n,n);               % inhibitory to excitatory
    EI    = eye(n,n);               % excitatory to inhibitory
    SE    = eye(n,n);               % self-inhibition (excitatory)
    SI    = eye(n,n);               % self-inhibition (inhibitory)
    
    
    % <<< intrinsic moduation >>>
    %----------------------------------------------------------------------
    IE    = DA;                     % self-inhibition (inhibitory)

    % <<< switch excitatory to excitatory -> excitatory to inhibitory >>>
    %----------------------------------------------------------------------
    in    = {};
    for i = 1:length(in)
        EI(in{i}(1),in{i}(2)) = EE(in{i}(1),in{i}(2));
        EE(in{i}(1),in{i}(2)) = 0;
    end
    

    % motion - excitatory and inhibitory: y = dx/dt
    %----------------------------------------------------------------------
    y(:,1) = EE*x(:,1) - SE*x(:,1) - IE*x(:,6) + P.C*u(:);
    y(:,6) = EI*x(:,1) - SI*x(:,6)*2;

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


% implement differential state equation y = dx/dt (hemodynamic)
%--------------------------------------------------------------------------
y(:,2)   = x(:,1) - sd.*x(:,2) - H(2)*(x(:,3) - 1);
y(:,3)   = x(:,2)./x(:,3);
y(:,4)   = (x(:,3) - fv)./(tt.*x(:,4));
y(:,5)   = (ff.*x(:,3) - fv.*x(:,5)./x(:,4))./(tt.*x(:,5));
y        = y(:);

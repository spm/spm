function [f,J] = spm_fx_wendling(x,u,P,M)
% State equations for the Wendling neural mass model (hippocampal CA1)
% FORMAT [f,J] = spm_fx_wendling(x,u,P,M)
% x      - state vector (10 states per source)
%   x(:,1)  - excitatory PSP on excitatory interneurons (y0)  [from pyramidal]
%   x(:,2)  - excitatory PSP on pyramidal cells (y1)          [from E + input]
%   x(:,3)  - slow inhibitory PSP on pyramidal cells (y2)     [from SDI]
%   x(:,4)  - fast inhibitory PSP on pyramidal cells (y3)     [from FSI]
%   x(:,5)  - excitatory PSP on SDI interneurons (y4)         [from pyramidal]
%   x(:,6)  - derivative of y0
%   x(:,7)  - derivative of y1
%   x(:,8)  - derivative of y2
%   x(:,9)  - derivative of y3
%   x(:,10) - derivative of y4
%
% u      - exogenous input
% P      - parameter structure
% M      - model structure
%
% f    = dx(t)/dt  = f(x(t))
% J    = df/dx
%
% Output (observed LFP) = y1 - y2 - y3 = x(:,2) - x(:,3) - x(:,4)
%
% Four populations:
%   Pyramidal cells (P)                - output
%   Excitatory interneurons (E)        - feedback excitation
%   Slow dendritic inhibitory (SDI)    - GABAa,slow (dendritic)
%   Fast somatic inhibitory (FSI)      - GABAa,fast (somatic)
%
% Three synaptic kernels:
%   h_e(t)  = A * a * t * exp(-a*t)   excitatory (AMPA), a = 100/s
%   h_i(t)  = B * b * t * exp(-b*t)   slow inhibition,   b = 50/s
%   h_g(t)  = G * g * t * exp(-g*t)   fast inhibition,   g = 500/s
%
% Connectivity C1-C7 (all active):
%   C1 = 135    Pyramidal -> Excitatory interneurons
%   C2 = 108    Excitatory interneurons -> Pyramidal
%   C3 = 33.75  Pyramidal -> SDI
%   C4 = 33.75  SDI -> Pyramidal
%   C5 = 40.5   Pyramidal -> FSI (excitatory input)
%   C6 = 13.5   SDI -> FSI (inhibitory input)
%   C7 = 108    FSI output gain
%
% Equations (Wendling et al., 2002, Eur J Neurosci 15:1499-1508):
%
%   y0'' = Aa*S(vP) - 2a*y0' - a^2*y0
%   y1'' = Aa*(p(t) + C2*S(C1*y0)) - 2a*y1' - a^2*y1
%   y2'' = Bb*C4*S(C3*y0) - 2b*y2' - b^2*y2
%   y3'' = Gg*C7*S(C5*y0 - C6*y4) - 2g*y3' - g^2*y3
%   y4'' = Bb*S(C3*y0) - 2b*y4' - b^2*y4
%
%   vP   = y1 - y2 - y3   (pyramidal membrane potential = output)
%   vFSI = C5*y0 - C6*y4  (FSI membrane potential)
%
%__________________________________________________________________________
%
% Wendling F, Bartolomei F, Bellanger JJ, Chauvel P (2002) Epileptic fast
% activity can be explained by a model of impaired GABAergic dendritic
% inhibition. Eur J Neurosci 15:1499-508
%__________________________________________________________________________

% get dimensions and configure state variables
%--------------------------------------------------------------------------
x    = spm_unvec(x,M.x);
n    = size(x,1);              % number of sources
s    = size(x,2);              % number of states (10)

% [default] fixed parameters (Wendling 2002, Table 1)
%--------------------------------------------------------------------------
E    = [32 16 4];              % extrinsic rates (forward, backward, lateral)
A0   = 3.25;                   % excitatory synaptic gain A (mV)
B0   = 22;                     % slow dendritic inhibitory gain B (mV)
G0   = 10;                     % fast somatic inhibitory gain G (mV)
T0   = [10 20 2];              % time constants (ms): 1/a, 1/b, 1/g
C0   = [135 108 33.75 33.75 40.5 13.5 108]; % connectivity C1-C7
R0   = [0.56 6];               % sigmoid parameters [slope, threshold]
e0   = 2.5;                    % half max firing rate
D0   = [2 4];                  % propagation delays (intrinsic, extrinsic) ms

% [specified] fixed parameters
%--------------------------------------------------------------------------
if isfield(M,'pF')
    try, E  = M.pF.E;  end
    try, A0 = M.pF.A;  end
    try, B0 = M.pF.B;  end
    try, G0 = M.pF.G;  end
    try, T0 = M.pF.T;  end
    try, C0 = M.pF.C;  end
    try, R0 = M.pF.R;  end
    try, e0 = M.pF.e0; end
    try, D0 = M.pF.D;  end
end

% exponential transform for extrinsic connectivity
%--------------------------------------------------------------------------
A{1} = exp(P.A{1})*E(1);
A{2} = exp(P.A{2})*E(2);
A{3} = exp(P.A{3})*E(3);
Ce   = exp(P.C);

% intrinsic gains (log-deviation scaling)
%--------------------------------------------------------------------------
EXC = A0*exp(P.E(:));         % excitatory gain per source (n x 1)
SDI = B0*exp(P.S(:));         % slow inhibitory gain per source
FSI = G0*exp(P.F(:));         % fast inhibitory gain per source

% time constants (ms -> seconds)
%--------------------------------------------------------------------------
Te  = T0(1)/1000*exp(P.T(:,1));   % excitatory (1/a)
Tsi = T0(2)/1000*exp(P.T(:,2));   % slow inhibitory (1/b)
Tfi = T0(3)/1000*exp(P.T(:,3));   % fast inhibitory (1/g)

% intrinsic connectivity: Ci(source, connection)
%--------------------------------------------------------------------------
Ci = zeros(n,7);
for k = 1:7
    Ci(:,k) = C0(k)*exp(P.H(:,k));
end

% sigmoid parameters
%--------------------------------------------------------------------------
Rp = R0.*exp(P.R);
r  = Rp(1);
v0 = Rp(2);
S0 = 2*e0/(1 + exp(r*v0));    % offset so S(0) = 0

% transpose states for column-vector operations
%--------------------------------------------------------------------------
x  = x';                       % now s x n

% membrane potentials
%--------------------------------------------------------------------------
% Pyramidal: vP = y1 - y2 - y3
vP   = x(2,:) - x(3,:) - x(4,:);

% Excitatory interneuron input: C1*y0
vE   = Ci(:,1)' .* x(1,:);

% SDI interneuron input: C3*y0
vSDI = Ci(:,3)' .* x(1,:);

% FSI interneuron input: C5*y0 - C6*y4
vFSI = Ci(:,5)' .* x(1,:) - Ci(:,6)' .* x(5,:);

% pre-synaptic firing rates: S(v) = 2*e0/(1+exp(r*(v0-v))) - S0
%--------------------------------------------------------------------------
SP   = 2*e0./(1 + exp(r*(v0 - vP)))   - S0;
SE   = 2*e0./(1 + exp(r*(v0 - vE)))   - S0;
SSDI = 2*e0./(1 + exp(r*(v0 - vSDI))) - S0;
SFSI = 2*e0./(1 + exp(r*(v0 - vFSI))) - S0;

% sigmoid derivatives: dS/dv = r*s*(1 - s/(2*e0)) where s = S + S0
%--------------------------------------------------------------------------
dSP   = r*(SP   + S0).*(1 - (SP   + S0)/(2*e0));
dSE   = r*(SE   + S0).*(1 - (SE   + S0)/(2*e0));
dSSDI = r*(SSDI + S0).*(1 - (SSDI + S0)/(2*e0));
dSFSI = r*(SFSI + S0).*(1 - (SFSI + S0)/(2*e0));

% input
%==========================================================================
if isfield(M,'u')
    U = u(:)*32;
else
    U = Ce*u(:);
end

% State: f(x) and Jacobian dfdx
%==========================================================================
F = cell(n,1);
J = cell(n,n);

for i = 1:n

    te  = Te(i);
    tsi = Tsi(i);
    tfi = Tfi(i);
    exc = EXC(i);
    sdi = SDI(i);
    fsi = FSI(i);
    c   = Ci(i,:);

    % rate constants for second-order systems
    %----------------------------------------------------------------------
    ke  = -2/te;     Ke  = -1/te^2;
    ksi = -2/tsi;    Ksi = -1/tsi^2;
    kfi = -2/tfi;    Kfi = -1/tfi^2;

    % input gains: Gain/T for each kernel
    %----------------------------------------------------------------------
    ge  = exc/te;    % A*a
    gsi = sdi/tsi;   % B*b
    gfi = fsi/tfi;   % G*g

    % linear dynamics Jacobian
    %----------------------------------------------------------------------
    dfdx_i = zeros(s,s);

    % f(1:5) = x(6:10)  (PSP derivatives)
    dfdx_i(1,6)   = 1;
    dfdx_i(2,7)   = 1;
    dfdx_i(3,8)   = 1;
    dfdx_i(4,9)   = 1;
    dfdx_i(5,10)  = 1;

    % second-order dynamics (spring + damping)
    dfdx_i(6,1)   = Ke;     dfdx_i(6,6)   = ke;     % y0: h_e kernel
    dfdx_i(7,2)   = Ke;     dfdx_i(7,7)   = ke;     % y1: h_e kernel
    dfdx_i(8,3)   = Ksi;    dfdx_i(8,8)   = ksi;    % y2: h_i kernel
    dfdx_i(9,4)   = Kfi;    dfdx_i(9,9)   = kfi;    % y3: h_g kernel
    dfdx_i(10,5)  = Ksi;    dfdx_i(10,10) = ksi;    % y4: h_i kernel

    % sigmoid-dependent Jacobian terms
    %----------------------------------------------------------------------
    dfdS_i = zeros(s,s);

    % f(6): y0'' += ge*S(vP), vP = x(2) - x(3) - x(4)
    dfdS_i(6,2) =  ge*dSP(i);
    dfdS_i(6,3) = -ge*dSP(i);
    dfdS_i(6,4) = -ge*dSP(i);

    % f(7): y1'' += ge*C2*S(vE), vE = C1*x(1)
    dfdS_i(7,1) = ge*c(2)*c(1)*dSE(i);

    % f(8): y2'' += gsi*C4*S(vSDI), vSDI = C3*x(1)
    dfdS_i(8,1) = gsi*c(4)*c(3)*dSSDI(i);

    % f(9): y3'' += gfi*C7*S(vFSI), vFSI = C5*x(1) - C6*x(5)
    dfdS_i(9,1) =  gfi*c(7)*c(5)*dSFSI(i);
    dfdS_i(9,5) = -gfi*c(7)*c(6)*dSFSI(i);

    % f(10): y4'' += gsi*S(vSDI), vSDI = C3*x(1)
    dfdS_i(10,1) = gsi*c(3)*dSSDI(i);

    % motion f(x) for this source
    %----------------------------------------------------------------------
    fi = dfdx_i*x(:,i);
    fi(6)  = fi(6)  + ge *SP(i);
    fi(7)  = fi(7)  + ge *(c(2)*SE(i) + U(i));
    fi(8)  = fi(8)  + gsi*c(4)*SSDI(i);
    fi(9)  = fi(9)  + gfi*c(7)*SFSI(i);
    fi(10) = fi(10) + gsi*SSDI(i);

    F{i}   = fi;
    J{i,i} = dfdx_i + dfdS_i;

    % extrinsic afferents
    %----------------------------------------------------------------------
    for j = 1:n
        if i ~= j
            % extrinsic input drives y1 (excitatory PSP on pyramidal)
            k_fwd = ge*(A{1}(i,j) + A{3}(i,j));   % forward + lateral
            k_bwd = ge*(A{2}(i,j) + A{3}(i,j));   % backward + lateral
            k_ext = k_fwd + k_bwd;

            F{i}(7) = F{i}(7) + k_ext*SP(j);

            % Jacobian: df_i(7)/dx_j via vP_j = x_j(2)-x_j(3)-x_j(4)
            dfdS_ext      = zeros(s,s);
            dfdS_ext(7,2) =  k_ext*dSP(j);
            dfdS_ext(7,3) = -k_ext*dSP(j);
            dfdS_ext(7,4) = -k_ext*dSP(j);

            J{i,j} = dfdS_ext;
        end
    end
end

% construct global motion and Jacobian
%--------------------------------------------------------------------------
f    = zeros(n*s,1);
dfdx = zeros(n*s,n*s);

for i = 1:n
    k = (1:n:s*n) + (i - 1);
    f(k,1) = F{i};
    for j = 1:n
        if ~isempty(J{i,j})
            l         = (1:n:s*n) + (j - 1);
            dfdx(k,l) = J{i,j};
        end
    end
end

% extrinsic and intrinsic delays
%--------------------------------------------------------------------------
De = D0(2).*exp(P.D)/1000;
Di = D0(1).*exp(P.I)/1000;
De = (eye(n,n) - 1).*De;
Di = (eye(s,s) - 1)*Di;
De = kron(ones(s,s),De);
Di = kron(Di,eye(n,n));

D  = Di + De;

% Implement: dx(t)/dt = f(x(t + d)) = inv(1 - D.*dfdx)*f(x(t))
%--------------------------------------------------------------------------
D  = spm_inv(speye(n*s,n*s) - D.*dfdx);
f  = D*f;
J  = D*dfdx;

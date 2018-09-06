function [f,J,Q] = spm_fx_bgt(x,u,P,M)
% state equations for a neural mass model of the basal ganglia & thalamus
% [striatum, gpe, stn, gpi, and thalamus] as a
% single source (no extrinsic connections)
%
% order           cells     states
% 1 = striatum  - ii        x(1,1:2)
% 2 = gpe       - ii        x(1,3:4)
% 3 = stn       - pyr       x(1,5:6)
% 4 = gpi       - ii        x(1,7:8)
% 5 = thalamus  - pyr       x(1,9:10)
%
% G(1,1) = str -> str (-ve self)
% G(1,2) = str -> gpe (-ve ext)
% G(1,3) = gpe -> gpe (-ve self)
% G(1,4) = gpe -> stn (-ve ext)
% G(1,5) = stn -> gpe (+ve ext)
% G(1,6) = str -> gpi (-ve ext)
% G(1,7) = stn -> gpi (+ve ext)
% G(1,8) = gpi -> gpi (-ve self)
% G(1,9) = gpi -> tha (-ve ext)
%__________________________________________________________________________
% Copyright (C) 2016 Wellcome Trust Centre for Neuroimaging
 
% Bernadette van Wijk
% $Id: spm_fx_bgt.m 7412 2018-09-06 10:12:18Z guillaume $


% check if intrinsic connections are free parameters
%--------------------------------------------------------------------------
try, P.G; catch, P.G = 0; end

% get dimensions and configure state variables
%--------------------------------------------------------------------------
x  = spm_unvec(x,M.x);               % neuronal states
n  = size(x,1);                      % number of sources

% [default] fixed parameters
%--------------------------------------------------------------------------
G  = [2 2 2 2 2 2 2 2 2]*200;   % synaptic connection strengths
T  = [8 8 4 8 8];               % synaptic time constants [str,gpe,stn,gpi,tha];
R  = 2/3;                       % slope of sigmoid activation function
% NB for more pronounced state dependent transfer functions use R  = 3/2;

if isfield(M,'BGT_G'); G = M.BGT_G; end
if isfield(M,'BGT_T'); T = M.BGT_T; end


% [specified] fixed parameters
%--------------------------------------------------------------------------
if isfield(M,'pF')
    try, E = M.pF.E; end
    try, G = M.pF.G; end
    try, T = M.pF.T; end
    try, R = M.pF.R; end
end

% input connections
%--------------------------------------------------------------------------
C    = exp(P.C);

% pre-synaptic inputs: s(V)
%--------------------------------------------------------------------------
R    = R.*exp(P.S);              % gain of activation function
F    = 1./(1 + exp(-R*x + 0));   % firing rate
S    = F - 1/(1 + exp(0));       % deviation from baseline firing (0)

% input
%==========================================================================
if isfield(M,'u')

    % endogenous input
    %----------------------------------------------------------------------
    U = u(:)*128;

else
    % exogenous input
    %----------------------------------------------------------------------
    U = C*u(:)*32;
end


% time constants and intrinsic connections
%==========================================================================
T     = T/1000;
for i = 1:size(P.T,2)
    T(:,i) = T(:,i).*exp(P.T(:,i));
end


% intrinsic/extrinsic connections to be optimised
%--------------------------------------------------------------------------
j     = 1:9;
for i = 1:size(P.G,2)
    G(:,j(i)) = G(:,j(i)).*exp(P.G(:,i));
end


% Motion of states: f(x)
%--------------------------------------------------------------------------

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 1 - Str: ii

% inhibitory interneurons: Hidden states - error
%--------------------------------------------------------------------------
u      =  U;
u      =  - G(:,1)*S(:,1) + u;
f(:,2) =  (u - 2*x(:,2) - x(:,1)./T(1,1))./T(1,1);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 2 - GPe: ii

% inhibitory interneurons: Hidden states - error
%--------------------------------------------------------------------------
u      =  0;  
u      =  - G(:,2)*S(:,1) + G(:,5)*S(:,5)  - G(:,3)*S(:,3) + u;
f(:,4) =  (u - 2*x(:,4) - x(:,3)./T(1,2))./T(1,2);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 3 - STN: pyr

% pyramidal cells: Hidden causes - error
%--------------------------------------------------------------------------
u      =  0;
u      =  - G(:,4)*S(:,3) + u;
f(:,6) =  (u - 2*x(:,6) - x(:,5)./T(1,3))./T(1,3);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 4 - GPi: ii

% inhibitory interneurons: Hidden states - error
%--------------------------------------------------------------------------
u      =  0; 
u      =  G(:,6)*S(:,5) - G(:,7)*S(:,1) - G(:,8)*S(:,7) + u;
f(:,8) =  (u - 2*x(:,8) - x(:,7)./T(1,4))./T(1,4);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 5 - Thalamus: pyr

% pyramidal cells: Hidden causes - error
%--------------------------------------------------------------------------
u      =  0; 
u      =  - G(:,9)*S(:,7) + u;
f(:,10) =  (u - 2*x(:,10) - x(:,9)./T(1,5))./T(1,5);




% Voltage
%==========================================================================
f(:,1) = x(:,2);
f(:,3) = x(:,4);
f(:,5) = x(:,6);
f(:,7) = x(:,8);
f(:,9) = x(:,10);
f      = spm_vec(f);

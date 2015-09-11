function [F] = mci_nmm_r2p2_dfdp (x,u,P,M)
% Parameter Jacobian for two region, two parameter NMM
% FORMAT [F] = mci_nmm_r2p2_dfdp (x,u,P,M)
%
% x         State
% u         Inputs
% P         Parameters
% M         Model structure
%
% F         F(i,j) = df(x)_i/dp_j
%__________________________________________________________________________
% Copyright (C) 2015 Wellcome Trust Centre for Neuroimaging

% Will Penny and Biswa Sengupta
% $Id: mci_nmm_r2p2_dfdp.m 6548 2015-09-11 12:39:47Z will $

% 18 state variables, 2 parameters
F=zeros(18,2);

curr_P=M.can_P; % Canonical parameter set

w21=P(1);
w12=P(2);

% 2 free parameters
curr_P.A{1}(2,1)=P(1); % Forward connection, w_21
curr_P.A{2}(1,2)=P(2); % Backward connection, w_12

P=curr_P;

% default parameters 
E = [1 1/2 1/8]*32;         % extrinsic rates (forward, backward, lateral)
D = [2 16];                 % propogation delays (intrinsic, extrinsic)
H = [4 32];                 % receptor densities (excitatory, inhibitory)
T = [8 16];                 % synaptic constants (excitatory, inhibitory)
R = [2 1]/3;                % parameters of static nonlinearity

% neuronal states into matrix form; x(r,:) for region r
x = spm_unvec(x,M.x);       
% extrinsic delays
De = D(2).*exp(P.D)/1000;                

% delayed pyramidal cell activity
pyr_ext = presynaptic (x(:,9)-De*(x(:,5)-x(:,6)),P,R);  

Te    = T(1)/1000*exp(P.T(:,1));         % excitatory time constants
He    = H(1)*exp(P.G(:,1));              % excitatory receptor density

HeTe=He/Te;

% Effect of forward connection on region 2 stellate cells
F(8,1)=HeTe*pyr_ext(1)*E(1)*w21*exp(w12); 
% Effect of backward connection on region 1 pyramidal cells
F(9,2)=HeTe*pyr_ext(2)*E(2)*w12*exp(w12); 
% Effect of backward connection on region 1 inhibitory cells
F(15,2)=HeTe*pyr_ext(2)*E(2)*w12*exp(w12); 
end

function [S] = presynaptic (x,P,R)

% pre-synaptic inputs: s(V)
%--------------------------------------------------------------------------
R     = R.*exp(P.S);
S     = 1./(1 + exp(-R(1)*(x - R(2)))) - 1./(1 + exp(R(1)*R(2)));

end

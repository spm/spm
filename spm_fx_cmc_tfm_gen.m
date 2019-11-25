function [ux, vx, wx] = spm_fx_cmc_tfm_gen(x,u,P,M,option)
% Generates pre synaptic signals for multimodal DCM for fMRI and M/EEG
% FORMAT [u,v,w]  = spm_fx_cmc_tfm_gen(x,u,P,M)
% FORMAT [u,v]    = spm_fx_cmc_tfm_gen(x,u,P,M)
% FORMAT [u]      = spm_fx_cmc_tfm_gen(x,u,P,M)
%Input
% -------------------------------------------------------------------------
% x      - state vector
%   x(:,1) - voltage     (spiny stellate cells)
%   x(:,2) - conductance (spiny stellate cells)
%   x(:,3) - voltage     (superficial pyramidal cells)
%   x(:,4) - conductance (superficial pyramidal cells)
%   x(:,5) - current     (inhibitory interneurons)
%   x(:,6) - conductance (inhibitory interneurons)
%   x(:,7) - voltage     (deep pyramidal cells)
%   x(:,8) - conductance (deep pyramidal cells)
% P        - parameters of canonical micro circuits
% u        - exogenous input
% M        - neural-mass model structure
% options  - options for calculation pre synaptic signals:
%   pre (pre synaptic),
%   int (only regional neuronal signals are taken to account for simulating signals),
%    d (different) , s (same) parameters of neurovascular scaling (this
%    option is only here for consistency in the multi model DCM - i.e., it is
%    not used inside this code
%   ext(external neuronal signals are taken into account for simulating presynaptic signals)
%   EX (a 4x1 matrix  with either 0 or 1 elements (order as follows :[ x(:,1)  x(:,3) x(:,5) x(:,7)])
%   to exclude some neuronal signals  from calculation of pre synaptic signals),
%   de(decomposed pre synaptic  activities are taken into account for simulating
%   presynaptic signals)
%
%  Examples of options               {'pre',   'd',    'int', EX},
%                                    {'pre',   's',    'int', EX},
%                                    {'pre',   'd',    'ext', EX},
%                                    {'pre',   's',    'ext', EX},
%                                    {'de',    's',    'int', EX},
%                                    {'de',    's',    'ext', EX},
%                                    {'de',    'd',    'int', EX},
%                                    {'de',    'd',    'ext', EX},
%
%-----------------------------Output---------------------------------------
%  ux = spm_fx_cmc_tfm(x,u,P,M,option)
%  ux  -  simulated presynaptic signal (including or exclude distal regions)
%
% [ux,vx] = spm_fx_cmc_tfm_gen(x,u,P,M,option)
% ux  - intrinsic presynaptic input, (inhibitory)-without external input
% vx  - intrinsic presynaptic input (excitatory)-without external input
%
% [ux,vx,wx] = spm_fx_cmc_tfm_gen(x,u,P,M,,option)
% ux  - intrinsic presynaptic input, (inhibitory)
% vx  - intrinsic presynaptic input (excitatory)
% wx  - extrinsic presynaptic input
%--------------------------------------------------------------------------
% Prior fixed parameter scaling [Defaults]
%
% E  = (forward and backward) extrinsic rates
% G  = intrinsic rates
% D  = propagation delays (intrinsic, extrinsic)
% T  = synaptic time constants
% R  = slope of sigmoid activation function
%
%__________________________________________________________________________
% David O, Friston KJ (2003) A neural mass model for MEG/EEG: coupling and
% neuronal dynamics. NeuroImage 20: 1743-1755
%__________________________________________________________________________
% Copyright (C) 2005 Wellcome Trust Centre for Neuroimaging

% Amirhossein Jafarian, and Karl Friston
% $Id: spm_fx_cmc_tfm_gen.m 7705 2019-11-22 15:06:38Z spm $

persistent in1 in2 in3 ;

% get dimensions and configure state variables
%--------------------------------------------------------------------------
x  = spm_unvec(x,M.x);                % neuronal states
n  = size(x,1);                       % number of sources


% [default] fixed parameters
%--------------------------------------------------------------------------
E  = [2 1 1 1]*512;                   % extrinsic (forward and backward)
T  = [16 8 1 2]*16;                   % synaptic rate constants
R  = 1;                               % gain of activation function
B  = 0;                               % baseline firing

% Extrinsic connections
%--------------------------------------------------------------------------
% ss = spiny stellate
% sp = superficial pyramidal
% dp = deep pyramidal
% ii = inhibitory interneurons
%--------------------------------------------------------------------------
A{1} = exp(P.A{1})*E(1);              % forward  connections (sp -> ss)
A{2} = exp(P.A{2})*E(2);              % forward  connections (sp -> dp)
A{3} = exp(P.A{3})*E(3);              % backward connections (dp -> sp)
A{4} = exp(P.A{4})*E(4);              % backward connections (dp -> ii)

% input connections
%--------------------------------------------------------------------------
C    = exp(P.C);

% pre-synaptic inputs: s(V)
%--------------------------------------------------------------------------
ig   = 2*(1:4);                       % indices of conductance
iv   = ig - 1;                        % indices of voltage
V    = x(:,iv);                       % Voltage
F    = 1./(1 + exp(B - R*V));         % firing rate
S    = F - 1/(1 + exp(B));            % deviation

% input
%==========================================================================
if isfield(M,'u')
    
    % endogenous input
    %----------------------------------------------------------------------
    U   = u(:)*128;
    M.m = size(U,1);
    
else
    % exogenous input
    %----------------------------------------------------------------------
    U   = C*u(:)*2;
    M.m = size(C,2);
    
end
clear u


% time constants and intrinsic connections
%==========================================================================
T      = ones(n,1)*T;
i      = 1:size(P.T,2);
T(:,i) = T(:,i).*exp(P.T);

% extrinsic connections
%--------------------------------------------------------------------------
% forward  A{1}  sp -> ss
% forward  A{2}  sp -> dp
% backward A{3}  dp -> sp
% backward A{4}  dp -> ii
%--------------------------------------------------------------------------

% Neuronal states (deviations from baseline)
%--------------------------------------------------------------------------
%   x(:,1) - voltage     (spiny stellate cells)
%   x(:,2) - conductance (spiny stellate cells)but
%   x(:,3) - voltage     (superficial pyramidal cells)
%   x(:,4) - conductance (superficial pyramidal cells)
%   x(:,5) - current     (inhibitory interneurons)
%   x(:,6) - conductance (inhibitory interneurons)
%   x(:,7) - voltage     (deep pyramidal cells)
%   x(:,8) - conductance (deep pyramidal cells)
%--------------------------------------------------------------------------
%     ss sp ii dp   % intrinsic connections
%--------------------------------------------------------------------------
g  = [-8 -4  -4  0;  % ss
    4 -8  -2  0;  % sp
    4  2  -4 -2;  % ii
    0  1  -2 -4]; % dp

g  = g*256*exp(P.S);

% intrinsic connections to be optimised (only the first is modulated)
%--------------------------------------------------------------------------
G       = ones(n,1)*diag(g)';
i       = 1:size(P.G,2);
G(:,i)  = G(:,i).*exp(P.G);


% Modulatory effects of sp depolarisation on recurrent inhibition
%--------------------------------------------------------------------------
if isfield(P,'M')
    G(:,2) = G(:,2).*exp(-P.M*32*S(:,2));
end


% Motion of states: f(x)
%--------------------------------------------------------------------------
% Afferents
%==========================================================================
u      = zeros(n,4);                % intrinsic - inhibitory
v      = zeros(n,4);                % intrinsic - excitatory
w      = zeros(n,4);                % extrinsic - excitatory
Pre_ext    = zeros(n,4); % pre-synaptics with extrinsic input;
Pre_no_ext = zeros(n,4); % pre-synaptics without extrinsic input;

% Granular layer (excitatory interneurons): spiny stellate: Hidden causes
%--------------------------------------------------------------------------
u(:,1) = G(:,1).*S(:,1) + g(1,3)*S(:,3) + g(1,2)*S(:,2);
w(:,1) = A{1}*S(:,2) + U;
Pre_ext(:,1)    = u(:,1);
Pre_no_ext(:,1) = u(:,1) + w(:,1);

if option{1,4}(1)== 0
    u(:,1) = 0;
    w(:,1) = 0;
    Pre_ext   (:,1)    = 0;
    Pre_no_ext(:,1)    = 0;
end

% Supra-granular layer (superficial pyramidal cells): Hidden causes - error
%--------------------------------------------------------------------------
u(:,2) = G(:,2).*S(:,2) + g(2,3)*S(:,3);
v(:,2) = g(2,1)*S(:,1);
w(:,2) = A{3}*S(:,4);
Pre_ext(:,2)    = u(:,2) + v(:,2);
Pre_no_ext(:,2) = u(:,2) + v(:,2) + w(:,2);

if option{1,4}(2)== 0
    u(:,2) = 0;
    v(:,2) = 0;
    w(:,2) = 0;
    Pre_ext(:,2)    = 0;
    Pre_no_ext(:,2) = 0;
end

% Supra-granular layer (inhibitory interneurons): Hidden states - error
%--------------------------------------------------------------------------
u(:,3) = G(:,3).*S(:,3);
v(:,3) = g(3,1)*S(:,1) + g(3,4)*S(:,4) + g(3,2)*S(:,2);
w(:,3) = A{4}*S(:,4);
Pre_ext(:,3)    = u(:,3) + v(:,3);
Pre_no_ext(:,3) = u(:,3) + v(:,3) + w(:,3);
if option{1,4}(3)==0
    u(:,3) = 0;
    v(:,3) = 0;
    w(:,3) = 0;
    Pre_ext(:,3)    = 0;
    Pre_no_ext(:,3) = 0;
end

% Infra-granular layer (deep pyramidal cells): Hidden states
%--------------------------------------------------------------------------
u(:,4) = G(:,4).*S(:,4) + g(4,3)*S(:,3);
v(:,4) = g(4,2)*S(:,2);
w(:,4) = A{2}*S(:,2);
Pre_ext(:,4)    = u(:,4) + v(:,4);
Pre_no_ext(:,4) = u(:,4) + v(:,4) + w(:,4);
if option{1,4}(4)==0
    u(:,4) = 0;
    v(:,4) = 0;
    w(:,4) = 0;
    Pre_ext(:,4)    = 0;
    Pre_no_ext(:,4) = 0;
end

% output
%--------------------------------------------------------------------------
t = 1 ;
if (strcmp(option{3}, 'ext')&& strcmp(option{1}, 'pre'))
    in1(:,:,t)  = Pre_ext;
    ux(:,:,t)    = in1;
end
if(strcmp(option{3}, 'int') && strcmp(option{1}, 'pre'))
    in1(:,:,t)  = Pre_no_ext;
    ux(:,:,t)    = in1;
end
if(strcmp(option{1}, 'de'))
    in1(:,:,t)  = u ;
    in2(:,:,t)  = v;
    in3(:,:,t)  = w;
    ux(:,:,t)  = in1;
    vx(:,:,t)  = in2;
    wx(:,:,t)  = in3;
end

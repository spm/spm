function [E,V] = spm_wendling_priors(A,B,C)
% Prior moments for a Wendling neural mass model of hippocampal LFP
% FORMAT [pE,pC] = spm_wendling_priors(A,B,C)
%
% A{3},B{m},C  - binary constraints on extrinsic connectivity
%
% pE - prior expectation
%
% synaptic parameters
%--------------------------------------------------------------------------
%    pE.E - excitatory synaptic gain (EXC, scales A = 3.25 mV)
%    pE.S - slow dendritic inhibitory gain (SDI, scales B = 22 mV)
%    pE.F - fast somatic inhibitory gain (FSI, scales G = 10 mV)
%    pE.T - synaptic time constants [excitatory, slow inh, fast inh]
%    pE.H - intrinsic connectivity C1-C7
%    pE.R - activation function parameters [slope, threshold]
%
% connectivity parameters
%--------------------------------------------------------------------------
%    pE.A - extrinsic coupling
%    pE.B - trial-dependent
%    pE.C - extrinsic stimulus input
%    pE.D - extrinsic delays
%    pE.I - intrinsic delays
%
% pC - prior (co)variances
%
% Intrinsic connectivity C1-C7 (Wendling convention):
%   C1: Pyramidal -> Excitatory interneurons
%   C2: Excitatory interneurons -> Pyramidal
%   C3: Pyramidal -> Slow dendritic inhibitory
%   C4: Slow dendritic inhibitory -> Pyramidal
%   C5: Pyramidal -> Fast somatic inhibitory
%   C6: Fast somatic inhibitory -> Pyramidal
%   C7: Slow inhibitory -> Fast inhibitory
%__________________________________________________________________________
%
% Wendling F, Bartolomei F, Bellanger JJ, Chauvel P (2002)
%__________________________________________________________________________

% defaults
%--------------------------------------------------------------------------
if nargin < 3
    A = {0 0 0};
    B = {0};
    C = 1;
end
n   = size(C,1);                                    % number of sources

% disable log zero warning
%--------------------------------------------------------------------------
warning('off','MATLAB:log:logOfZero');

% parameters for Wendling neural mass model
%==========================================================================

% sigmoid parameters
%--------------------------------------------------------------------------
E.R = [0 0];               V.R = [1 1]/8;

% synaptic gains (log-deviations from default)
%--------------------------------------------------------------------------
E.E = zeros(n,1);          V.E = ones(n,1)/8;      % excitatory gain (EXC)
E.S = zeros(n,1);          V.S = ones(n,1)/8;      % slow inhibitory gain (SDI)
E.F = zeros(n,1);          V.F = ones(n,1)/8;      % fast inhibitory gain (FSI)

% time constants [excitatory, slow inhibitory, fast inhibitory]
%--------------------------------------------------------------------------
E.T = zeros(n,3);          V.T = ones(n,3)/8;      % time constants

% intrinsic connectivity C1-C7
%--------------------------------------------------------------------------
E.H = zeros(n,7);
V.H = ones(n,7)/16;                                % C1-C7: all free

% extrinsic connectivity
%--------------------------------------------------------------------------
Q     = sparse(n,n);
for i = 1:length(A)
      A{i} = ~~A{i};
    E.A{i} = log(A{i} + eps);                      % forward
    V.A{i} = A{i}/2;                               % backward
    Q      = Q | A{i};                             % and lateral connections
end

for i = 1:length(B)
      B{i} = ~~B{i};
    E.B{i} = 0*B{i};                               % input-dependent scaling
    V.B{i} = B{i}/2;
    Q      = Q | B{i};
end
C      = ~~C;
E.C    = C*32 - 32;                                % where inputs enter
V.C    = C/32;

% delays
%--------------------------------------------------------------------------
E.D    = sparse(n,n);     V.D = Q/16;              % extrinsic delays
E.I    = 0;               V.I = 1/32;              % intrinsic delays

warning('on','MATLAB:log:logOfZero');

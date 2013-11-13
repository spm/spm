function [pE,pC] = spm_cmm_NMDA_priors(A,B,C)
% prior moments for a canonical neural-mass model of ERPs
% FORMAT [pE,pC] = spm_cmm_priors(A,B,C)
%
% A{3},B{m},C  - binary constraints on extrinsic connections
%
% pE - prior expectation - f(x,u,P,M)
%
% population variance
%--------------------------------------------------------------------------
%     E.S    - variance
%
% synaptic parameters
%--------------------------------------------------------------------------
%    pE.T    - synaptic time constants
%    pE.G    - intrinsic connectivity
%
% connectivity parameters
%--------------------------------------------------------------------------
%    pE.A    - extrinsic
%    pE.B    - trial-dependent
%    pE.C    - stimulus input
%
% stimulus and noise parameters
%--------------------------------------------------------------------------
%    pE.R    - onset and dispersion
%    pE.D    - delays
%    pE.U    - exogenous background activity
%
% pC - prior (co)variances
%
% Because priors are specified under log normal assumptions, most
% parameters are simply scaling coefficients with a prior expectation
% and variance of one.  After log transform this renders pE = 0 and
% pC = 1;  The prior expectations of what they scale are specified in
% spm_erp_fx
%__________________________________________________________________________
%
% David O, Friston KJ (2003) A neural mass model for MEG/EEG: coupling and
% neuronal dynamics. NeuroImage 20: 1743-1755
%__________________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging
 
% Karl Friston
% $Id: spm_cmm_NMDA_priors.m 5741 2013-11-13 12:10:48Z guillaume $
 
 
% disable log zero warning
%--------------------------------------------------------------------------
warning('off','MATLAB:log:logOfZero');
n     = size(C,1);                                % number of sources
u     = size(C,2);                                % number of inputs
p     = 4;                                        % number of populations
 
% parameters for neural-mass forward model
%==========================================================================
 
% population variance (NMDA mediates the effects of BN(1,1)
%--------------------------------------------------------------------------
pE.S      = 0;             pC.S    = 1/64;

 
% intrinic [excitatory AMPA, inhibitory GABAa, and excitatory NMDA ] time constants (H the effects of B)
%--------------------------------------------------------------------------
pE.T  = zeros(n,3);        pC.T =  ones(n,3)/64;
pE.G  = zeros(n,1);        pC.G = zeros(n,1)/16;

% Capacitance and backround activity
%--------------------------------------------------------------------------
pE.CV = zeros(1,p);     pC.CV = ones(1,p)/16;
pE.E  = 0;              pC.E  = 1/64;


% extrinsic connectivity (n x n)
%==========================================================================

% restructure adjacency matrices
%--------------------------------------------------------------------------
A{1}  = A{1} | A{3};                              % forward
A{2}  = A{2} | A{3};                              % backward
for i = 1:2
    pE.A{i} = A{i}*32 - 32;
    pC.A{i} = A{i}/8;
end

for i = 1:2
    pE.AN{i} = A{i}*32 - 32;
    pC.AN{i} = A{i}/8;
end
% input-dependent scaling
%--------------------------------------------------------------------------
for i = 1:length(B)
    B{i} = ~~B{i};
    pE.B{i} = B{i} - B{i};
    pC.B{i} = B{i}/8;
end


for i = 1:length(B)
    B{i} = ~~B{i};
    pE.BN{i} = B{i} - B{i};
    pC.BN{i} = B{i}/8;
end

% exogenous inputs
%--------------------------------------------------------------------------
C     = ~~C;
pE.C  = C*32 - 32;
pC.C  = C/8;



% intrinsic connectivity (p x p x n)
%==========================================================================
% 1 - excitatory spiny stellate cells (granular input cells)
% 2 - superficial pyramidal cells     (forward  output cells)
% 3 - inhibitory interneurons         (intrisic interneuons)
% 4 - deep pyramidal cells            (backward output cells)
%--------------------------------------------------------------------------
gC    = [1   0   1   0;
         1   1   1   0;
         1   0   1   1;
         0   1   1   1]/32;
     
pE.H  = repmat(zeros(p,p),[1 1 n]);
pC.H  = repmat(gC        ,[1 1 n]);

% Exogenous inputs: onset and dispersion
%--------------------------------------------------------------------------
pE.R  = zeros(u,2);     pC.R  = ones(u,1)*[1/16 1/16];

% Delays (intrinsic and extrinsic)
%--------------------------------------------------------------------------
pE.D  = [0 0];          pC.D  = [1 1]/64;


warning('on','MATLAB:log:logOfZero');

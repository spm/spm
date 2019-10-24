function DEM = DEM_evidence_accumulation
% Saccadic eye movements under active inference:
%__________________________________________________________________________
% This demo illustrates evidence accumulation (and responses) using a very
% simple generative model. In this model, there are three hidden states
% corresponding to right motion, no motion and left motion - as registered
% uniformly over 16 visual channels. Motion is slowly introduced, which
% moves the hidden states to one of the unstable fixed points; thereby
% inducing proprioceptive predictions that cause a motor response. The
% generative model is as minimal as possible and is based on generalised
% Lotka-Volterra dynamics to emulate a dynamic form of winner takes all. In
% other words, the only prior beliefs of this generative model are that the
% world can be in one of a number of (unstable) states. Evidence is
% accumulated slowly because the input is noisy (and is assigned a low
% precision). This reveals the evidence accumulation dynamics that drive
% action, when inference is sufficiently confident. These dynamics are
% formally equivalent to the race or drift diffusion dynamics in normative
% (descriptive) formulations of evidence accumulation.
%__________________________________________________________________________
% Copyright (C) 2011 Wellcome Trust Centre for Neuroimaging
 
% Karl Friston
% $Id: DEM_evidence_accumulation.m 7679 2019-10-24 15:54:07Z spm $
 
 
% hidden causes and states
%==========================================================================
% x    - hidden states:
%   x(1) - log likelihood of hypothesis 1 (e.g. moving right)
%   x(2) - log likelihood of hypothesis 2 (e.g. not moving)
%   x(3) - log likelihood of hypothesis 3 (e.g. moving left)
%
% v    - hidden causes (none for model and moving left for process)
%
% g    - sensations:
%   g(1) - looking to the left (proprioception)
%   g(3) - retinal input - channel 1
%   g(4) - retinal input - channel 2
%   g(5) - ...
%--------------------------------------------------------------------------
 
 
% generative model
%==========================================================================
M(1).E.s = 1/2;                                % smoothness
M(1).E.n = 3;                                  % order of
M(1).E.d = 1;                                  % generalised motion
 
% level 1: mmulti-stable dynamics generating signal over all channels
%--------------------------------------------------------------------------
M(1).f  = @(x,v,P)(1 - sum(exp(x)) - (x + log(3))/32)/64;
M(1).g  = @(x,v,P)[[-1 0 1]*spm_phi(16*(exp(x) - 1)); [-1 0 1]*exp(x) + zeros(16,1)];
M(1).x  = -log(3)*[1; 1; 1];                   % hidden states
M(1).V  = [0 (zeros(1,16) + exp(-4))];         % error precision (g)
M(1).W  = exp(8);                              % error precision (f)
 
 
% generative process
%==========================================================================
 
% first level
%--------------------------------------------------------------------------
G(1).g  = @(x,v,a,P)[a; v + zeros(16,1)]';
G(1).V  = [exp(8) zeros(1,16) + exp(2)];       % error precision
G(1).U  = [exp(0) zeros(1,16)];                % motor gain
 
% second level
%--------------------------------------------------------------------------
G(2).v  = 0;                                   % exogenous forces
G(2).a  = 0;                                   % action forces
G(2).V  = exp(16);
 
 
% generate and invert
%==========================================================================
N     = 32;                                    % length of data sequence
pst   = (1:N);                                 % peristimulus time (bins)
 
DEM.G = G;
DEM.M = M;
DEM.C = spm_phi(((1:N) - N/4)*32/N);           % stimulus motion
 
DEM   = spm_ADEM(DEM);
spm_DEM_qU(DEM.qU,DEM.pU)
 
% plot exponential transform of posterior estimates of hidden states
%--------------------------------------------------------------------------
subplot(2,2,2)
spm_plot_ci(DEM.qU.x{1},DEM.qU.S,pst,(1:3),'exp'), hold on
axis([1 N 0 2]), axis square

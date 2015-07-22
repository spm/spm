function DEM_demo_ontology
% Demo for a bird songs: In this example, we simulate local field potential
% using the prediction error from the song-bird example below. We look at
% these responses under natural stimuli and after removing the second
% level of the hierarchy to show it is necessary for veridical perception.
% We then repeat but omitting dynamical priors by forsaking generalised 
% coordinates
%__________________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging
 
% Karl Friston
% $Id: DEM_demo_ontology.m 6502 2015-07-22 11:37:13Z karl $
 
 
% hierarchical non-linear generative model (dynamic & chaotic)
%==========================================================================

rng('default')

% timing
%--------------------------------------------------------------------------
N        = 128;                      % length of stimulus (bins)
 
% correlations
%--------------------------------------------------------------------------
M(1).E.s = 1;
M(1).E.K = exp(-2);


% level 1
%--------------------------------------------------------------------------
% P(1): Prandtl number
% P(2): 8/3
% P(3): Rayleigh number
 
P.A      = rand(4,2)/128;
P.B      = randn(3,2)/8;
x        = [0.9; 0.8];
M(1).f   = @(x,v,P) (v - x*8)/128;
M(1).g   = @(x,v,P)[spm_phi(P.A*exp(x)); spm_softmax(P.B*x)];
M(1).x   = x;
M(1).pE  = P;
M(1).V   = exp(1);
M(1).W   = exp(8);
 
 
% level 2
%--------------------------------------------------------------------------
P        = [10; 8/3];
x        = [0.9; 0.8; 30];
M(2).f   = @(x,v,P)[-P(1) P(1) 0; (32 - x(3)) -1 0; x(2) 0 -P(2)]*x/128;
M(2).g   = @(x,v,P) x([2 3]);
M(2).x   = x;
M(2).pE  = P;
M(2).V   = exp(8);
M(2).W   = exp(8);
 
 
% create data
%==========================================================================
 
% create innovations & add causes
%--------------------------------------------------------------------------
DEM   = spm_DEM_generate(M,N);
 

% deconvolve
%--------------------------------------------------------------------------
DEM   = spm_DEM(DEM);



% show songs and prediction error (ERP)
%==========================================================================
spm_DEM_qU(DEM.qU,DEM.pU)
colormap('pink')





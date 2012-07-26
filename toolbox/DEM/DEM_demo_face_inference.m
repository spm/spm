function DEM_demo_face_inference
% Recognising facial expressions: This demo uses the linear convolution
% model with two hidden states and one cause to generate coefficients of
% visual basis functions that produce a moving face. The basis functions are
% images have been chosen so that the appropriate nonlinear mixture creates
% a smile. The coefficients of the i-th basis image is 
%
% cos((i - 1)*pi*g(x))
%
% where g(x) is some none linear mixture of hidden sates that lies in the
% range [0,1]. (neutral to smiling). Inversion of this model corresponds to
% nonlinear Bayesian de-convolution of visual input to recognise the dynamic
% expressions. The associated (roving MMN) demonstration uses this
% generative model to illustrate perceptual learning and repetition suppression
% when we repeat the stimulus.  Clicking on the images will display the
% movies entailed by the true and estimated causes.
%__________________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging
 
% Karl Friston
% $Id: DEM_demo_face_inference.m 4804 2012-07-26 13:14:18Z karl $


% temporal smoothness - s.d. of kernel
%--------------------------------------------------------------------------
G(1).E.s = 1;
 
% level 1
%--------------------------------------------------------------------------
G(1).m  = 1;                                % 1 input or cause
G(1).n  = 2;                                % 2 hidden states
G(1).l  = 3;                                % 3 outputs (coefficients for face)

P       = struct;
P.f     = [-1  4 ;
           -2 -1]/16;                        % The Jacobian
P.g     = [0 -1];                           % The mixing parameters
 
G(1).f  = inline('P.f*x + [v; 0]','x','v','P');
G(1).g  = inline('cos([0:2]*pi*spm_phi(P.g*x))','x','v','P');
G(1).pE = P;                                % The prior expectation
G(1).V  = exp(8);                           % error precision (data)
G(1).W  = exp(16);                          % error precision (states)                                          % with a low level of noise
 
% level 2
%--------------------------------------------------------------------------
G(2).l  = 1;                                % 1 output 
G(2).V  = exp(16);                          % error precision (cause)
 
% The data: Data [stimuli] are created by integrating the model for some
% input.  The input here is simply a bump [Gaussian] function.
%==========================================================================
M       = G;                                % make M the canonical model
N       = 64;                               % length of data sequence
 
% create innovations & add causes
%--------------------------------------------------------------------------
c       = exp(-((1:N) - 16).^2/(2.^2));     % this is the Gaussian cause
 
% integrate G to obtain causal (v) and hidden states (x)
%--------------------------------------------------------------------------
DEM     = spm_DEM_generate(G,c,P);
 
 
% invert
%==========================================================================
DEM.M(1).V = exp(8);
DEM.M(1).W = exp(8);
DEM.M(2).V = exp(2);
DEM        = spm_DEM(DEM);
 
% render true and perceived stimuli in move format
%--------------------------------------------------------------------------

% plot causal and hidden states
%--------------------------------------------------------------------------
spm_DEM_qU(DEM.qU,DEM.pU)

% render true and perceived stimuli in move format
%--------------------------------------------------------------------------
subplot(4,4,11)
spm_DEM_movie(DEM.pU);
title('true','FontSize',16)
 
subplot(4,4,12)
spm_DEM_movie(DEM.qU);
title('perceived','FontSize',16)

subplot(4,4,15)
spm_DEM_movie([1 0 0 0]');
title('basis 1','FontSize',16)
 
subplot(4,4,16)
spm_DEM_movie([0 1 0 0]');
title('basis 2','FontSize',16)

disp(' '),disp('Click faces to play movies'),disp(' ')

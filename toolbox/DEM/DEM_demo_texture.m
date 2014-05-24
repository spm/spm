function LAP = DEM_demo_texture
% This demonstration considers the figure-ground segregation problem where, 
% crucially, a figure is defined texturally – in terms of its second order 
% statistics; in other words, a visual object is manifest in terms of its 
% texture or spectral power density in the spatial domain. This definition 
% precludes any first-order attributes; such as increased luminance. This 
% sort of problem is common in the inverse literature and is usually solved 
% using a prior on the [co]variance of random fluctuations generating data. 
% Here, we simulate a contiguous object, whose texture is determined by the 
% variance of random fluctuations in luminance – and the variance (or 
% precision) is modulated by Gaussian basis functions. The resulting signal 
% is mixed with uniform Gaussian noise to produce sensory data. These 
% (one-dimensional) data are then subject to Bayesian inversion using 
% generalized predictive coding – (as implemented in spm_LAP) – to recover 
% the underlying object.
% 
% Technically, this scheme optimizes expectations of the hidden causes of 
% the data, which are the amplitude of radial basis functions controlling 
% the precision of retinotopic signals. Heuristically, the figure is 
% recognized by selectively attending to sensory input from the figure. 
% This enables sensory noise to be suppressed in unattended parts of the 
% visual field. However, this form of attention is distinct from simply 
% boosting sensory precision (the precision of sensory prediction errors) 
% as in simulations of the Posner paradigm or biased competition. Here, 
% the hidden causes are optimized in a way that renders them less precise 
% and therefore more sensitive to ascending (prediction error) sensory 
% input. This illustrates the functional importance of the relative 
% precision of sensory and extrasensory prediction errors in modulating 
% the influence of ascending sensory information that competes to influence 
% posterior expectations.
%
%__________________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging
 
% Karl Friston
% $Id: DEM_demo_texture.m 6017 2014-05-24 14:36:02Z karl $
 

% Create a generative model:
%==========================================================================
rng('default')

% tell the model precision depends on hidden causes
%--------------------------------------------------------------------------
G(1).E.method.v = 1;
G(1).E.linear   = 1;
                                       
% level 1; textured stimulus with noise (log precision of four)
%--------------------------------------------------------------------------
G(1).v  = zeros(128,1);                   % output channels (stimuli)
G(1).V  = exp(4);                         % error variances (noise)
G(1).g  = inline('v','x','v','P');


% level 2; underlying causes (three Gaussian patches)
%--------------------------------------------------------------------------
ph      = '6 - exp(-((1:128)''*ones(1,3) - ones(128,1)*[48 64 80]).^2/32)*v(1:3)';
G(2).v  = zeros(128,1);                   % textured stimulus
G(2).V  = [];
G(2).ph = inline(ph,'x','v','h','M');
G(2).g  = inline('zeros(128,1)','x','v','P');


% level 2; amplitude and size
%--------------------------------------------------------------------------
U       = [8;8;0];                        % amplitude and size
G(3).v  = U;
G(3).V  = 1/2;


% evaluate G to generate stimulus
%--------------------------------------------------------------------------
LAP   = spm_DEM_generate(G,U);
 
% invert to simulate  predictive coding
%==========================================================================
LAP   = spm_LAP(LAP);

% plot results
%--------------------------------------------------------------------------
spm_DEM_qU(LAP.qU,LAP.pU)

f  = 6 - G(2).ph([],U,[],G);

subplot(3,2,1), title('prediction and error','FontSize',16)
subplot(3,2,2), title('signal plus noise','FontSize',16)
subplot(3,2,3), title('prediction','FontSize',16)
subplot(3,2,4), title('true signal','FontSize',16)
subplot(3,2,6), plot(U), spm_axis tight, axis square, box off
title('true causes','FontSize',16)
subplot(3,2,4), hold on, plot(f/4,':r')


 

function DCM = spm_meta_model(DCM)
% Meta-modelling of Bayes-optimal responses (Newton's method)
% FORMAT DCM = spm_meta_model(DCM)
%
% store estimates in DCM
%--------------------------------------------------------------------------
% DCM.M   - meta-model specification
%      M: [1 x m struct]   - hierarchical inference model (cf DEM.M)
%      G: [1 x s struct]   - generative process (for spm_ADEM or spm_ALAP)
%      U: [n x N double]   - n prior beliefs over N samples
%     pE: [1 x 1 struct]   - prior expectation of meta-model parameters
%     pC: [1 x 1 struct]   - prior variance of meta-model parameters
%
% DCM.xY  - data structure
%      y: [N x p double]   - N samples of a p-variate response
%     X0: [N x q double]   - q-variate confounds
%     dt: [1 x 1 double]   - size of time bin for each sample
%      Q: {[N x N double]} - precision component[s]
%
% DCM.xU  - input structure
%      u: [r x N double]   - r-variate input (hidden causes G in DEM)
%
% Computes (and stores in DCM_MM_???)
%--------------------------------------------------------------------------
% DCM.DEM - Inference (with MAP parameters)
% DCM.Ep  - conditional expectation
% DCM.Cp  - conditional covariances
% DCM.Eh  - conditional log-precision
% DCM.Ey  - conditional response
% DCM.F   - log-evidence
%
% This routine illustrates Bayesian meta modelling - the Bayesian inversion
% of a model of a Bayesian observer. This requires the specification of two
% models: An inference model used by the subject (specified by a DEM
% structure) and a meta-model (specified by a DCM structure). The inference
% model is completed by a response model to furnish the meta-model; where
% the response model takes the output of the (active) inference scheme
% specified by the DEM and generates an observed (behavioural or
% neurophysiological) response. Crucially either the inference model or
% the response model or both can have free parameters - that are optimised
% using Bayesian nonlinear system identification in the usual way.
%
% Although this routine is a function, it is expected that people will fill
% in the model-specific parts in a local copy, before running it.  The
% current example uses a model of slow pursuit and generates synthetic data
% (responses) to illustrate how it works. To replace these simulated data
% with real data, simply specify the DCM.xY (and xU fields) with
% empirical values. If other fields do not exist, exemplar fields will be filled in.
%
% The conditional density of the parameters and F values (log-evidence) can
% be used in the usual way for inference on parameters or Bayesian model
% comparison (as for other DCMs)
%__________________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging
 
% Karl Friston
% $Id: spm_meta_model.m 5788 2013-12-06 20:08:57Z karl $
 
 
 
% INFERENCE MODEL: DEM.M and DEM.G
%==========================================================================
try
    M = DCM.M;
    
catch
    
    % create illustrative (M,G) if not specified (EXAMPLE)
    %----------------------------------------------------------------------

    % hidden causes and states
    %======================================================================
    % x    - hidden states:
    %   x.o(1) - oculomotor angle
    %   x.o(2) - oculomotor velocity
    %   x.x(1) - target angle - extrinsic coordinates
    %
    % v    - causal states: force on target
    %
    % g    - sensations:
    %   g(1) - oculomotor angle (proprioception)
    %   g(2) - oculomotor velocity
    %   g(:) - visual input - intrinsic coordinates
    %----------------------------------------------------------------------
    
    
    % Set-up
    %======================================================================
    M(1).E.s = 1/2;                               % smoothness
    M(1).E.n = 4;                                 % order of
    M(1).E.d = 1;                                 % generalised motion
    
    % angular frequency of target motion (per time bin)
    %----------------------------------------------------------------------
    w  = 2*pi/32;
    
    
    % occlusion and visual mapping functions
    %----------------------------------------------------------------------
    occ = @(x) x < 1/2;
    vis = @(x) exp(-((-8:8)' - x.x(1) + x.o(1)).^2)*occ(x.x);
    
    
    % Generative model (M) (sinusoidal movement)
    %======================================================================
    % Endow the model with internal dynamics (a simple oscillator) so that is
    % recognises and remembers the trajectory to anticipate jumps in rectified
    % sinusoidal motion.
    
    % slow pursuit following with (second order) generative model
    %----------------------------------------------------------------------
    x.o = [0;0];                                  % motor angle & velocity
    x.x = 0;                                      % target location
    
    % level 1: Displacement dynamics and mapping to sensory/proprioception
    %----------------------------------------------------------------------
    M(1).f = @(x,v,P) [x.o(2); (v - x.o(1))/4 - x.o(2)/2; v - x.x];
    M(1).g = @(x,v,P) [x.o; vis(x)];
    M(1).x = x;                                   % hidden states
    M(1).V = exp(4);                              % error precision
    M(1).W = exp(4);                              % error precision
    
    
    % level 2: With hidden (memory) states
    %----------------------------------------------------------------------
    M(2).f  = @(x,v,P)[x(2); -x(1)]*v;
    M(2).g  = @(x,v,P) x(1);
    M(2).x  = [0; 0];                             % hidden states
    M(2).V  = exp(4);                             % error precision
    M(2).W  = exp(4);                             % error precision
    
    % level 3: Encoding frequency of memory states (U)
    %----------------------------------------------------------------------
    M(3).v = 0;
    M(3).V = exp(4);
    
    
    % generative model (G) (with no noise)
    %======================================================================
    
    % first level
    %----------------------------------------------------------------------
    G(1).f = @(x,v,a,P) [x.o(2); a/4 - x.o(2)/8; v - x.x];
    G(1).g = @(x,v,a,P) [x.o; vis(x)];
    G(1).x = x;                                  % hidden states
    G(1).U = sparse(1,[1 2],[1 1],1,19)*exp(4);  % motor gain
    
    % second level
    %----------------------------------------------------------------------
    G(2).v = 0;                                  % exogenous force
    G(2).a = 0;                                  % action force

end
 
% Check generative model
%--------------------------------------------------------------------------
M      = spm_DEM_M_set(M);
 
 
% Experimental input (C)
%==========================================================================
try
    xU    = DCM.xU;
    N     = size(xU.u,2);
   
catch
    
    % Sine wave cause (EXAMPLE)
    %----------------------------------------------------------------------
    N     = 64;                                 % length of data sequence
    xU.dt = 16;                                 % time step (ms)
    xU.u  = sin((1:N)*w).*((1:N) > 16);         % sinusoidal target motion
end
 
 
% Priors (U)
%==========================================================================
try
    U = DCM.U;

catch
    
    % (EXAMPLE) prior beliefs about target frequency (w)
    %----------------------------------------------------------------------
    U = zeros(1,N) + w;                   
end
 
 
% Data and confounds
%==========================================================================
try
    xY = DCM.xY;
    
catch
    
    % Generate simulated data (EXAMPLE)
    %----------------------------------------------------------------------
    
    % meta-model and TRUE parameters
    %----------------------------------------------------------------------
    MM.M    = M;
    MM.G    = G;
    MM.U    = U;
    
    P.G     = 2;
    P.W     = 4;
    [y DEM] = IS(P,MM,xU);
    
    
    % and simulate observation noise
    %----------------------------------------------------------------------
    xY.y    = y + spm_Q(1/2,N,1)*randn(N,1)/16;
    
    spm_figure('GetWin','Simulated respsone');
    spm_DEM_qU(DEM.qU,DEM.pU)
 
end
 
% DCT confounds
%==========================================================================
xY.X0  = spm_dctmtx(N,1);
 
% and [serial] correlations (precision components) AR model
%--------------------------------------------------------------------------
xY.Q   = {spm_Q(1/2,N,1)};
 
 
% META-MODEL parameters
%==========================================================================
% These parameterise the inference scheme function (IS) at the bottom of
% this script - this functions defines how the parameters are used and
% therefore optimised. Change this function to specify the meta-model
 
try
    pE    = DCM.pE;
    pC    = DCM.pC;
    
catch
    
    % prior expectations (EXAMPLE)
    %----------------------------------------------------------------------
    pE.G  = 1;
    pE.W  = 2;
    
    % prior covariance (EXAMPLE)
    %----------------------------------------------------------------------
    pC.G  = 1;
    pC.W  = 1;
    
end


% This completes the meta-model specification. The fields are now assembled 
% and passed to spm_nlsi_GN (nonlinear system identification using Gauss-
% Newton-like gradient ascent).
%==========================================================================

 
% META-MODEL (MM)
%==========================================================================
MM.M   = M;
MM.G   = G;
MM.U   = U;
 
% hyperpriors (assuming about 99% signal to noise)
%--------------------------------------------------------------------------
hE     = 6 - log(var(spm_vec(xY.y)));
hC     = exp(-4);
 
% Meta-model
%--------------------------------------------------------------------------
MM.IS  = @IS;
MM.pE  = pE;
MM.pC  = pC;
MM.hE  = hE;
MM.hC  = hC;
 
 
% model inversion
%==========================================================================
MM.Nmax      = 16;
[Ep,Cp,Eh,F] = spm_nlsi_GN(MM,xU,xY);
 
 
% integrate (A)DEM scheme with MAP etimates
%--------------------------------------------------------------------------
[Ey,DEM] = IS(Ep,MM,xU);
 
 
% store estimates in DCM
%--------------------------------------------------------------------------
DCM.DEM = DEM;                  % Inference scheme (with MAP parameters)
DCM.M   = MM;                   % meta-model
DCM.xY  = xY;                   % data structure
DCM.xU  = xU;                   % input structure
DCM.Ep  = Ep;                   % conditional expectation
DCM.Cp  = Cp;                   % conditional covariances
DCM.Eh  = Eh;                   % conditional log-precision
DCM.Ey  = Ey;                   % conditional response
DCM.F   = F;                    % log-evidence
 
 
% save
%--------------------------------------------------------------------------
save(sprintf('DCM_MM_%s',date),'DCM', spm_get_defaults('mat.format'));


% and plot
%--------------------------------------------------------------------------
spm_figure('GetWin','Parameter estimates');
subplot(2,1,1)
spm_plot_ci(Ep,Cp), hold on
Pp       = spm_vec(P);
pE       = spm_vec(pE);
plot(Pp(1),Pp(2),'.r' ,'MarkerSize',32), hold on
plot(pE(1),pE(2),'.g','MarkerSize',32), hold off
axis([0 4 0 8])
xlabel('Gain','FontSize',12)
ylabel('Precision','FontSize',12)
title('Posterior (blue), prior (green) and true (red) value','FontSize',16)

 
return
 
 
% inference scheme x = IS(p,M,U)
%==========================================================================
function [y,DEM] = IS(P,M,U)
 
% parameterise inference model (EXAMPLE)
%--------------------------------------------------------------------------
M.M(2).W = exp(P.W);
 
% reset random number generator (for action schemes)
%--------------------------------------------------------------------------
rng('default')
 
 
% Bayesian inference (spm_DEM, spm_LAP, spm_ADEM or spm_ALAP)
%--------------------------------------------------------------------------
DEM.M  = M.M;
DEM.G  = M.G;
DEM.U  = M.U;
DEM.C  = U.u;
DEM    = spm_ADEM(DEM);
 
% (Hidden) behavioural or electrophysiological response (EXAMPLE)
%--------------------------------------------------------------------------
x      = DEM.pU.x{1}(1,:);
y      = (P.G*x)';
 
return
 
 
 



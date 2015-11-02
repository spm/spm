function DEM_evidence_accumulation
% Saccadic eye movements under active inference:
%__________________________________________________________________________
% This demo illustrates exploration or visual search in terms of optimality
% principles based on straightforward ergodic or allostatic principles.
% In other words, to maintain the constancy of our external milieu, it is
% sufficient to expose ourselves to predicted and predictable stimuli.
% Being able to predict what is currently seen also enables us to predict
% fictive sensations that we will experience from another viewpoint. This
% provides a principled way in which to explore and sample the world for
% example with visual searches using saccadic eye movements. These
% theoretical considerations are remarkably consistent with a number
% of compelling heuristics; most notably the Infomax principle or the
% principle of minimum redundancy, signal detection theory and recent
% formulations of salience in terms of Bayesian surprise. The example
% here uses saliency (the posterior precision associated with fictive
% sampling of sensory data) to simulate saccadic eye movements under
% active inference.
%__________________________________________________________________________
% Copyright (C) 2011 Wellcome Trust Centre for Neuroimaging

% Karl Friston
% $Id: DEM_evidence_accumulation.m 6587 2015-11-02 10:29:49Z karl $


% hidden causes and states
%==========================================================================
% x    - hidden states:
%   o(1) - oculomotor angle
%   o(2) - oculomotor angle
%   x(1) - relative amplitude of visual hypothesis 1
%   x(2) - relative amplitude of visual hypothesis 2
%   x(3) - ...
%
% v    - hidden causes
%
% g    - sensations:
%   g(1) - oculomotor angle (proprioception - x)
%   g(2) - oculomotor angle (proprioception - y)
%   g(3) - retinal input - channel 1
%   g(4) - retinal input - channel 2
%   g(5) - ...
%--------------------------------------------------------------------------


% generative model
%==========================================================================
M(1).E.s = 1/2;                                % smoothness
M(1).E.n = 4;                                  % order of
M(1).E.d = 1;                                  % generalised motion

% level 1: Displacement dynamics and mapping to sensory/proprioception
%--------------------------------------------------------------------------
M(1).f  = inline('([1 - sum(exp(x)) - (x + log(3))/32 - [0; 1; 0]*v])/32','x','v','P');
M(1).g  = inline('[ [-1 0 1]*spm_phi(16*(exp(x) - 1)); [-1 0 1]*exp(x) + zeros(16,1)]','x','v','P');
M(1).x  = -log(3)*[1; 1; 1];                   % hidden states
M(1).V  = [0 (zeros(1,16) + exp(-4))];         % error precision (g)
M(1).W  = exp(32);                             % error precision (f)
% M(1).xP = exp(2);                            % error precision (f)

% level 2:
%--------------------------------------------------------------------------
M(2).v  = 0;                                   % priors
M(2).V  = exp(8);


% generative process
%==========================================================================


% first level
%--------------------------------------------------------------------------
G(1).g  = inline('[a; v + zeros(16,1)]','x','v','a','P');
G(1).V  = [exp(8) zeros(1,16) + exp(-2)];      % error precision
G(1).U  = [1 zeros(1,16)];                     % motor gain

% second level
%--------------------------------------------------------------------------
G(2).v  = 0;                                   % exogenous forces
G(2).a  = 0;                                   % action forces
G(2).V  = exp(16);


% generate and invert
%==========================================================================
N     = 32;                                    % length of data sequence
pst   = (1:N);                                 % perstimulus time (bins)

DEM.G = G;
DEM.M = M;
DEM.C = spm_phi(((1:N) - N/4)*32/N);
DEM.U = sparse(1,N);

DEM     = spm_ADEM(DEM);
spm_DEM_qU(DEM.qU,DEM.pU)

subplot(2,2,2)
spm_plot_ci(DEM.qU.x,DEM.qU.S,pst,3,'exp'), hold on
plot(pst,exp(DEM.qU.x{1}),'LineWidth',1),  hold off
axis([1 N 0 2])
axis square







return
    


    % (k) trial
    %----------------------------------------------------------------------
    for k = 1:8
        
        % solve and save saccade
        %------------------------------------------------------------------
        DEM     = spm_ADEM(DEM);
        DEM     = spm_ADEM_update(DEM);
        
        % overlay true values
        %------------------------------------------------------------------
        spm_DEM_qU(DEM.qU,DEM.pU)
        
        % compute salience
        %------------------------------------------------------------------
        [S L]   = spm_salience_map(DEM.M,nr);
        
        % optimise prior belief
        %------------------------------------------------------------------
        S       = (S - min(S)).*(1 - R);
        [i j]   = max(S);
        DEM.U   = L(:,j)*ones(1,N);
        DEM.S   = reshape(S,nr,nr);
        
        % inhibition of return (IOR)
        %------------------------------------------------------------------
        D       = exp(-sum((L - L(:,j)*ones(1,nr*nr)).^2)/(2*s*s))';
        R       = a*R + D;
        
        % store
        %------------------------------------------------------------------
        ADEM{k} = DEM;
        
    end
    
    % save
    %----------------------------------------------------------------------
    save ADEM_saccades ADEM
    

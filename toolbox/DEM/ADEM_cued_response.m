function ADEM_cued_response
% Cued responses under active inference: 
%__________________________________________________________________________
% This demo illustrates cued sequential movements. It uses active inference
% under a hierarchal generative model of sequential cues and consequent
% movements. The agent has a (second level) model of (two) contingencies or
% contexts; these correspond to the sequential appearance of targets in a
% clockwise direction. The other context has no sequential aspect. The
% first level model is contextually modulated to produce the appropriate
% sequence of (location – specific) affordances, which predict both
% visual and proprioceptive consequences. This is sufficient to engender
% cued reaching movements, which are slightly anticipatory if the agent
% infers the correct probabilistic context. However, if we reverse the
% order of the stimuli there is an accuracy and reaction time cost, due to
% the fact that the sequence is unpredictable.  Furthermore, there is a
% set switching cost as the hidden states at the second (contextual) level
% are inferred. This provides a simple but very rich model of cued reaching
% movements and set switching that is consistent with notions of salience
% and affordance. Furthermore, we can simulate Parkinsonism by
% reducing the precision of affordance – based cues. These are the visual
% attributes that confer saliency on the current target. Reducing this
% precision (for example, dopamine) delays and can even preclude set
% switching, with associated costs in pointing accuracy. By completely
% removing the precision of the salience or affordance cues, we obtain
% autonomous behaviour that is prescribed by the itinerant expectations of
% the agent. This can be regarded as perseveration in a pathological
% setting or the emission of autonomous behaviour in the absence of any
% precise sensory information
%__________________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging
 
% Karl Friston
% $Id: ADEM_cued_response.m 6290 2014-12-20 22:11:50Z karl $
 
 
% hidden causes and states
%==========================================================================
% x    - hidden states:
%   x.o  - motor angle (intrinsic frame of reference)
%   x.a  - target salience
%
% v    - hidden causes
%
% P    - parameters
%
% g    - sensations:
%   g.o  - motor angle (proprioception)
%   g.p  - finger location (visual)
%   g.c  - target contrast (visual)
%--------------------------------------------------------------------------

% parameters mapping from (unstable) point attractors to visual space
%--------------------------------------------------------------------------
P.x    = [-1  1  1 -1;                       % target locations (extrinsic)
          -1 -1  1  1];
 
P.s    = 1;                                  % softmax parameter
n      = size(P.x,2);                        % number of attractors
N      = 128;                                % length of data sequence
 
% hidden states (M)
%--------------------------------------------------------------------------
x.o    = sparse(2,1);                        % finger location
x.a    = sparse(n,1) - 4;                    % attractor (SHC) states
x.a(4) = 0;
 
 
% Recognition model
%==========================================================================
M(1).E.s = 1;                                % smoothness
M(1).E.n = 3;                                % order of 
M(1).E.d = 2;                                % generalised motion
 
 
% precisions: sensory input
%--------------------------------------------------------------------------
V.o = exp(4) + sparse(2,1);                  % motor (proprioceptive)
V.p = exp(4) + sparse(2,1);                  % target locations (visual)
V.c = exp(4) + sparse(n,1);                  % target salience  (visual)

W.o = exp(4) + sparse(2,1);                  % motor (proprioceptive)
W.a = exp(4) + sparse(n,1);                  % target salience  (visual)
 
 
% level 1: Displacement dynamics and mapping to sensory/proprioception
%--------------------------------------------------------------------------
M(1).f  = 'spm_fx_dem_cue';                  % plant dynamics
M(1).g  = 'spm_gx_dem_cue';                  % prediction
M(1).pE = P;                                 % target locations
M(1).x  = x;                                 % hidden states
M(1).V  = diag(spm_vec(V));                  % error precision
M(1).W  = diag(spm_vec(W));                  % error precision
 
% level 2
%--------------------------------------------------------------------------
M(2).f  = inline('spm_lotka_volterra(x,1/16)', 'x', 'v', 'P');
M(2).g  = inline('spm_softmax(x)',             'x', 'v', 'P');
 
M(2).x  = [2; -2];                           % hidden states
M(2).v  = [1;  0];                           % hidden causes
M(2).V  = exp(4);                            % error precision
M(2).W  = exp(4);                            % error precision
 
% generative model
%==========================================================================
 
% hidden states (G)
%--------------------------------------------------------------------------
x.o = [0;0];                                 % finger location
x.a = sparse(n,1);                           % attractor (SHC) states
 
% first level
%--------------------------------------------------------------------------
G(1).f  = 'spm_fx_adem_cue';
G(1).g  = 'spm_gx_adem_cue';
G(1).pE = P.x;
G(1).x  = x;                                 % hidden states
G(1).V  = exp(16);                           % error precision
G(1).W  = exp(16);                           % error precision
G(1).U  = [1 1 0 0 0 0 0 0]*exp(2);          % error precision
 
% second level
%--------------------------------------------------------------------------
G(2).v  = sparse(n,1);                       % exogenous cues
G(2).a  = [0; 0];                            % action forces
G(2).V  = exp(16);
 
% Create cues
%==========================================================================
isi   = 12;                                  % inter stimulus interval
dur   = 4;                                   % stimulus duration
j     = [1 2 3 4 1 4 3 2 1 4 3 2 1];         % stimulus order
 
C     = sparse(n,N);
for i = 1:length(j)
    C(j(i),:) = C(j(i),:) + exp(-([1:N] - i*isi).^2/(2*dur^2));
end
 
% generate and invert of over different DA levels (D = V.c)
%==========================================================================
DEM.G = G;
DEM.M = M;
DEM.C = C;

SIM   = 'salience'
D     = 2 + (1:6)/2;
for i = 1:length(D)
    
    switch SIM
        
        case {'salience'}
            % set precision of salience or affordance - V(1).c
            %--------------------------------------------------------------
            ADEM{i}        = DEM;
            V.c            = exp(D(i)) + sparse(n,1);
            ADEM{i}.M(1).V = diag(spm_vec(V));
            
        case {'proprioception'}
            % set precision of salience or affordance - W(1).o
            %--------------------------------------------------------------
            ADEM{i}        = DEM;
            W.o            = exp(D(i)) + sparse(2,1);
            ADEM{i}.M(1).W = diag(spm_vec(W));
            
        case {'affordance'}
            % set precision of salience or affordance - W(1).a
            %--------------------------------------------------------------
            ADEM{i}        = DEM;
            W.a            = exp(D(i)) + sparse(n,1);
            ADEM{i}.M(1).W = diag(spm_vec(W));
            
        case {'switching'}
            % set precision of salience or affordance - V(2)
            %--------------------------------------------------------------
            ADEM{i}        = DEM;
            ADEM{i}.M(2).V = exp(D(i));
            
    end
    
    
    % solve
    %----------------------------------------------------------------------
    ADEM{i} = spm_ADEM(ADEM{i});
    
end


 
% dynamics and movies
%==========================================================================
 
% normal
%--------------------------------------------------------------------------
spm_figure('GetWin','Figure 1');i  = 6;
spm_DEM_qU(ADEM{i}.qU)
 
spm_figure('GetWin','Figure 4');
subplot(3,1,1)
spm_dem_cue_movie(ADEM{i})
title(sprintf('High: (%-0.1f)',D(i)),'FontSize',16)
 
% delayed set switching
%--------------------------------------------------------------------------
spm_figure('GetWin','Figure 2'); i = 3;
spm_DEM_qU(ADEM{i}.qU)
 
spm_figure('GetWin','Figure 4');
subplot(3,1,2)
spm_dem_cue_movie(ADEM{i})
title(sprintf('Intermediate: (%-0.1f)',D(i)),'FontSize',16)
 
% failure to switch set
%--------------------------------------------------------------------------
spm_figure('GetWin','Figure 3'); i = 1;
spm_DEM_qU(ADEM{2}.qU)
 
spm_figure('GetWin','Figure 4');
subplot(3,1,3)
spm_dem_cue_movie(ADEM{i})
title(sprintf('Low: (%-0.1f)',D(i)),'FontSize',16)
 
 
% reaction times and accuracy
%==========================================================================
spm_figure('GetWin','Figure 5'); clf
 
C     = {'r','g','b','c','m','y'};
for i = 1:length(D)
 
    % get behavioral performance
    %----------------------------------------------------------------------
    [on,rt,ac] = spm_ADEM_cue_rt(ADEM{i});
 
    % graphics
    %----------------------------------------------------------------------
    subplot(2,1,1), hold on
    plot(on,rt,[C{i} '-'],on,rt,[C{i} '.'],'MarkerSize',24)
    xlabel('cue onset (sec)','FontSize',12)
    ylabel('milliseconds','FontSize',12)
    title('reaction times','FontSize',16)
    axis square, box off, hold off
 
    subplot(2,1,2), hold on
    plot(on,ac,[C{i} '-'],on,ac,[C{i} '.'],'MarkerSize',24)
    xlabel('cue onset (sec)','FontSize',12)
    ylabel('inverse distance','FontSize',12)
    title('spatial error','FontSize',16)
    axis square, box off, hold off
    drawnow
 
end
 

return


% demonstrate autonomous movement (by removing salient information)
%==========================================================================
LDEM        = DEM;
V.c         = exp(-2) + sparse(n,1);
LDEM.M(1).V = diag(spm_vec(V));
LDEM        = spm_ADEM(LDEM);
 
% perseverative dynamics and trajectories
%--------------------------------------------------------------------------
spm_figure('GetWin','Figure 6'); clf
spm_DEM_qU(LDEM.qU)
 
spm_figure('GetWin','Figure 8');
subplot(3,1,1)
spm_dem_cue_movie(LDEM)
title('Autonomous movements','FontSize',16)
 
 
% now remove precision on empirical priors of motion
%==========================================================================
LDEM.M(1).W = exp(-3);
LDEM        = spm_ADEM(LDEM);
 
% perseverative dynamics and trajectories
%--------------------------------------------------------------------------
spm_figure('GetWin','Figure 7'); clf
spm_DEM_qU(LDEM.qU)
 
spm_figure('GetWin','Figure 8');
subplot(3,1,2)
spm_dem_cue_movie(LDEM)
title('Confusion','FontSize',16)

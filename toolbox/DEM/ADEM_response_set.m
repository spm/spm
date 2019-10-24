% Cued responses under active inference: 
%__________________________________________________________________________
% This demo illustrates cued sequential movements. It uses active inference
% under a hierarchal generative model of sequential cues and consequent
% movements. The agent has a (second level) model of (two) contingencies or
% contexts; these correspond to the sequential appearance of targets in a
% clockwise direction. The other context has no sequential aspect. The
% first level model is contextually modulated to produce the appropriate
% sequence of (location - specific) affordances, which predict both
% visual and proprioceptive consequences. This is sufficient to engender
% cued reaching movements, which are slightly anticipatory if the agent
% infers the correct probabilistic context. However, if we reverse the
% order of the stimuli there is an accuracy and reaction time cost, due to
% the fact that the sequence is unpredictable.  Furthermore, there is a
% set switching cost as the hidden states at the second (contextual) level
% are inferred. This provides a simple but very rich model of cued reaching
% movements and set switching that is consistent with notions of salience
% and affordance. Furthermore, we can simulate Parkinsonic lesions by
% reducing the precision of affordance - based cues. These are the visual
% attributes that confer saliency on the current target. Reducing this
% precision (for example, dopamine) delays and can even preclude set
% switching, with associated costs in pointing accuracy. By completely
% removing the precision of the salience or affordance cues, we obtain
% autonomous behavior that is prescribed by the itinerant expectations of
% the agent. This can be regarded as perseveration in a pathological
% setting or the emission of autonomous behavior in the absence of any
% precise sensory information
%__________________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging
 
% Karl Friston
% $Id: ADEM_response_set.m 7679 2019-10-24 15:54:07Z spm $
 
 
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
%   g.p  - target locations (visual) - intrinsic coordinates (polar)
%   g.c  - target contrast  (visual)
%--------------------------------------------------------------------------
M      = struct;
G      = struct;

% parameters mapping from (unstable) point attractors to visual space
%--------------------------------------------------------------------------
P      = [-1  1  1 -1;                       % target locations (extrinsic)
          -1 -1  1  1];
 
n      = size(P,2);                          % number of attractors
N      = 128;                                % length of data sequence
 
% Recognition model
%==========================================================================
M(1).E.s = 1;                                % smoothness
M(1).E.n = 3;                                % order of 
M(1).E.d = 2;                                % generalised motion

% precisions (sensory)
%--------------------------------------------------------------------------
V.o     = exp(8)  + sparse(2,1);
V.p     = exp(8)  + sparse(2,1);
V.c     = exp(8)  + sparse(n,1);


% level 1: Displacement dynamics and mapping to sensory/proprioception
%--------------------------------------------------------------------------
M(1).f  = inline('(P*v - tan(x))/4', 'x', 'v', 'P');
M(1).g  = inline('[x; tan(x); exp(-sum((P - tan(x)*ones(1,size(P,2))).^2)/2)'']', 'x', 'v', 'P');
M(1).pE = P;                                 % target locations
M(1).x  = [0; 0];                            % hidden states (movement)
M(1).V  = diag(spm_vec(V));                  % error precision
M(1).W  = exp(8);                            % error precision


% level 2: Displacement dynamics and mapping to sensory/proprioception
%--------------------------------------------------------------------------
M(2).f  = inline('spm_lotka_volterra(x,v(1)/2)', 'x', 'v', 'P');
M(2).g  = inline('spm_softmax(x,2)',             'x', 'v', 'P');
M(2).pE = P;                                 % target locations
M(2).x  = [0; -1; -6; 2];                    % hidden states (affordance)
M(2).V  = exp(6);                            % error precision
M(2).W  = exp(4);                            % error precision
M(2).xP = exp(0);

 
% level 3
%--------------------------------------------------------------------------
M(3).f  = inline('spm_lotka_volterra(x,1/64)', 'x', 'v', 'P');
M(3).g  = inline('spm_softmax(x)',             'x', 'v', 'P');
 
M(3).x  = [2; -2];                           % hidden states (set)
M(3).V  = exp(2);                            % error precision
M(3).W  = exp(0);                            % error precision

 
% generative model
%==========================================================================

% first level
%--------------------------------------------------------------------------
G(1).f  = inline('tanh(a) - x/8',  'x', 'v', 'a', 'P');
G(1).g  = inline('[x; tan(x); v]', 'x', 'v', 'a', 'P');
G(1).pE = P;
G(1).x  = [0; 0];                            % hidden states (movement)
G(1).V  = exp(16);                           % error precision
G(1).W  = exp(16);                           % error precision
 
 
% second level
%--------------------------------------------------------------------------
G(2).v  = sparse(n,1);                       % exogenous cues (salience)
G(2).a  = [0; 0];                            % action
G(2).V  = exp(16);
 
% Create cues
%==========================================================================
dt    = 64;                                  % bin size (ms)
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


DEM   = spm_ADEM(DEM);
subplot(4,2,7)
spm_dem_set_movie(DEM)

return


% graphics
%----------------------------------------------------------------------
[on,rt,ac] = spm_ADEM_set_rt(DEM);

subplot(2,1,1), hold on
plot(on,rt,'r:',on,rt,'r.','MarkerSize',24)
xlabel('cue onset (sec)','FontSize',12)
ylabel('milliseconds','FontSize',12)
title('reaction times','FontSize',16)
axis square, box off, hold off

% set precision of salience or affordance
%--------------------------------------------------------------------------
D     = (-3:3)/2;
D     = -2:2; 
for i = 1:length(D)
 
    ADEM{i}        = DEM;
    ADEM{i}.M(3).V = ADEM{i}.M(3).V*exp(D(i));
    ADEM{i}        = spm_ADEM(ADEM{i});   
    
end

 
% reaction times and accuracy
%==========================================================================
spm_figure('GetWin','Figure 4'); clf
 
C     = {'r','g','b','c','m'};
for i = 1:length(C)
 
    % get behavioral performance
    %----------------------------------------------------------------------
    [on,rt,ac] = spm_ADEM_set_rt(ADEM{i});
 
    % graphics
    %----------------------------------------------------------------------
    subplot(2,1,1), hold on
    plot(on,rt,[C{i} ':'],on,rt,[C{i} '.'],'MarkerSize',24)
    xlabel('cue onset (sec)','FontSize',12)
    ylabel('milliseconds','FontSize',12)
    title('reaction times','FontSize',16)
    axis square, box off, hold off
 
    subplot(2,1,2), hold on
    plot(on,ac,[C{i} ':'],on,ac,[C{i} '.'],'MarkerSize',24)
    xlabel('cue onset (sec)','FontSize',12)
    ylabel('inverse distance','FontSize',12)
    title('spatial error','FontSize',16)
    axis square, box off, , hold off, drawnow
 
end

return
 
% dynamics and movies
%==========================================================================
 
% normal
%--------------------------------------------------------------------------
spm_figure('GetWin','Figure 1');
spm_DEM_qU(ADEM{6}.qU)
 
subplot(3,2,5)
spm_dem_set_movie(ADEM{6})
 
% delayed set switching
%--------------------------------------------------------------------------
spm_figure('GetWin','Figure 2');
spm_DEM_qU(ADEM{3}.qU)
 
subplot(3,2,5)
spm_dem_set_movie(ADEM{3})
 
% failure to switch set
%--------------------------------------------------------------------------
spm_figure('GetWin','Figure 3');
spm_DEM_qU(ADEM{2}.qU)
 
subplot(3,2,5)
spm_dem_set_movie(ADEM{2})


% demonstrate autonomous movement (by removing salient information)
%==========================================================================
DEM        = DEM;
V.c        = exp(0) + sparse(n,1);
DEM.M(1).V = diag(spm_vec(V));
DEM        = spm_ADEM(DEM);
 
% perseverative dynamics and trajectories
%--------------------------------------------------------------------------
spm_figure('GetWin','DEM');
spm_DEM_qU(DEM.qU)
 
subplot(3,2,5)
spm_dem_set_movie(DEM)

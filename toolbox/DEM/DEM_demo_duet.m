function DEM = DEM_demo_duet
% Demo for bird songs: In this example, we show that DEM can not only
% estimate the hidden states of an autonomous system but can also
% deconvolve dynamic changes in its control parameters.  We illustrate
% this using a slow Lorentz attractor to drive a fast one; both showing
% deterministic chaos.  We endow the simulations with a little ethological
% validity by using the states of the fast Lorentz attractor as control
% variables in the syrinx of a song bird (usually these would control a van
% der Pol oscillator model). We will look at the true and inferred songs
% with and without the last chirps missing.  The sonograms displayed
% can be played by a mouse click on the image.  Subsequent plots show
% simulated event-related potential to show that there is a marked
% responses (prediction error) of the system when an expected ‘syllable’ is
% omitted. This demonstrates the implicit sequence-decoding of input
% streams, using generative models based upon attractors.
% Having simulated normal omission-related responses, we then reduce the
% precision at the second level (on both hidden causes and states) and
% repeat the simulation. The result is an attenuation of the omission-
% related response or mismatch negativity. If we try to compensate by
% reducing the sensory precision, then the autonomous dynamics predicting
% the sequence of chirps supervenes, producing false inference. This
% can be thought of as a – crude - model of hallucinosis.
%__________________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging
 
% Karl Friston
% $Id: DEM_demo_duet.m 6039 2014-06-04 18:50:28Z karl $
 
rng('default')

A   = 2;                                   % number of agents (birds)
dt  = 1/64;                                % time bin (seconds)

% Hierarchical non-linear generative model (dynamic & chaotic)
%==========================================================================
M(1).E.d        = 1;                       % approximation order
M(1).E.n        = 3;                       % embedding order
M(1).E.s        = 1/4;                     % smoothness
 
M(1).E.method.x = 0;                       % state-dependent noise
M(1).E.method.v = 0;                       % state-dependent noise
M(1).E.method.h = 0;                       % suppress optimisation
M(1).E.method.g = 0;                       % suppress optimisation

P     = cell(A,1);
x     = cell(A,1);
a     = cell(A,1);
U     = cell(1,A);
V     = cell(1,A);

x0    = [1; 1];

Up    = exp([ 8  8 -8 -8]);                % sensory attenutation
Uq    = exp([-8 -8 -8 -8]);                % sensory attention
Vp    = exp([-8 -8  0  0]);                % attenutation
Vq    = exp([-8 -8  2  2]);                % attention

for i = 1:A
    if rem(i,2), U{i} = Up; else, U{i} = Uq; end
    if rem(i,2), V{i} = Vp; else, V{i} = Vq; end
end

for i = 1:A
    x{i} = x0;                             % states (of syrinx)
end
G(1).f  = @(x,v,a,P) a - spm_vec(x)/4;
G(1).g  = @(x,v,a,P) Gg1(x,v,a,P);         % SOUND PRODUCTION
G(1).x  = x;                               % hidden state
G(1).V  = exp(8);                          % precision (noise)
G(1).W  = exp(8);                          % precision (states)
G(1).U  = spm_cat(U);                      % precision (action)
 
 
% level 2; causes
%--------------------------------------------------------------------------
for i = 1:A
    a{i} = [0; 0];                         % action (fequency and volume)
end
G(2).v  = 0;                               % exogenous  cause
G(2).a  = spm_vec(a);                      % endogenous cause (action)
G(2).V  = exp(8);

% level 1
%--------------------------------------------------------------------------
for i = 1:A
    P{i} = [10; 8/3];                      % parameters
    x{i} = [1; x0];                        % hidden states
end
M(1).f  = @(x,v,P) Mf1(x,v,P);
M(1).g  = @(x,v,P) Mg1(x,v,P);
M(1).x  = x;
M(1).pE = P;
M(1).W  = exp(2);
M(1).V  = spm_cat(V); 

% level 2
%--------------------------------------------------------------------------
for i = 1:A
    P{i} = [10; 8/3];                      % parameters
    x{i} = [1; 1; 30];                     % hidden states
   
    if i > 1,
        x{i} = [1; 1; 1];                 % hidden states
    end
end
M(2).f  = @(x,v,P) Mf2(x,v,P);
M(2).g  = @(x,v,P) Mg2(x,v,P);
M(2).x  = x;
M(2).v  = 0;
M(2).pE = P;
M(2).V  = exp(8);
M(2).W  = exp(8);
 
 
% hidden cause and prior expectations
%--------------------------------------------------------------------------
T     = 2;                                 % number of trials
N     = 128;                               % length of stimulus (bins)
C     = zeros(1,N);

% assemble model structure
%--------------------------------------------------------------------------
DEM.M = M;
DEM.G = G;
DEM.C = C;

% reset initial hidden states and invert
%==========================================================================
for t = 1:T
    
    DEM    = spm_ALAP(DEM);
    LAP{t} = DEM;
    spm_DEM_qU(LAP{t}.qU,LAP{t}.pU)
    
    % update percicions and switch roles
    %----------------------------------------------------------------------
    if rem(t,2)
        for i = 1:A
            if rem(i,2), U{i} = Uq; else, U{i} = Up; end
            if rem(i,2), V{i} = Vq; else, V{i} = Vp; end
        end
    else
        for i = 1:A
            if rem(i,2), U{i} = Up; else, U{i} = Uq; end
            if rem(i,2), V{i} = Vp; else, V{i} = Vq; end
        end
    end
    DEM.G(1).U  = spm_cat(U);
    DEM.M(1).V  = spm_cat(V);
    
    
    % update states
    %----------------------------------------------------------------------
    DEM   = spm_ADEM_update(DEM);
    
end


% illustrate responses with sonogram (and sound file)
%--------------------------------------------------------------------------
spm_figure('GetWin','Figure 1'); clf

qU    = [];
for t = 1:T
    v = LAP{t}.qU.v{1}([1 2],:);
    if rem(t - 1,2)
        v(1,:) = v(1,:)/4;
    end
    qU = [qU v];
end

subplot(3,1,1)
colormap('pink')
spm_DEM_play_song(qU,T*N*dt);
title('percept','Fontsize',16)

subplot(3,1,2)
for i = 1:A
    x = [];
    for t = 1:T
        x = [x LAP{t}.qU.x{1}([1 2 3] + (i - 1)*3,:)];
    end
    if i > 1
        plot((1:size(x,2))*dt,x,'c'),hold on
    else
        plot((1:size(x,2))*dt,x,'r'),hold on
    end
end
xlabel('time (seconds)')
title('First level expectations (hidden states)','Fontsize',16)

subplot(3,1,3)
for i = 1:A
    x = [];
    for t = 1:T
        x = [x LAP{t}.qU.x{2}([1 2 3] + (i - 1)*3,:)];
    end
    if i > 1
        plot((1:size(x,2))*dt,x,'c'),hold on
    else
        plot((1:size(x,2))*dt,x,'r'),hold on
    end
end
xlabel('time (seconds)')
title('Second level expectations (hidden states)','Fontsize',16)



% Equations of motion and observer functions
%==========================================================================
function g = Gg0(x,v,a,P)
for i = 1:length(x)
    g{i,1} = [x{i}; x{i}];
end

function g = Gg1(x,v,a,P)
for i = 1:length(x)
    s(i)   = abs(x{i}(1));
end
[a j] = max(s);
for i = 1:length(x)
    g{i,1} = [x{i}; x{j}];
end

function g = Mg1(x,v,P)
for i = 1:length(x)
    g{i,1} = x{i}([2 3 2 3]);
end

function g = Mg2(x,v,P)
for i = 1:length(x)
    g{i,1} = x{i}(3);
end

function f = Mf1(x,v,P)
for i = 1:length(x)
    f{i,1} = [-P{i}(1) P{i}(1) 0; (v{i}(1) - 4 - x{i}(3)) -1 0; x{i}(2) 0 -P{i}(2)]*x{i}/16;
end

function f = Mf2(x,v,P)
for i = 1:length(x)
    f{i,1} = [-P{i}(1) P{i}(1) 0; (32 - x{i}(3)) -1 0; x{i}(2) 0 -P{i}(2)]*x{i}/128;
end



 
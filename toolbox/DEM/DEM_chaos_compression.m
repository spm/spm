function MDP = DEM_chaos_compression
% Structure learning from pixels: the Lorenz attractor
%__________________________________________________________________________
%
% This demonstration routine uses a renormalisable generative model to
% learn, recognise and generate stochastic chaos. It uses the same
% procedures used to learn paths and orbits from image sequences, apt for
% video compression. Effectively, it learns stochastic chaos based upon
% exponential divergence of trajectories in terms of a switching dynamical
% system, where successive paths are selected based upon a probability
% transition matrix. Here, an image of a white ball moving on the
% attracting manifold of a Lorenz system is generated. This is then
% quantised and learned in the usual way by fast structure learning.
% Subsequent (active) learning then fills in probability transitions, which
% are subsequently refined using (simple) Bayesian model reduction. The
% ensuing generative model can then be used to filter images (i.e.,
% recognise the phase of the underlying dynamics), and generate stochastic
% chaos when the stimulus is removed.
%__________________________________________________________________________
% Copyright (C) 2019 Wellcome Trust Centre for Neuroimaging

% Karl Friston
% $Id: DEM_dSprites.m 8447 2023-06-17 16:32:24Z karl $

%% set up and preliminaries
%--------------------------------------------------------------------------
rng(1)

% get model of stochastic chaos
%==========================================================================
M    = spm_DEM_M('Lorenz');

% create innovations & add causes
%--------------------------------------------------------------------------
N    = 1024;
U    = sparse(1,N);
DEM  = spm_DEM_generate(M,U);

% show realisation
%==========================================================================
spm_figure('GetWin','Lorenz'); clf
spm_DEM_qU(DEM.pU);

% create a pretty image of Lorenzian trajectories
%--------------------------------------------------------------------------
subplot(2,3,4), hold off
x     = DEM.pU.x{1};
plot(x(1,:),x(2,:),':','Color',[.9 .7 0]), hold on
h = plot(x(1,1),x(2,1),'.w','MarkerSize',64); set(gca,'Color','k')
axis image, axis([-24 24 -24 24])
for t = 1:N
    set(h,'Xdata',x(1,t),'Ydata',x(2,t))
    I(:,:,:,t) = frame2im(getframe(gca));
end

% Map from image to discrete state space (c.f., Amortisation) 
%--------------------------------------------------------------------------
RGB.nd    = 32;                    % Diameter of tiles in pixels
RGB.nb    = 5;                     % Number of discrete singular variates 
RGB.mm    = 16;                    % Maximum number of singular modes
RGB.su    = 16;                    % Variance threshold
RGB.R     = 2;                     % temporal resampling

T         = N/RGB.R;               % number of voxels             
i         = (1:(N/2));             % training set
[O,L,RGB] = spm_rgb2O(I(:,:,:,i),RGB);

% And show the images generated from a discrete representation
%--------------------------------------------------------------------------
subplot(2,2,3)
spm_imshow(I(:,:,:,1:16))
axis on, title('Original image','FontSize',12)
subplot(2,2,4)
spm_imshow(spm_O2rgb(O(:,1:8),RGB))
axis on, title('Discretised image','FontSize',12)

% Use the ensuing sequence for (RG) structure learning
%--------------------------------------------------------------------------
t     = 1:(T/2);
MDP   = spm_MB_structure_learning(O(:,t),L);

% learn from subsequent episodes
%==========================================================================
Nm    = numel(MDP);
FIX.A = 1;                             %  enable likelihood learning
FIX.B = 0;                             % disable transition learning

% active learning
%--------------------------------------------------------------------------
i     = (1:N/2) + N/2 - 64;
S     = spm_rgb2O(I(:,:,:,i),RGB);

RDP   = spm_mdp2rdp(MDP,0,1/512,2,FIX);
RDP.T = fix((T/2)/(2^(Nm - 1)));

RDP   = spm_RDP_O(RDP,S);
PDP   = spm_MDP_VB_XXX(RDP);
RDP   = spm_RDP_update(RDP,PDP,'SIMPLE');

% report
%--------------------------------------------------------------------------
spm_figure('GetWin','Active learning'); clf
spm_show_RGB(PDP,RGB,8,0);


% Recognition and generation
%==========================================================================

% Create partial stimulus
%--------------------------------------------------------------------------
S     = O(:,1:(T/4));
for g = 1:size(S,1)
    for t = 1:size(S,2)
        if t < T/8
            S{g,t} = O{g,t + T/8};
        else
            S{g,t} = spm_dir_norm(ones(size(S{g,t})));
        end
    end
end

RDP   = spm_RDP_O(RDP,S);
PDP   = spm_MDP_VB_XXX(RDP);

% Illustrate stimulus and posterior predictons
%--------------------------------------------------------------------------
spm_figure('GetWin','Evoked responses'); clf
[I,J] = spm_show_RGB(PDP,RGB,8,1);

I = max(I,[],2);
I = permute(I,[1 4 3 2]);
J = max(J,[],2);
J = permute(J,[1 4 3 2]);

spm_figure('GetWin','Posterior predictions'); clf
subplot(4,1,1);
spm_imshow(I), axis on
title('Predictive posterior','FontSize',12)

subplot(4,1,2);
spm_imshow(J), axis on
title('Obeservaton likelihood','FontSize',12)
xlabel('time (bins)')

return





function spm_dem_search_plot(DEM)
% plots visual search in extrinsic and intrinsic coordinates
% FORMAT spm_dem_search_plot(DEM)
%
% DEM - {DEM} structures from visual search simulations
%
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
%__________________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% Karl Friston
% $Id: spm_dem_search_plot.m 4580 2011-12-02 20:22:19Z karl $


% Preliminaries
%--------------------------------------------------------------------------
clf, global STIM
N  = length(DEM);
S  = spm_read_vols(STIM.V);

% Stimulus
%======================================================================
Dx = STIM.V.dim(1)/2;
Dy = STIM.V.dim(2)/2;
a  = [];
q  = [];
c  = [];

for i = 1:N
    
    % i-th saccade - position
    %----------------------------------------------------------------------
    pU = DEM{i}.pU.x{1}(1:2,:)*16;
    qU = DEM{i}.qU.x{1}(1:2,:)*16;
    T  = length(pU);
    
    % conditional confidence
    %----------------------------------------------------------------------
    qC = DEM{i}.qU.S;
    for t = 1:length(qC)
        qV(t) = qC{t}(3,3);
    end
    
    % accumulate responses
    %----------------------------------------------------------------------
    a  = [a DEM{i}.qU.a{2}];                % action
    q  = [q DEM{i}.qU.x{1}(3:end,:)];       % hidden perceptual states
    c  = [c qV];                            % conditional variance
    
    
    % eye movements in extrinsic coordinates
    %======================================================================
    subplot(6,N,i)
    
    image((S + 1)*32), axis image off, hold on
    plot(qU(2,T) + Dy,qU(1,T) + Dx,'.g','Markersize',8)
    plot(pU(2,T) + Dy,pU(1,T) + Dx,'.r','Markersize',16)
    drawnow, hold off
    
    
    % i-th saccade - sensory samples
    %----------------------------------------------------------------------
    pU = DEM{i}.pU.v{1}(3:end,:);
    
    % sensory input
    %======================================================================
    subplot(6,N,i + N*2)
    s = spm_unvec(pU(:,T),STIM.R);
    imagesc(s), axis image off, drawnow
    
    
    % i-th saccade - percept
    %----------------------------------------------------------------------
    qU = DEM{i}.qU.x{1}(3:end,:);
    
    % percept
    %======================================================================
    subplot(6,N,i + N*4)
    
    % hypotheses (0 < H < 1 = normalised neg-entropy)
    %--------------------------------------------------------------------------
    h     = spm_softmax(qU(:,T),2);
    H     = 1 + h'*log(h)/log(length(h));
    
    % retinotopic predictions
    %--------------------------------------------------------------------------
    s     = 0;
    for j = 1:length(h)
        s = s + h(j)*spm_read_vols(STIM.H{j});
    end
    image(s*H*64), axis image off, drawnow
    
end

% set ButtonDownFcn
%--------------------------------------------------------------------------
t  = (1:length(a));
subplot(6,1,2)
plot(a')
title('Action (EOG)','FontSize',16)
xlabel('time')

subplot(6,1,4)
spm_plot_ci(t,q(1,:),c); hold on
plot(t,q), hold off
set(gca,'YLim',[-1 1]*8)
title('Posterior beleif (EOG)','FontSize',16)
xlabel('time')

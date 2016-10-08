function spm_dem_mdp_movie(DEM)
% creates a movie of visual search in extrinsic and intrinsic coordinates
% FORMAT spm_dem_mdp_movie(DEM)
%
% DEM - {DEM} structures from visual search simulations
%
% hidden causes and states
%==========================================================================
% x    - hidden states:
%   x(1) - oculomotor angle
%   x(2) - oculomotor angle
% v    - hidden causes
%   v(1) - location of object
%   v(2) - location of object
%   v(3) - relative amplitude of visual hypothesis 1...
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
% $Id: spm_dem_mdp_movie.m 6901 2016-10-08 13:21:41Z karl $


% Preliminaries
%--------------------------------------------------------------------------
clf, global STIM
if ~isfield(STIM,'W'), STIM.W = 1/6; end
if ~isfield(STIM,'B'), STIM.B = 1;   end
B   = STIM.B;
W   = STIM.W;
R   = STIM.R;


% Stimulus
%==========================================================================
N     = length(DEM);
a     = [];
for j = 1:numel(DEM(1).X)
    X{j} = [];
end
for i = 1:N
    
    % get stimulus and position
    %--------------------------------------------------------------------------
    h   = find(DEM(i).C(3:end,1));
    P   = DEM(i).C(1:2,1);
    S   = spm_read_vols(STIM.H{h});
    
    Dx  = STIM.H{h}.dim(1)/2 - P(1)*16;
    Dy  = STIM.H{h}.dim(2)/2 - P(2)*16;
    
    dim = size(STIM.R);
    vox = STIM.H{h}.dim(1);
    dx  = vox/dim(1)*STIM.W;
    di  = dx*([1 dim(1)] - dim(1)/2) + Dx;
    dj  = dx*([1 dim(2)] - dim(2)/2) + Dy;
    di  = [di;di]; di = di(:);
    dj  = [dj;dj]';dj = dj(:);
    ax  = [(Dy - vox) (Dy + vox) (Dx - vox) (Dx + vox)];
    
    
    
    % i-th saccade - position
    %----------------------------------------------------------------------
    pU = DEM(i).pU.x{1}(1:2,:)*16;
    qU = DEM(i).qU.x{1}(1:2,:)*16;
    T  = length(pU);
    a  = [a DEM(i).qU.a{2}];
    for j = 1:numel(DEM(i).X)
        X{j} = [X{j} DEM(i).X{j}];
    end
    
    % eye movements in extrinsic coordinates
    %======================================================================
    subplot(3,2,1)
    for t = 1:T
        image((S + 1)*32), axis image, hold on
        plot(qU(2,t) + Dy,qU(1,t) + Dx,'.g','Markersize',8)
        plot(pU(2,t) + Dy,pU(1,t) + Dx,'.r','Markersize',16)
        plot(pU(2,t) + dj,pU(1,t) + di,'+' ,'Markersize',16)
        
        % show location of image ccentre
        %------------------------------------------------------------------
        plot(Dy,Dx,'+r','Markersize',32); axis(ax);
        
        drawnow, hold off
        
        % save
        %------------------------------------------------------------------
        Me(t + T*(i - 1)) = getframe(gca);
        
    end
    
    % sensory input
    %======================================================================
    subplot(3,2,2)
    STIM.B = B*B';
    STIM.W = W;
    STIM.R = R;
    for t = 1:T
        
        x   = DEM(i).pU.x{1}(:,t);
        v   = DEM(i).pU.v{2}(1:2,t);
        h   = DEM(i).pU.v{2}(3:end,t);
        s   = ADEM_sample_image(x - v,h);
        image(s*64), axis image, drawnow
        
        % save
        %----------------------------------------------------------------------
        Mi(t + T*(i - 1)) = getframe(gca);
        
    end
    
    % percept
    %======================================================================
    subplot(3,2,4)
    STIM.B = 1;
    STIM.W = W;
    STIM.R = ones(128,128);
    for t = 1:T
        
        x   = DEM(i).qU.x{1}(:,t);
        v   = DEM(i).qU.v{2}(1:2,t);
        h   = DEM(i).qU.v{2}(3:end,t);
        s   = ADEM_sample_image(x - v,h);
        image(s*64), axis image, drawnow
        
        % save
        %------------------------------------------------------------------
        Mq(t + T*(i - 1)) = getframe(gca);
        
    end
    
end

% reset feild of STIM
%--------------------------------------------------------------------------
STIM.B = B;
STIM.W = W;
STIM.R = R;

% set ButtonDownFcn
%--------------------------------------------------------------------------
subplot(3,2,3)
plot(a')
title('Action (EOG)','FontSize',16)
xlabel('time')
axis square

for j = 1:numel(X)
    subplot(3,2,4 + j)
    imagesc(X{j})
    title('Accumlated evidence','FontSize',16)
    xlabel('time')
    axis square
end

subplot(3,2,1)
set(gca,'Userdata',{Me,16})
set(gca,'ButtonDownFcn','spm_DEM_ButtonDownFcn')
title('saccades (click axis for movie)','FontSize',16)

subplot(3,2,2)
set(gca,'Userdata',{Mi,16})
set(gca,'ButtonDownFcn','spm_DEM_ButtonDownFcn')
title('samples (click axis for movie)','FontSize',16)

subplot(3,2,4)
set(gca,'Userdata',{Mq,16})
set(gca,'ButtonDownFcn','spm_DEM_ButtonDownFcn')
title('percept (click axis for movie)','FontSize',16)

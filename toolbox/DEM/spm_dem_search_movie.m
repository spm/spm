function spm_dem_search_movie(DEM)
% creates a movie of visual search in extrinsic and intrinsic coordinates
% FORMAT spm_dem_search_movie(DEM)
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

% Karl Friston
% Copyright (C) 2008-2022 Wellcome Centre for Human Neuroimaging


% Preliminaries
%--------------------------------------------------------------------------
clf, global STIM
if ~iscell(DEM), DEM = {DEM};           end
if ~isfield(STIM,'W'), STIM.W = 1/6;    end
if ~isfield(STIM,'P'), STIM.P = [0;0];  end
if ~isfield(STIM,'U'), STIM.U = STIM.V; end
if ~isfield(STIM,'S'), STIM.S = STIM.H; end

N   = length(DEM);
S   = spm_read_vols(STIM.U);

% Stimulus
%==========================================================================
Dx  = STIM.P(1)*16 + STIM.U.dim(1)/2;
Dy  = STIM.P(2)*16 + STIM.U.dim(2)/2;

dim = size(STIM.R);
vox = STIM.U.dim(1);
dx  = vox/dim(1)*STIM.W;
di  = dx*([1 dim(1)] - dim(1)/2) + Dx;
dj  = dx*([1 dim(2)] - dim(2)/2) + Dy;
di  = [di;di]; di = di(:);
dj  = [dj;dj]';dj = dj(:);
ax  = [(Dy - vox) (Dy + vox) (Dx - vox) (Dx + vox)]/2;


a     = [];
for i = 1:N
    
    % i-th saccade - position
    %----------------------------------------------------------------------
    pU = DEM{i}.pU.x{1}(1:2,:)*16;
    qU = DEM{i}.qU.x{1}(1:2,:)*16;
    T  = length(pU);
    a  = [a DEM{i}.qU.a{2}];
    
    % eye movements in extrinsic coordinates
    %======================================================================
    subplot(2,2,1)
    for t = 1:T
        image((S + 1)*32), axis image, hold on
        plot(qU(2,t) + Dy,qU(1,t) + Dx,'.g','Markersize',8)
        plot(pU(2,t) + Dy,pU(1,t) + Dx,'.r','Markersize',16)
        plot(pU(2,t) + dj,pU(1,t) + di,'+' ,'Markersize',16)
        
        % show location of image ccentre
        %------------------------------------------------------------------
        plot(Dy,Dx,'+r','Markersize',32);  % axis(ax);

        drawnow, hold off
        
        % save
        %------------------------------------------------------------------
        Me(t + T*(i - 1)) = getframe(gca);
        
    end
    
    % sensory input
    %======================================================================
    subplot(2,2,2)
    for t = 1:T
        
        o   = DEM{i}.pU.x{1}(:,t);
        s   = ADEM_sample_image(STIM.U,o,STIM.R);
        image(s*64), axis image, drawnow
        
        % save
        %----------------------------------------------------------------------
        Mi(t + T*(i - 1)) = getframe(gca);
        
    end
    
    % i-th saccade - percept
    %----------------------------------------------------------------------
    qU = DEM{i}.qU.x{1}(3:end,:);
    
    % percept
    %======================================================================
    subplot(2,2,4)
    for t = 1:T
        
        % hypotheses
        %------------------------------------------------------------------
        h     = spm_softmax(qU(:,t),2);
        H     = 1 + h'*log(h)/log(length(h));
    
        
        % retinotopic predictions
        %------------------------------------------------------------------
        s     = 0;
        for j = 1:length(h)
            s = s + h(j)*spm_read_vols(STIM.S{j});
        end
        image(s*H*64), axis image, drawnow
        
        % save
        %------------------------------------------------------------------
        Mq(t + T*(i - 1)) = getframe(gca);
        
    end
    
end

% set ButtonDownFcn
%--------------------------------------------------------------------------
subplot(2,2,3)
plot(a')
title('Action (EOG)','FontSize',16)
xlabel('time')
axis square

subplot(2,2,1)
set(gca,'Userdata',{Me,16})
set(gca,'ButtonDownFcn','spm_DEM_ButtonDownFcn')
title('saccades (click axis for movie)','FontSize',16)

subplot(2,2,2)
set(gca,'Userdata',{Mi,16})
set(gca,'ButtonDownFcn','spm_DEM_ButtonDownFcn')
title('samples (click axis for movie)','FontSize',16)

subplot(2,2,4)
set(gca,'Userdata',{Mq,16})
set(gca,'ButtonDownFcn','spm_DEM_ButtonDownFcn')
title('percept (click axis for movie)','FontSize',16)

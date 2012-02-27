function spm_dem_occlusion_movie(DEM)
% creates a movie of visual pursuit with occlusion
% FORMAT spm_dem_occlusion_movie(DEM)
%
% DEM - DEM structure from simulations
%
% hidden causes and states
%==========================================================================
% x    - hidden states:
%   x.o(1) - oculomotor angle
%   x.o(2) - oculomotor velocity
%   x.x(1) - target location - extrinsic coordinates
%
% v    - causal states: force on target
%
% g    - sensations:
%   g(1) - oculomotor angle (proprioception)
%   g(2) - oculomotor velocity
%   g(:) - visual input - intrinsic coordinates
%--------------------------------------------------------------------------
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging
 
% Karl Friston
% $Id: spm_dem_occlusion_movie.m 4663 2012-02-27 11:56:23Z karl $
 
 
 
% movie
%--------------------------------------------------------------------------
for i = 1:length(DEM.pU.v{2})
    
    
    % stimulus and tracking in extrinsic coordinates
    %======================================================================
    subplot(2,2,1), cla, axis image ij, hold on
    px  = spm_unvec(DEM.pU.x{1}(:,i),DEM.M(1).x);   % extrinsic target
    qx  = spm_unvec(DEM.qU.x{1}(:,i),DEM.M(1).x);   % conditional estimate
    
    
    % target location
    %----------------------------------------------------------------------
    plot(0,0,'+','MarkerSize',32,'Color',[1 1 1]*(1 - 1/8))
    c   = max(DEM.pU.v{1}(3:end,i));
    c   = 1 - min(max(0,c),1);
    plot(px.x(1),0,'.','MarkerSize',32,'Color',[1 c c])
    c   = max(DEM.qU.v{1}(3:end,i));
    c   = 1 - min(max(0,c),1);
    plot(qx.x(1),0,'.','MarkerSize',32,'Color',[1 c c])
    plot(px.o(1),0,'+','MarkerSize',32,'Color',[1 .4 .4])
    axis image, grid on, box off
    axis([-1 1 -1 1]*(1 + 1/2))
    drawnow
    
    % save
    %----------------------------------------------------------------------
    Me(i) = getframe(gca);
    
end
 
% set ButtonDownFcn
%--------------------------------------------------------------------------
subplot(2,2,1)
set(gca,'Userdata',{Me,16})
set(gca,'ButtonDownFcn','spm_DEM_ButtonDownFcn')
title('Extrinsic (left click for movie)','FontSize',16)
xlabel('displacement (cm)')
hold off
 
 
% stimulus and tracking in intrinsic coordinates
%==========================================================================
subplot(2,2,2)
dt  = 16;                                            % time step (ms)
dr  = 2;                                             % degrees
x   = size(DEM.pU.v{1},1) - 2;
t   = size(DEM.pU.v{1},1);
imagesc([1 t]*dt,[-x x]/2,DEM.pU.v{1}(3:end,:))
title('Intrinsic (retinal input)','FontSize',16)
ylabel('displacement (degrees)')
xlabel('time (ms)')
axis square
 
subplot(2,2,3)
t    = size(DEM.pU.v{1},2);
t    = (1:t)*dt;
do   = DEM.pU.x{1}(1,:)*dr;
dx   = DEM.pU.x{1}(3,:)*dr;
plot(t,do,'-.k',t,dx,'k')
title('target and oculomotor angles','FontSize',16)
xlabel('time (ms)')
ylabel('displacement (degrees)')
axis square
spm_axis tight
legend('eye', 'target');
 
subplot(2,2,4)
occ  = max(DEM.pU.v{1}(3:end,:)) < 1/32;
dodt = diff(DEM.pU.x{1}(1,:))*1000*dr/dt;
dxdt = diff(DEM.pU.x{1}(3,:))*1000*dr/dt;
plot(t(1:end - 1),dodt,'-.k',t(1:end - 1),dxdt,'k',t,occ*128 - 64,':k')
title('target and oculomotor velocities','FontSize',16)
xlabel('time (ms)')
ylabel('velocity (degrees per second)')
axis square
spm_axis tight 



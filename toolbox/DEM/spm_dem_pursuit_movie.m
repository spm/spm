function spm_dem_pursuit_movie(DEM,c)
% creates a movie of visual prusuit in extrinsic and intrinsic coordinates
% FORMAT spm_dem_pursuit_movie(DEM)
%
% DEM - DEM structure from reaching simulations
%
% hidden causes and states
%--------------------------------------------------------------------------
% x    - hidden states:
%   x.o(1) - oculomotor angle
%   x.o(2) - oculomotor angle
%   x.x(1) - target location (visual) - extrinsic coordinates (Cartesian)
%   x.x(2) - target location (visual) - extrinsic coordinates (Cartesian)
%   x.a(:) - attractor (SHC) states
%
% v    - causal states
%   v(1) - not used
%
% g    - sensations:
%   g(1) - oculomotor angle (proprioception)
%   g(2) - oculomotor angle (proprioception)
%   g(3) - target location (visual) - intrinsic coordinates (polar)
%   g(4) - target location (visual) - intrinsic coordinates (polar)
%__________________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging
 
% Karl Friston
% $Id: spm_dem_pursuit_movie.m 4322 2011-05-04 15:28:08Z karl $


% check subplot specifier
%--------------------------------------------------------------------------
try, c; catch, c = 0; end
 
% movie
%--------------------------------------------------------------------------
N      = length(DEM.pU.v{2});
p      = -pi:pi/64:pi;
circle = [sin(p);cos(p)];
 
for i = 1:N
 
    
    % stimulus and tracking in extrinsic coordinates
    %======================================================================
    subplot(2,2,1 + c), cla, axis image ij, hold on
    px  = spm_unvec(DEM.pU.x{1}(:,i),DEM.M(1).x);   % extrinsic target
    qx  = spm_unvec(DEM.qU.x{1}(:,i),DEM.M(1).x);   % conditional estimate
    
 
    % target location
    %----------------------------------------------------------------------
    plot(0,0,'+','MarkerSize',32,'Color',[1 1 1]*(1 - 1/8))
    plot(px.x(1),px.x(2),'.','MarkerSize',32,'Color',[1 .0 .0])
    plot(qx.x(1),px.x(2),'.','MarkerSize',32,'Color',[1 .4 .4])
    plot(px.o(1),px.o(2),'+','MarkerSize',32,'Color',[1 .4 .4])
    axis image, grid on, box off
    axis([-1 1 -1 1]*(1 + 1/8))
    drawnow
    
    % save
    %----------------------------------------------------------------------
    Me(i) = getframe(gca);
    
    
    % stimulus and tracking in intrinsic coordinates
    %======================================================================
    subplot(2,2,2 + c), cla, hold on
    pv  = DEM.pU.v{1}(:,i);   % intrinsic stimulus
    
    % target location
    %----------------------------------------------------------------------
    for j = 1:3
       plot(circle(1,:)*j/2,circle(2,:)*j/2,':k');
    end
    plot(p,p*0,':k');plot(p*0,p,':k');
    plot(pv(3),pv(4),'.','MarkerSize',32,'Color',[1 .0 .0]);     
    plot(-tan(pv(1)),-tan(pv(2)),'+','MarkerSize',32,'Color',[1 1 1]*(1 - 1/8));
    plot(0,0,'+','MarkerSize',32,'Color',[1 .4 .4]);
    axis image, box off
    axis([-1 1 -1 1]*pi/2)
    drawnow
    
    % save
    %----------------------------------------------------------------------
    Mi(i) = getframe(gca);
 
end
 
 
% set ButtonDownFcn
%--------------------------------------------------------------------------
subplot(2,2,1 + c)
set(gca,'Userdata',{Me,16})
set(gca,'ButtonDownFcn','spm_DEM_ButtonDownFcn')
title('Extrinsic (left click for movie)','FontSize',16)
xlabel('displacement (cm)')
hold off
 
subplot(2,2,2 + c)
set(gca,'Userdata',{Mi,16})
set(gca,'ButtonDownFcn','spm_DEM_ButtonDownFcn')
title('Intrinsic (left click for movie)','FontSize',16)
xlabel('displacement (radians)')
hold off

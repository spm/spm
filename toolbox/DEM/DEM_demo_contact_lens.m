% This demo illustrates tracking under the contact lens problem:
% The contact lens refers to the non-Gaussian uncertainty induced by
% nonlinear measurements. Here it is illustrated in terms of tracking the
% motion of a target in Cartesian coordinates, given the distance to target
% (range) and direction as measurements. The problem is to accumulate
% information over time about the target location under random fluctuations
% on the velocity (technically this is a constant acceleration model).
% Comparative evaluations are made with Extended Kalman filtering.
%
% See: X. Tian, Y. Bar-Shalom, Coordinate Conversion and Tracking for 
% Very Long Range Radars. IEEE Transactions on Aerospace and Electronic
% Systems, AES-45(3):1073–1088, July 2009.
%__________________________________________________________________________
% Copyright (C) 2011 Wellcome Trust Centre for Neuroimaging
 
% Karl Friston
% $Id: DEM_demo_contact_lens.m 4230 2011-03-07 20:58:38Z karl $
 
 
% non-linear generative model
%==========================================================================
clear
 
% The problem: states = [x(1) x(2)]; causes = [v(1) v(2)]
%--------------------------------------------------------------------------
x        = [1e4 1e4]';         % initial states: location = 10km,  10km
v        = [-8   -2]';         %                 velocity = -8m/s, -2m/s
V        = [1e-2  4];          % observation precision (inverse variance) 
                               % standard deviation of range: 10m
                               % standard deviation of angle: .5 mrad
s        = 128;                % precision of fluctuations in motion
                               % standard deviation of velocity: 0.08 m/s^2
                               
 % preliminaries
%--------------------------------------------------------------------------
N        = 128;                % length of sequence
t        = 1:N;                % time (seconds)
M(1).E.s = 1/sqrt(2*s);        % smoothness of fluctuations
M(1).E.n = 4;                  % order of generalised coordinates
 
 
% level 1 (simple slow motion of a distant radar target)
%--------------------------------------------------------------------------
g      = '[sqrt(x(1)^2 + x(2)^2); 1000*atan(x(2)/x(1))]';
M(1).f = inline('v','x','v','P');
M(1).g = inline( g ,'x','v','P');
M(1).x = x;
M(1).V = diag(V);
M(1).W = s^2;
 
% level 2 (fluctuating velocity)
%--------------------------------------------------------------------------
M(2).v = v;
M(2).V = s;
 
% create data
%==========================================================================
DEM    = spm_DEM_generate(M,N);
 
spm_figure('Getwin','DEM');
spm_DEM_qU(DEM.pU)
 
 
% Comparative inversions (variants of generalised filtering)
%==========================================================================
DEM.M(1).x = [1e4; 2e3];
DEM.M(2).v = [0; 0];
DEM.M(1).V = exp(8);
DEM.M(1).W = exp(4);
DEM.M(2).V = exp(0);
 
 
% DEM
%--------------------------------------------------------------------------
DEM        = spm_DEM(DEM);
 
spm_figure('Getwin','DEM');
spm_DEM_qU(DEM.qU,DEM.pU)
 
 
% Extened Kalman filtering: states = [x(1) x(2) v(1) v(2)];
%--------------------------------------------------------------------------
f      = '[x(3); x(4); 0; 0]';
g      = '[sqrt(x(1)^2 + x(2)^2); 1000*atan(x(2)/x(1))]';
G(1).f = inline(f,'x','v','P');
G(1).g = inline(g,'x','v','P');
G(1).x = [1e4; 2e3; 0; 0];
G(1).V = diag(V);
G(1).W = diag([s^2 s^2 s s]);
 
% EKF  = spm_ekf(DEM.M,DEM.Y);
EKF    = spm_ekf(G,DEM.Y);
 
% show prediction error
%==========================================================================
spm_figure('Getwin','Graphics');
 
D(1,:)  = sqrt(sum((DEM.pU.x{1}([1 2],:) - DEM.qU.x{1}([1 2],:)).^2));
D(2,:)  = sqrt(sum((DEM.pU.x{1}([1 2],:) - EKF([1 2],:)).^2));
 
% plot
%--------------------------------------------------------------------------
subplot(2,1,1)
plot(t,log(D(:,:)))
title('log(location error (m))','Fontsize',16)
xlabel('time (secs)')
ylabel('log distance from target')
axis square
legend('DEM','EKF')
drawnow
 
subplot(2,1,2)
plot(DEM.pU.x{1}(1,:),DEM.pU.x{1}(2,:)     ),hold on
plot(DEM.qU.x{1}(1,:),DEM.qU.x{1}(2,:),'-.'),hold on
plot(EKF(1,2:end),EKF(2,2:end),':'),hold off
title('location (m)','Fontsize',16)
axis([90 105 95 105]*1e2)
xlabel('x')
ylabel('y')
axis square

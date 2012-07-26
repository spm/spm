function DEM_demo_contact_lens
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
% $Id: DEM_demo_contact_lens.m 4804 2012-07-26 13:14:18Z karl $
 
 
% non-linear generative model
%==========================================================================
 
% The problem: states = [x(1) x(2)]; causes = [v(1) v(2)]
%--------------------------------------------------------------------------
x        = [1e4 1e4 -8 -2]';   % initial states: location = 10km,  10km
                               %                 velocity = -8m/s, -2m/s
V        = [1e-2 4];           % observation precision (inverse variance) 
                               % standard deviation of range: 1/sqrt(V(1)) = 10m
                               % standard deviation of angle: 1/sqrt(V(2)) = .5 mrad
s        = 1;                  % smoothness of fluctuations
w        = 8;                  % precision of fluctuations in motion
W        = [128 128 w w];      % standard deviation of velocity: 1/sqrt(w) = .3536 m/s^2
                               
% preliminaries
%--------------------------------------------------------------------------
N        = 256;                % length of sequence
t        = 1:N;                % time (seconds)
M(1).E.s = s;                  % smoothness of fluctuations
M(1).E.n = 4;                  % order of generalised coordinates
M(1).E.K = 128;                % rate of generalized gradient descent
 
% Level 1: Hidden states = [x(1) x(2) v(1) v(2)];
%--------------------------------------------------------------------------
f        = '[x(3); x(4); 0; 0]';
g        = '[sqrt(x(1)^2 + x(2)^2); 1000*atan(x(2)/x(1))]';
M(1).f   = inline(f,'x','v','P');
M(1).g   = inline(g,'x','v','P');
M(1).x   = x;
M(1).V   = diag(V);            % precision of observation noise
M(1).W   = diag(W);            % precision of fluctuations on hidden states
 
   
% create data
%==========================================================================
DEM    = spm_DEM_generate(M,N);
 

% Comparative inversions (variants of generalised filtering)
%==========================================================================

% reset initial position and bearing
%--------------------------------------------------------------------------
DEM.M(1).x = [1e3; 1e3; 0; 0];

% DEM and EKF
%--------------------------------------------------------------------------
DEM     = spm_DEM(DEM);
[EKF S] = spm_ekf(DEM.M,DEM.Y);
 
spm_figure('Getwin','DEM');
spm_DEM_qU(DEM.qU,DEM.pU)
 
 
% show prediction errors
%==========================================================================
spm_figure('GetWin','Figure 1');
 
D(1,:) = sqrt(sum((DEM.pU.x{1}([1 2],:) - DEM.qU.x{1}([1 2],:)).^2));
D(2,:) = sqrt(sum((DEM.pU.x{1}([1 2],:) - EKF([1 2],:)).^2));
D(3,:) = sqrt(sum((DEM.pU.x{1}([3 4],:) - DEM.qU.x{1}([3 4],:)).^2));
D(4,:) = sqrt(sum((DEM.pU.x{1}([3 4],:) - EKF([3 4],:)).^2));

for i = 1:N
    E         = DEM.pU.x{1}(:,i) - DEM.qU.x{1}(:,i);
    NEES(1,i) = E'*spm_inv(DEM.qU.S{i})*E;
    E         = DEM.pU.x{1}(:,i) - EKF(:,i);
    NEES(2,i) = E'*spm_inv(S{i})*E;
end


% plot errors
%--------------------------------------------------------------------------
subplot(2,2,1)
plot(t,log(D(1:2,:)))
title('log(location error (m))','Fontsize',16)
xlabel('time (secs)')
ylabel('log distance from target')
axis square
legend('DEM','EKF')
 
subplot(2,2,2)
plot(t,log(D(3:4,:)))
title('log(speed error (m/s))','Fontsize',16)
xlabel('time (secs)')
ylabel('log distance from target')
axis square
legend('DEM','EKF')


% plot trajectories
%--------------------------------------------------------------------------
subplot(2,2,3)
plot(DEM.pU.x{1}(1,:),DEM.pU.x{1}(2,:)     ),hold on
plot(DEM.qU.x{1}(1,:),DEM.qU.x{1}(2,:),'-.'),hold on
plot(EKF(1,:),EKF(2,:),':'),hold off
title('location (m)','Fontsize',16)
axis([8 12 8 12]*1000)
xlabel('x')
ylabel('y')
axis square
legend('true','DEM','EKF')

% plot NEES
%--------------------------------------------------------------------------
subplot(2,2,4)
plot(t,log(NEES))
title('log(NEES)','Fontsize',16)
xlabel('time (secs)')
axis square
legend('DEM','EKF')
drawnow


return

% prediction errors as a function of K
%==========================================================================
K     = -0:1/2:12;
for i = 1:length(K)
    
    % DEM and EKF (linear)
    %----------------------------------------------------------------------
    DEM.M(1).E.K = exp(K(i));
    DEM          = spm_DEM(DEM);
    
    % mean square error
    %----------------------------------------------------------------------
    kse(i,1) = mean(mean((DEM.pU.x{1} - DEM.qU.x{1}).^2));
    
end

% show prediction errors
%==========================================================================
spm_DEM_qU(DEM.qU,DEM.pU)
 
% plot trajectories
%--------------------------------------------------------------------------
subplot(2,1,2)
plot(K,log(kse))
title('log(MSE)','Fontsize',16)
xlabel('log(K)','Fontsize',12)
axis square



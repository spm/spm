function DEM_path_integrals
%--------------------------------------------------------------------------
% Illustration of approximations to path integrals. This routine generates
% a path from dynamics whose Fokker Planck solution corresponds to a
% Gaussian with a given (diagonal) precision. It then samples random
% segments (after scaling and smoothing) and evaluates their action. This
% evaluation is in terms of the sum of squares residuals between realised
% and predicted flow plus path dependent and path-independent terms based
% upon the surprisal associated with the solution of the Fokker Planck
% equation. The point being made here is that the terms based upon the
% surprisal (cyan dots) upper bound the action (blue dots).
%__________________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% Karl Friston
% $Id: DEM_path_integrals.m 7447 2018-10-13 15:32:10Z karl $


% set up
%--------------------------------------------------------------------------
spm_figure('GetWin','Figure 1'); clf

dt    = 2^(-8);                            % time step
N     = 2^16;                              % number of time steps
G     = 256;                               % random fluctuations
P     = inv([8,0;0,1]);                    % precision of NESS density
Q     = [0 -1; 1 0]*G/4;                   % solenoidal flow
J     = @(x,P) (-1/2)*sum(P*x.^2);         % log NESS density
f     = @(x,G,Q,P) (Q - G*eye(2))*P*x;     % flow
xt    = [0;0];                             % initial state

% Integrate a path
%--------------------------------------------------------------------------
% more accurate integration scheme:
% dx    = spm_sde_dx(Q - G*eye(2),sqrt(G)*eye(2),f(xt,G,Q),dt);
%--------------------------------------------------------------------------
for t = 1:N
    x(:,t) = xt;
    dx     = f(xt,G,Q,P)*dt + sqrt(2*G*dt)*randn(2,1);
    xt     = xt + dx;
end

% check sample covariance corresponds to the NESS density
%--------------------------------------------------------------------------
disp('numerical covariance')
disp(cov(x'))


% shows states visited
%--------------------------------------------------------------------------
subplot(2,2,1),  plot(x(1,:),x(2,:),'.','MarkerSize',1)
title('States visited','FontSize',16)
xlabel('state space'), ylabel('state space')
axis square, axis xy


n     = 256;                            % length of each path
s     = 512;                            % number of paths
for i = 1:s
    
    % sample a smooth path
    %----------------------------------------------------------------------
    x0   = ceil(rand*(N - n));
    xt   = x(:,x0 + (1:n))*randn*2;
    xt   = spm_conv(xt,0,n/8);
    
    % observed unpredicted flow
    %----------------------------------------------------------------------
    dx   = diff(xt,1,2)/dt;
    qx   = f(xt(:,1:end - 1),G,Q,P);
    
    % evaluate action
    %----------------------------------------------------------------------
    A(i,1) = -sum(sum((dx - qx).^2)/(2*G))/2*dt;  % path integral form
    L(i,1) = -sum(sum(dx.^2)/(2*G))/2*dt;         % Kinetic
    L(i,2) = (J(xt(:,end),P) - J(xt(:,1),P))/2;   % Path-independent
    L(i,3) = G*sum(J(sqrt(P)*xt,P))/2*dt;         % Path-dependent
    
    
    % illustrate a subset of sampled paths
    %----------------------------------------------------------------------
    if i < 32
        col  = spm_softmax(randn(3,1).^2);
        subplot(2,2,3),  plot(xt(1,:),xt(2,:),'Color',col),hold on
    end
    
end

% actions of paths
%----------------------------------------------------------------------
title('Sample paths','FontSize',16)
xlabel('state space'), ylabel('state space')
axis square, axis xy

subplot(2,2,2),  plot(A,sum(L,2),'.',A,sum(L(:,[2 3]),2),'c.',A,A,'-')
title('Estimated action','FontSize',16)
xlabel('action'), ylabel('Estimates')
axis square, axis xy, hold off
legend({'all components','surprisal terms','true value'},'Location','SouthEast')

subplot(2,2,4),  plot(A,L,'.')
title('Components','FontSize',16)
xlabel('action paths'), ylabel('action components')
axis square, axis xy
legend({'kinetic','path-independent','path-dependent'},'Location','SouthEast')



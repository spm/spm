function [Ey,Cy,Vy] = spm_NESS_forecastness(DEM)
% Analytic predictive density under the Laplace assumption
% FORMAT [Ey,Cy,Vy] = spm_NESS_forecastness(DEM)
%--------------------------------------------------------------------------
% DEM - inverted DEM structure
%
% Ey  - posterior expectation
% Cy  - posterior covariances
% Vy  - posterior covariances
%
% This routine integrates the Fokker Planck equation under a Laplace
% (Gaussian) approximation to the evolving density over states. It uses the
% Hessian at the current point in state space to approximate the (inverse)
% covariance and assumes an initial density with a high precision. This
% application involves some housekeeping because the parameters of the flow
% pertain to a first-order approximation to the coefficients of a Helmholtz
% decomposition (with the surprisal specified in terms of kernels), whereas
% the Gaussian approximation requires a second-order polynomial expansion
% (i.e., K = 3).
%__________________________________________________________________________

% state-space for (Laplace) solution 
%==========================================================================

% current time (weeks)
%--------------------------------------------------------------------------
T     = size(DEM.Y,2);               % current time
nT    = size(DEM.U,2);               % number of future time points

% scale responses, Y to obtain latent states, X
%--------------------------------------------------------------------------
P     = DEM.M(1).pE;                 % flow parameters (order K)
g     = DEM.M(1).g;                  % observer g
gy    = DEM.M(1).gy;                 % inverse g
W     = DEM.M(1).W;                  % precision
f     = DEM.M(1).f;                  % flow 
u     = DEM.U(:,end);                % current cause 
x     = gy(DEM.Y(:,end)',u',P);      % current state
[n,m] = size(P.scale);               % number of states and responses      

% parameters of equilibrium density (Ep: Laplace approximation)
%==========================================================================

% surprisal coefficients for (Gaussian) NESS 
%--------------------------------------------------------------------------
E     = f(x,u,pE,[],'E');            % expectation (NESS)
H     = f(E,u,pE,[],'H');            % Hessian (NESS) at expectation
C     = inv(H + eye(n,n)/4);         % upper bound covariance (2 s.d.)
Ep.Sp = spm_ness_cond_inv(E(:),C);   % Laplace approximation

% solenoidal coefficients (for a K = 3 basis set)
%--------------------------------------------------------------------------
nb    = numel(Ep.Sp);
Ep.Qp = [];
for i = 1:n
    for j = i:n
        for k = 1:nb
            Ep.Qp(end + 1) = 0;
        end
    end
end

% amplitude of random fluctuations
%--------------------------------------------------------------------------
Ep.G  = inv(W)/2;

% Suprisal coeffcients of initial (Gaussian) density
%==========================================================================
x     = x(:);
Sp    = spm_ness_cond_inv(x,Ep.G);

% Run system forwards in time to nonequilibrium steady-state density
%--------------------------------------------------------------------------
U     = DEM.B;
Ey    = zeros(m,nT);                 % expectation of density
Vy    = ones(m,nT);                  % variance of density
Cy    = cell(1,nT);                  % covariance of density
dt    = 1;
for t = 1:nT

    % momments of Laplace approximation
    %----------------------------------------------------------------------
    [~,C]   = spm_ness_cond(n,3,Sp);
    Sp      = spm_ness_cond_inv(x,C);
    [m,C]   = spm_ness_cond(n,3,Sp);

    Ey(:,t) = g(m,u,P);
    Cy{t}   = P.scale'*C*P.scale;
    Vy(:,t) = diag(Cy{t});

    % update surprisal coefficients of current density
    %----------------------------------------------------------------------
    [dfds,fs] = spm_diff(@spm_NESS_ds,Sp,P,U,1);
    fS        = U.b\fs;
    dfdS      = U.b\dfds;
    dS        = spm_dx(dfdS,fS,dt);
    Sp        = Sp + dS;

    % update path of least action
    %----------------------------------------------------------------------
    for i = 1:8
        [dfdx,fx] = spm_diff(f,x,u,P,1);
        dx        = spm_dx(dfdx,fx,dt/8);
        x         = x + dx;
    end

end


% plot the past and forecast
%==========================================================================

% time points
%--------------------------------------------------------------------------
s     = (1:T);
t     = (1:nT) + T - 1;              % future time points

% Density dynamics
%--------------------------------------------------------------------------
for i = 1:n
    subplot(6,2,i),  hold on, set(gca,'ColorOrderIndex',i)
    spm_plot_ci(Ey(i,:),Vy(i,:),t,[],'plot')
    hold on, set(gca,'ColorOrderIndex',i),
    plot(s,DEM.Y(i,:),'LineWidth',2), plot(get(gca,'XLim'),[0,0],':k')
    title('predictive density','FontSize',14)
end

return


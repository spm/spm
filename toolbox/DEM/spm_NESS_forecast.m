function [Ey,Cy,Vy] = spm_NESS_forecast(DEM)
% Analytic predictive density under the Laplace assumption
% FORMAT [Ey,Cy,Vy] = spm_NESS_forecast(DEM)
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
f     = DEM.M(1).f;                  % flow 
x     = DEM.M(1).x;                  % current state
u     = DEM.M(2).v;                  % current cause
[n,m] = size(P.scale);               % number of states and responses

% amplitude of random fluctuations
%--------------------------------------------------------------------------
G     = f(x,u,P,[],'G');

% moments of initial (Gaussian) density
%==========================================================================
x     = x(:);
C     = G;

% Run system forwards in time to nonequilibrium steady-state density
%--------------------------------------------------------------------------
Ey    = zeros(m,1);                  % expectation of density
Vy    = ones(m,1);                   % variance of density
Cy    = {eye(n,n)};                  % covariance of density
dt    = 4;                           % integration steps
for t = 1:nT

    % momments of Laplace approximation
    %----------------------------------------------------------------------
    u       = DEM.U(:,t);            % causes 
    Ey(:,t) = g(x,u,P);              % prediction: expectation
    Cy{t}   = P.scale'*C*P.scale;    % prediction: covariance
    Vy(:,t) = diag(Cy{t});           % prediction: variance

    % density dynamics
    %----------------------------------------------------------------------
    dfdxx = spm_diff(f,x,u,P,[1,1]);
    for i = 1:dt

        % mean
        %------------------------------------------------------------------
        [dfdx,fx] = spm_diff(f,x,u,P,1);
        for j = 1:numel(x)
            fx(j) = fx(j) + trace(C*dfdxx{j})/2;
        end
        dx = spm_dx(dfdx,fx,1/dt);
        x  = x + dx;

        % check for numerical divergence
        %------------------------------------------------------------------
        if any(isnan(x)), break, end

        % covariance
        %------------------------------------------------------------------
        dC = dfdx*C + C*dfdx' + G + G';
        C  = C + dC/dt;
        C  = sqrtm(C'*C);

    end

    % check for numerical divergence
    %----------------------------------------------------------------------
    if any(isnan(x)) || any(abs(x) > 16)
        nT = t; break
    end

end

% plot the past and forecast
%==========================================================================
if isfield(DEM,'nograph'), return, end

% time points
%--------------------------------------------------------------------------
s     = (1:nT) + T - nT;
t     = (1:nT) + T - 1;

% Density dynamics
%--------------------------------------------------------------------------
for i = 1:n
    subplot(6,2,i),  hold on, set(gca,'ColorOrderIndex',i)
    spm_plot_ci(Ey(i,:),Vy(i,:),t,[],'plot')
    hold on, set(gca,'ColorOrderIndex',i),
    plot(s,DEM.Y(i,s),'LineWidth',2), plot(get(gca,'XLim'),[0,0],':k')
    title('predictive density','FontSize',14)
end

return


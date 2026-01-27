function [Ey,Cy,Vy] = spm_NESS_forecast(DEM,nT)
% Analytic predictive density under the Laplace assumption
% FORMAT [Ey,Cy,Vy] = spm_NESS_forecast(DEM,nT)
%--------------------------------------------------------------------------
% DEM - inverted DEM structure
% nT  - number of future time points
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

% scale responses, Y to obtain latent states, X
%--------------------------------------------------------------------------
scale = diag(DEM.M(1).pE.scale);     % scaling
x     = DEM.Y(:,end)'/scale;         % current state
u     = DEM.U(:,end);                % current cause 
pE    = DEM.M(1).pE;                 % flow parameters (order K) 
W     = DEM.M(1).W;                  % precision
n     = numel(x);                    % number of states

% parameters of equilibrium density (Ep: Laplace approximation)
%==========================================================================

% surprisal coefficients for (Gaussian) NESS 
%--------------------------------------------------------------------------
E     = DEM.M(1).f(x,u,pE,[],'E');   % expectation (NESS)
H     = DEM.M(1).f(E,u,pE,[],'H');   % Hessian (NESS) at expectation
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
Sp    = spm_ness_cond_inv(x(:),Ep.G);

% Run system forwards in time to nonequilibrium steady-state density
%--------------------------------------------------------------------------
U     = DEM.B;
Ey    = zeros(n,nT);                 % expectation of density
Vy    = ones(n,nT);                  % variance of density
Cy    = cell(1,nT);                  % covariance of density
dt    = 1;
for t = 1:nT

    % momments of Laplace approximation
    %----------------------------------------------------------------------
    [m,C]   = spm_ness_cond(n,3,Sp);
    Ey(:,t) = scale*m;
    Cy{t}   = scale*C*scale;
    Vy(:,t) = diag(Cy{t});

    % update surprisal coefficients of current density
    %----------------------------------------------------------------------
    ds  = spm_NESS_ds(Sp,Ep,U);
    dds = spm_diff(@spm_NESS_ds,Sp,Ep,U,1);
    dS  = U.b\ds;
    ddS = U.b\dds;
    dS  = spm_dx(ddS,dS,dt);
    Sp  = Sp + dS;

end


% plot the past and forecast
%==========================================================================

% current time (weeks)
%--------------------------------------------------------------------------
if isfield(DEM,'T'), T = DEM.T;  else, T = size(DEM.Y,2);  end

% time points
%--------------------------------------------------------------------------
r     = size(DEM.Y,2);
r     = (1:r)  + T - r;
t     = (1:nT) + T - 1;              % future time points

% Density dynamics
%--------------------------------------------------------------------------
for i = 1:n
    subplot(6,2,i),  hold on, set(gca,'ColorOrderIndex',i)
    spm_plot_ci(Ey(i,:),Vy(i,:),t,[],'plot')
    hold on, set(gca,'ColorOrderIndex',i),
    plot(r,DEM.Y(i,:),'LineWidth',2), plot(get(gca,'XLim'),[0,0],':k')
    title('predictive density','FontSize',14)
end

return


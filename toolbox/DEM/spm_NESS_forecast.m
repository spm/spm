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
X     = DEM.Y'/scale;                % past states
x     = X(end,:);                    % current state                
pE    = DEM.M(1).pE;                 % flow parameters (order K) 
W     = DEM.M(1).W;                  % precision
n     = numel(x);                    % number od states

% parameters of equilibrium density (Ep: Laplace approximation)
%==========================================================================

% surprisal coefficients for (Gaussian) NESS 
%--------------------------------------------------------------------------
M.K   = 2;
M.L   = 2;
M.W   = W;
[~,~,~,~,H,~,E] = spm_NESS_gen_lap(pE,M,x);
C     = inv(H + eye(n,n)/4);
Ep.Sp = spm_ness_cond_inv(E(:),C);

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
Sp    = spm_ness_cond_inv(x(:),C/4);

% Run system forwards in time to nonequilibrium steady-state density
%--------------------------------------------------------------------------
U     = DEM.B;
dt    = 1/32;
N     = nT/dt;
Ey    = zeros(n,N);               % expectation of density
Vy    = ones(n,N);                % variance of density
Cy    = cell(1,N);                % covariance of density
for t = 1:N

    % momments of Laplace approximation
    %----------------------------------------------------------------------
    [m,C]   = spm_ness_cond(n,3,Sp);
    Ey(:,t) = scale*m;
    Cy{t}   = scale*C*scale;
    Vy(:,t) = diag(Cy{t});

    % update surprisal coefficients of current density
    %----------------------------------------------------------------------
    dS = spm_NESS_ds(Sp,Ep,U);
    dS = U.b\dS;
    Sp = Sp + dS*dt;

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
t     = (1:nT) + T - 1;            % future time points
i     = ((1:nT) - 1)/dt + 1;

Ey    = Ey(:,i);
Cy    = Cy(i);
Vy    = Vy(:,i);

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

% Notes
%==========================================================================

% update mean and (positive definite) covariance
%--------------------------------------------------------------------------
[dm,dC] = spm_ness_cond(n,3,Sp + dS/1e6);

dm      = (dm - m)*1e6;
m       = m + dm*dt;
[e,v]   = eig(C);
v       = diag(v);
dv      = (diag(e'*dC*e)./v - 1)*1e6;
v       = log(v) + dv*dt;
C       = e*diag(exp(v))*e';
Sp      = spm_ness_cond_inv(m,C);

% method of moments: based on sampled surprisal
%--------------------------------------------------------------------------
X     = DEM.B.X;
for i = 1:size(X,1)
    [~,s]  = spm_NESS_gen_lap(pE,M,X(i,:) + E);
    S(i,1) = s;
end
p = spm_softmax(-S);
e = p'*X;
c = 0;
for i = 1:size(X,1)
    s = X(i,:) - e;
    c = c + p(i)*(s'*s);
end

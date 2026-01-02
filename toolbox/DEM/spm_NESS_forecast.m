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

% parameters of equilibrium density (Ep: Laplace approximation)
%==========================================================================

% surprisal coefficients for (Gaussian) NESS 
%--------------------------------------------------------------------------
M.K   = 2;
M.L   = 2;
M.W   = W;
[~,~,~,~,H,~,E] = spm_NESS_gen_lap(pE,M,x);
C      = inv(H);
[Sp,o] = spm_ness_cond_inv(E(:),C);
Ep.Sp  = Sp;

% solenoidal coefficients (for a K = 3 basis set)
%--------------------------------------------------------------------------
n     = numel(x);
K     = 2;
p     = (1:K) - 1;
for i = 2:n
    p = repmat(p,1,K);
    p = [p; kron((1:K) - 1,ones(1,K^(i - 1)))];
end
k     = sum(p) < K;
p     = p(:,k);
no    = size(o,2);
nb    = size(p,2);
Qp    = zeros(no,n,n);
for k = 1:nb
    r(k) = find(ismember(o',p(:,k)','rows'));
end
d     = 1;
for i = 1:n
    for j = i:n
        for k = 1:nb
            Qp(r(k),i,j) = pE.Qp(d);
            d = d + 1;
        end
    end
end
Ep.Qp = [];
for i = 1:n
    for j = i:n
        for k = 1:no
            Ep.Qp(end + 1) = Qp(k,i,j);
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

    % time-derivative of surprisal coefficients
    %----------------------------------------------------------------------
    dS = spm_NESS_ds(Sp,Ep,U);
    dS = U.b\dS;

    % update mean and (positive definite) covariance
    %----------------------------------------------------------------------
    [dm,dC] = spm_ness_cond(n,3,Sp + dS*dt);

    dm = dm - m;
    dC = dC - C;
    m  = m + dm;
    C  = expm(logm(C) + C\dC);
    Sp = spm_ness_cond_inv(m,C);

    % update surprisal coefficients
    %----------------------------------------------------------------------
    % [~,~,~,~,H,~,E] = spm_NESS_gen_lap(pE,M,m);
    % Ep.Sp           = spm_ness_cond_inv(E(:),inv(H));
    % fprintf('Fokker–Planck: %i of %i\n',t,N)

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
xLim  = [-4 1]*nT + T;
for i = 1:n
    subplot(6,2,i),  hold on, set(gca,'ColorOrderIndex',i)
    spm_plot_ci(Ey(i,:),Vy(i,:),t,[],'plot')
    hold on, set(gca,'ColorOrderIndex',i),
    plot(r,DEM.Y(i,:),'LineWidth',2), plot(xLim,[0,0],':k')
    title('predictive density','FontSize',14)
    set(gca,'XLim',xLim)
end


return
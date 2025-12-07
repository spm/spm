function [Ey,Cy] = spm_NESS_forecasting(DEM,nT)
% priors for a NESS genertive model
% FORMAT [Ey,Cy] = spm_NESS_forecasting(DEM,nT)
%--------------------------------------------------------------------------
% DEM - inverted DEM structure
% nT  - number of future time points
%
% Ey  - posterior expectation
% Cy  - posterior covariances
%
% This routine takes a DEM structure that has been inverted given some data
% Y and generates trajectories into the future over nT time points. The
% sample density is then recorded in terms of the first Ey and second order
% Cy moments. The ensuing trajectories are plotted; along with the
% preceding timeseries data (and posteriors over latent states).
%__________________________________________________________________________


% number of states and samples
%--------------------------------------------------------------------------
N     = 32;                        % number of paths
n     = numel(DEM.M(1).x);         % number of states
T     = size(DEM.Y,2);             % number of time points
t     = (1:nT) + T;                % future time points

% calibrate
%--------------------------------------------------------------------------
DEM.M(1).W  = DEM.M(1).W/(4^2);

% plot the past
%--------------------------------------------------------------------------
spm_DEM_qU(DEM.qU)
subplot(2,2,1),      set(gca,'ColorOrderIndex',1), hold on
plot(DEM.Y','.'),    set(gca,'ColorOrderIndex',1)


% future projections
%--------------------------------------------------------------------------
Ep    = DEM.qP.P(1);               % posterior expectation of parameters
pY    = zeros(n,nT,N);             % predictive realisations
for i = 1:32

    % generate trajectory from posterior initial conditions
    %----------------------------------------------------------------------
    Cx         = full(DEM.qU.S{end});
    x0         = DEM.qU.x{1}(:,end) + sqrtm(Cx)*randn(n,1);
    DEM.M(1).x = x0;
    FEM        = spm_DEM_generate(DEM.M,nT,Ep);
    pY(:,:,i)  = full(FEM.Y);

    % projections
    %----------------------------------------------------------------------
    subplot(2,2,1),          set(gca,'ColorOrderIndex',1), hold on
    plot(t,FEM.Y',':'),      set(gca,'ColorOrderIndex',1)
    set(gca,'XLim',[(T - 2*nT),t(end)])

    % hidden states
    %----------------------------------------------------------------------
    subplot(2,2,2),          set(gca,'ColorOrderIndex',1), hold on
    plot(t,FEM.pU.x{1},':'), set(gca,'ColorOrderIndex',1)
    set(gca,'XLim',[(T - 2*nT),t(end)])
    drawnow

end

% credible intervals
%--------------------------------------------------------------------------
t     = (1:nT) + T - 1;
for i = 1:nT
    y       = squeeze(pY(:,i,:))';
    Ey(:,i) = mean(y);
    Cy(:,i) = var(y);
end
subplot(2,1,2), set(gca,'ColorOrderIndex',1)
for i = 1:n
    spm_plot_ci(Ey(i,:),Cy(i,:),t), hold on
end
set(gca,'ColorOrderIndex',1)
plot(1:T,DEM.Y), hold off


return
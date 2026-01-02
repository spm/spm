function [Ey,Cy,Vy,pY] = spm_NESS_forecasting(DEM,nT,n,m)
% Generates paths into the future over nT time points
% FORMAT [Ey,Cy,Vy,pY] = spm_NESS_forecasting(DEM,nT,n,m)
%--------------------------------------------------------------------------
% DEM - inverted DEM structure
% nT  - number of future time points
% n   - number of enslaving states
% m   - number of enslaved states
%
% Ey  - posterior expectation
% Cy  - posterior covariances
% Vy  - posterior covariances
%
% pY(:,:,i) - i-th path
%
% This routine takes a DEM structure that has been inverted given some data
% Y and generates trajectories into the future over nT time points. The
% sample density is then recorded in terms of the first Ey and second order
% Cy moments. The ensuing trajectories are plotted; along with the
% preceding timeseries data (and posteriors over latent states).
%__________________________________________________________________________


% current time (weeks)
%--------------------------------------------------------------------------
if isfield(DEM,'T'), T = DEM.T;  else, T = size(DEM.Y,2);  end

% number of states and samples
%--------------------------------------------------------------------------
V     = 0;                         % calibration
N     = 32;                        % number of paths
nx    = spm_length(DEM.M(1).x);    % number of states
np    = spm_length(DEM.M(1).pE);   % number of parameters
nY    = size(DEM.Y,1);             % number of outcomes
r     = size(DEM.Y,2);
r     = (1:r)  + T - r;            % past time points
t     = (1:nT) + T - 1;            % future time points


% plot the past
%==========================================================================
if isfield(DEM,'qP')

    % use posteriors
    %----------------------------------------------------------------------
    pE    = DEM.M(1).pE;            % prior     expectation of parameters
    Ep    = DEM.qP.P(1);            % posterior expectation of parameters
    Cp    = DEM.qP.C;               % posterior covariances of parameters
    Ex    = DEM.qU.x{1}(:,end);     % expectation of inital states
    Cx    = DEM.qU.S{end};          % covariances of inital states

    spm_DEM_qU(DEM.qU)
    subplot(2,2,1),        set(gca,'ColorOrderIndex',1), hold on
    plot(r,DEM.Y','.'),    set(gca,'ColorOrderIndex',1)

else

    % use posteriors
    %----------------------------------------------------------------------
    pE    = DEM.M(1).pE;            % expectation of parameters
    Ep    = DEM.M(1).pE;            % expectation of parameters
    Cp    = DEM.M(1).pC;            % covariances of parameters
    Ex    = DEM.M(1).x;             % inital states
    Cx    = 0;                      % covariances of inital states

    % scale responses Y to obtain latent states X
    %----------------------------------------------------------------------
    scale = diag(DEM.M(1).pE.scale);
    X     = DEM.Y'/scale;

    subplot(2,2,1),       set(gca,'ColorOrderIndex',1)
    plot(r,DEM.Y','.'),   set(gca,'ColorOrderIndex',1), hold on
    title('Response variables','FontSize',14), axis square 

    subplot(2,2,2),       set(gca,'ColorOrderIndex',1)
    plot(r,X'),           set(gca,'ColorOrderIndex',1), hold on
    title('Latent states','FontSize',14), axis square 

end

Ep  = spm_vec(Ep);
if isstruct(Cp)
    Cp = diag(spm_vec(Cp));
end
Cx  = full(Cx);
Cp  = full(Cp);


% future projections
%==========================================================================
yLim  = [min(DEM.Y,[],'all'), max(DEM.Y,[],'all')]*2;
xLim  = [-4 1]*nT + T;

pY    = zeros(nY,nT,N);
for i = 1:N

    % generate trajectory from posterior initial conditions
    %----------------------------------------------------------------------
    for j = 1:8
        x0         = Ex + spm_sqrtm(Cx)*randn(nx,1)*V;
        P          = Ep + spm_sqrtm(Cp)*randn(np,1)*V;
        P          = spm_unvec(P,pE);
        DEM.M(1).x = x0;
        FEM        = spm_DEM_generate(DEM.M,nT,P);
        if all(all(isfinite(FEM.Y)))
            break
        end
    end
    pY(:,:,i)  = full(FEM.Y);

    % projections
    %----------------------------------------------------------------------
    subplot(2,2,1),           set(gca,'ColorOrderIndex',1), hold on
    plot(t,FEM.Y',':'),       set(gca,'ColorOrderIndex',1)
    set(gca,'XLim',xLim)
    set(gca,'YLim',yLim)

    % hidden states
    %----------------------------------------------------------------------
    subplot(2,2,2),           set(gca,'ColorOrderIndex',1), hold on
    plot(t,FEM.pU.x{1}',':'), set(gca,'ColorOrderIndex',1)
    set(gca,'XLim',xLim)
    drawnow,  fprintf('% i of %i\n',i,N)

end

% credible intervals
%--------------------------------------------------------------------------
for i = 1:nT
    y       = squeeze(pY(:,i,:))';
    Ey(:,i) = mean(y);
    Cy{i}   = cov(y);
    Vy(:,i) = var(y);
end

% enslaving outcomes
%==========================================================================
if nargin < 3, return; end

% scaled legacy (training) data
%--------------------------------------------------------------------------
subplot(4,1,3), set(gca,'ColorOrderIndex',1), hold on
in    = 1:n;
for i = in
    spm_plot_ci(Ey(i,:),Vy(i,:),t), hold on
end
set(gca,'ColorOrderIndex',1)
plot(r,DEM.Y(in,:),'LineWidth',2), hold off
title('enslaving forecast','FontSize',14)
set(gca,'XLim',xLim)

% enslaved outcomes
%==========================================================================
im    = (1:m) + n;
for i = im
    subplot(4,1,4),  hold on,       set(gca,'ColorOrderIndex',i)
    spm_plot_ci(Ey(i,:),Vy(i,:),t), set(gca,'ColorOrderIndex',i)
    hold on, plot(r,DEM.Y(i,:),'LineWidth',2)
    title('enslaved forecast','FontSize',14)
    set(gca,'XLim',xLim)
end

return


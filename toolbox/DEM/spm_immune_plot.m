function y = spm_immune_plot(P,c,U,Y)
% Plotting for immune model
% FORMAT y = spm_immune_plot(P,c,M,U)
% P - Priors
% c - covariance
% U - inputs (timing of measurements)
% Y - data
%__________________________________________________________________________
% Copyright (C) 2020 Wellcome Centre for Human Neuroimaging
 
% Thomas Parr
% Copyright (C) 2008-2022 Wellcome Centre for Human Neuroimaging


% Simulate model
%--------------------------------------------------------------------------
[y,X] = spm_immune_gen(P);

% Timing
%--------------------------------------------------------------------------
T = size(X{1},1);                               % Maximum time in hours
t = (1:T)/24;                                   % Converted to days

% Create covariance matrix
%--------------------------------------------------------------------------
c = diag(spm_vec(c));

% changes in outcomes with respect to parameters
%--------------------------------------------------------------------------
try, M.T = max(U)*24; catch, M.T = 1920;             end
try, U;               catch, U   = [1:M.T]/24;       end

% Plotting
%--------------------------------------------------------------------------
subplot(3,2,1)
plot(t,[X{1}(:,2:3),X{5}(:,2:3)])
title('Antibody response')
legend('Non-neutralising (IgM)','neutralising (IgM)','Non-neutralising (IgG)','neutralising (IgG)')
ylabel('Proportion')
xlabel('Time (days)')
axis square

subplot(3,2,2)
plot(t,[X{2}(:,2:5) X{3}(:,2:3)])
title('Cell responses')
legend('Plasma (IgM)','Immature Memory','Mature Memory','Plasma (IgG)','CD4+ T-cells','CD8+ T-cells')
ylabel('Proportion')
xlabel('Time (days)')
axis square

subplot(3,2,3)
plot(t,X{4}(:,2:3))
title('Viral load')
legend('Extracellular','Intracellular')
ylabel('Proportion')
xlabel('Time (days)')
axis square


dYdP = spm_diff(@(P,M,U)spm_immune_gen(P,M,U),P,M,U,1);

% conditional covariances
%--------------------------------------------------------------------------
Ny    = size(y,2);
for i = 1:Ny
    if iscell(dYdP)
        for j = 1:size(dYdP,2)
            D{j} = dYdP{j}(:,i);
        end
        dydp{i}  = spm_cat(D);
    else
        dydp{i}  = dYdP;
    end
    C{i}     = dydp{i}*c*dydp{i}';
end

subplot(3,2,4)
spm_plot_ci(y(:,2)',C{2},t), hold on
spm_plot_ci(y(:,1)',C{1},t)

set(gca,'ColorOrderIndex',1)
if nargin >3
    plot(U,Y(U*24,1:2),'.')
    legend('IgM','IgM','IgG','IgG','IgM','IgM','IgG','IgG')
else
    legend('IgM','IgM','IgG','IgG')
end
title('Antibodies')
ylabel('Antibody level')
xlabel('Time (days)')
yl = ylim;
ylim([0,yl(2)]);

subplot(3,2,5)
spm_plot_ci(50-4*log(y(:,3)'),16*diag(C{3})./(y(:,3).^2),t) 
% The - log relates this to the threshold cycle, and the division of the covariance uses the chain rule to rexpress this as a log under first order approximations
% The 50 is chosen as a heuristic baseline under the assumption that more than 40 cycles implies negligible concentration. 
hold on
set(gca,'ColorOrderIndex',1)
if nargin >3
    plot(U,Y(U*24,3),'.')
end
set(gca,'Ydir','reverse')
title('Viral load')
ylabel('Threshold cycle (Ct)')
xlabel('Time (days)')
yl = ylim;
ylim([0,42]);

subplot(3,2,6)
spm_plot_ci(y(:,5)',C{5},t), hold on
spm_plot_ci(y(:,4)',C{4},t)

hold on
set(gca,'ColorOrderIndex',1)
if nargin >3
    plot(U,Y(U*24,3),'.')
    legend('IFN-\gamma CD4+','IFN-\gamma CD4+','IFN-\gamma CD8+','IFN-\gamma CD8+','IFN-\gamma CD4+','IFN-\gamma CD4+','IFN-\gamma CD8+','IFN-\gamma CD8+')
else
    legend('IFN-\gamma CD4+','IFN-\gamma CD4+','IFN-\gamma CD8+','IFN-\gamma CD8+')
end
title('Cell mediated immunity')
ylabel('IFN-\gamma + /100 cells','Interpreter','Tex')
xlabel('Time (days)')
yl = ylim;
ylim([0,yl(2)]);

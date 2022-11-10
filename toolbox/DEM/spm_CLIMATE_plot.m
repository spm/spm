function spm_CLIMATE_plot(Y,X,U,T,A)
% Graphics for climate simulations
% FORMAT spm_CLIMATE_plot(Y,X,U,T,A)
% Y      - expected timeseries
% X      - latent states
% U      - indices of outcome
% T      - dates (date numbers)
% A      - data structure
%
% This auxiliary routine plots the trajectory of outcome variables and
% underlying latent or hidden states. The top panel corresponds to the
% posterior predicted expectation of the requested outcome while the
% subsequent panels show the (posterior expectations of) latent states over
% time, in groups of three. If a data structure is supplied, the
% appropriate empirical data will be superimposed over the predicted
% outcomes.
%__________________________________________________________________________
% Copyright (C) 2020 Wellcome Centre for Human Neuroimaging

% Karl Friston
% Copyright (C) 2008-2022 Wellcome Centre for Human Neuroimaging

% Plot outcomes
%==========================================================================
global CHOLD, if isempty(CHOLD); CHOLD = 1; end

% plot hidden states and outcomes
%==========================================================================

% factors and names
%--------------------------------------------------------------------------
[~,~,str] = spm_CLIMATE_priors;

% graphics
%--------------------------------------------------------------------------
subplot(3,1,1), set(gca,'ColorOrderIndex',1); hold on
t = (1:size(Y,1))/12;
if exist('T','var')
    plot(T,Y)
    datetick('x','yyyy')
    if exist('A','var')
        set(gca,'ColorOrderIndex',1);
        for j = 1:numel(U)
            plot(A(U(j)).date,A(U(j)).Y,'o')
        end
    end
else
    plot(t,Y)
end

xlabel('time (years)'),ylabel('outcome')
title('Outcomes','FontSize',16)
legend(str.outcome(U)), legend('boxoff'), box off

% Latent states
%--------------------------------------------------------------------------
k = 1;
l = {};
subplot(3,2,2 + k), set(gca,'ColorOrderIndex',1);
for i = 1:size(X,2)
    
    plot(t,X(:,i)), hold on
    ylabel('state')
    title('Latent states','FontSize',12), set(gca,'XLim',[0, t(end)])
    l = [l, str.states(i)];
    
    if ~rem(i,3)
        legend(l), legend('boxoff'), box off
        k = k + 1;
        l = {};
        if i < size(X,2)
            subplot(3,2,2 + k), set(gca,'ColorOrderIndex',1);
        end
    end
end

drawnow

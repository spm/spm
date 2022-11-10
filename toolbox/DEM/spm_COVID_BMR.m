function spm_COVID_BMR(DCM)
% Bayesian model reduction for COVID models
% FORMAT spm_COVID_BMR(DCM)
% DCM - dynamic causal model for covid outbreak
%
% This subroutine applies Bayesian model reduction to a DCM for the corona
% virus outbreak, asking whether any parameters can be treated as fixed
% parameters by reducing its prior variance to 0. Finally, the optimum
% priors are identified by applying discrete levels of shrinkage priors to
% each parameter.
%__________________________________________________________________________
% Copyright (C) 2020 Wellcome Centre for Human Neuroimaging

% Karl Friston
% Copyright (C) 2008-2022 Wellcome Centre for Human Neuroimaging

% setup
%==========================================================================

% can any parameters be removed?
%--------------------------------------------------------------------------
pE   = DCM.M.pE;
pC   = DCM.M.pC;
qE   = DCM.Ep;
qC   = DCM.Cp;

% range of shrinkage priors (V) to consider
%--------------------------------------------------------------------------
V     = 4:-1:-16; 
name  = fieldnames(pE);
nP    = numel(name);
nV    = numel(V);
L     = zeros(1,nP);                % likelihood
J     = zeros(nV,nP);
pV    = zeros(1,nP);
qV    = zeros(1,nP);

for i = 3:numel(name)
    
    % if a free parameter
    %----------------------------------------------------------------------
    if pC.(name{i})(1)
        
        % could this be a fixed parameter? (i.e.  vanishing prior variance)
        %------------------------------------------------------------------
        rE   = pE;
        rC   = pC;
        rC.(name{i}) = rC.(name{i})/exp(8);
        
        F    = spm_log_evidence(qE,qC,pE,pC,rE,rC);
        L(i) = F;
        
        % and what would have been the best prior variance?
        %------------------------------------------------------------------
        for j = 1:numel(V)   
            rE     = pE;
            rC     = pC;
            rC.(name{i})(:) = exp(V(j));
            F      = spm_log_evidence(qE,qC,pE,pC,rE,rC);
            J(j,i) = F;
        end
        
        % current log variance
        %------------------------------------------------------------------
        pV(i) = log(pC.(name{i})(1));
        
        % optimal log variance (for these data)
        %------------------------------------------------------------------
        [d,j] = max(J(:,i));
        qV(i) = V(j);
        
    end
end

% plot and display results of Bayesian model comparison
%--------------------------------------------------------------------------
spm_figure('GetWin','BMR-COVID'); clf;
subplot(2,1,1)
bar(L), set(gca,'XLim',[1 nP]), set(gca,'YLim',[-32 32]), hold on
plot([0 nP],[-3,-3],'-.')
title('A fixed effect?','Fontsize',16)
xlabel('parameter'),ylabel('change in log evidence')
axis square, box off

disp('consider fixing the following parameters'), disp(' ')
disp(name(L > -1))

subplot(2,1,2)
imagesc(1:nP,V,1 - spm_softmax(J))
title('Optimal prior variances','Fontsize',16)
xlabel('parameter'),ylabel('log prior variance')
hold on, plot(1:nP,pV,'.g','Markersize',32)
hold on, plot(1:nP,qV,'.r','Markersize',32)
axis square, box off
legend({'current','optimal'})







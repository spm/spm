function spm_COVID_BMR(DCM)
% state dependent probability transition matrices
% FORMAT spm_COVID_BMR(DCM)
% DCM - dynamic causal model for covid outbreak
%
% This subroutine applies Bayesian model reduction to a DCM for the corona
% virus outbreak, asking whether any parameters can be removed - or can be
% treated as fixed parameters by reducing their prior covariance to 0.
% Finally, the optimum priors are identified by applying discrete levels of
% shrinkage priors to each parameter.
%__________________________________________________________________________
% Copyright (C) 2020 Wellcome Centre for Human Neuroimaging

% Karl Friston
% $Id: spm_COVID_BMR.m 8175 2021-11-01 15:46:35Z guillaume $

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
V     = -(6:2:16); 
name  = fieldnames(pE);
nP    = numel(name);
nV    = numel(V);
L     = zeros(1,nP);                % likelihood
K     = zeros(1,nP);
J     = zeros(nV,nP);
pV    = zeros(1,nP);
for i = 3:numel(name)
    
    % if a free parameter
    %----------------------------------------------------------------------
    if pC.(name{i})(1)
        
        % would this be a fixed parameter?
        %------------------------------------------------------------------
        rE   = pE;
        rC   = pC;
        rC.(name{i}) = rC.(name{i})/exp(8);
        
        F    = spm_log_evidence(qE,qC,pE,pC,rE,rC);
        L(i) = F;
        
        % could it be removed?
        %------------------------------------------------------------------
        rE.(name{i}) = rE.(name{i}) - 8;
        F    = spm_log_evidence(qE,qC,pE,pC,rE,rC);
        K(i) = F;
        
        % and what would have been the best prior variance?
        %------------------------------------------------------------------
        pV(i) = log(rC.(name{i})(1));
        for j = 1:numel(V)   
            rE     = pE;
            rC     = pC;
            rC.(name{i})(:) = exp(V(j));
            F      = spm_log_evidence(qE,qC,pE,pC,rE,rC);
            J(j,i) = F;
        end       
    end
end

% plot and display results of Bayesian model comparison
%--------------------------------------------------------------------------
spm_figure('GetWin','BMR-COVID'); clf;
subplot(2,2,1)
bar(L), set(gca,'YLim',[-32 32])
title('A fixed effect?','Fontsize',16)
xlabel('parameter'),ylabel('change in log evvidence')
disp(name(L > 0))

subplot(2,2,2)
bar(K), set(gca,'YLim',[-32 32])
title('A redundant effect?','Fontsize',16)
xlabel('parameter'),ylabel('change in log evvidence')
disp(name(K > 0))

subplot(2,1,2)
imagesc(1:nP,V,spm_softmax(J))
title('Optimal prior variances','Fontsize',16)
xlabel('parameter'),ylabel('log prior variance')
hold on, plot(1:nP,pV,'.','Markersize',32)










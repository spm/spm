function [BMA] = spm_dcm_peb_bmc(PEB,models)
% hierarchical (PEB) dynamic model comparison and averaging
% FORMAT [BMA] = spm_dcm_peb_bmc(PEB,models)
%
% PEB -  between subject (second level) effects (from spm_dcm_peb)
% ------------------------------------------------------------
%     PEB.Snames - string array of Ns first level model names
%     PEB.Pnames - string array of Np parameters of interest
%
%     PEB.M.X  -   second level (between subject) design matrix
%     PEB.M.W  -   second level (within  subject) design matrix
%     PEB.M.Q  -   precision components of second level random effects
%     PEB.M.pE -   prior expectation of second level parameters
%     PEB.M.pC -   prior covariance  of second level parameters
%     PEB.Ep   -   posterior expectation of second level parameters
%     PEB.Cp   -   posterior covariance  of second level parameters
%
% models - parameter field in HCM.Ep to compare [default: 'B']
%          or logical (Nm x Np) matrix of Nm (parameteric) model space
%          or an array of DCMs specifying Nm (parameteric) model space
%
% BMA    - DCM structure of Bayesian model average
% -------------------------------------------------------------
%     BMA.Snames - string array of first level model names
%     BMA.Pnames - string array of parameters of interest
%     BMA.Pind
%
%     BMA.SUB  - first level (within subject)
%         SUB(i).Ep - posterior expectation under BMA of empirical priors
%         SUB(i).Cp - posterior covariances under BMA of empirical priors
%         SUB(i).F  - (reduced) free energy under BMA of empirical priors
%
%     BMA.Ep   - BMA expectation of second level parameters
%     BMA.Cp   - BMA covariances of second level parameters
%     BMA.F    - free energy over model space
%     BMA.P    - posterior probability over models
%     BMA.Px   - posterior probability over parameters (differences)
%     BMA.Pw   - posterior probability over parameters (common)
%     BMA.M    - second level model
%     BMA.K    - model space
%
%--------------------------------------------------------------------------
% This routine performs Bayesian model comparison & averaging at the second
% level of a hierarchical (PEB) model. The model space is defined either
% in terms of fields (e.g. 'A' or 'B') or as a logical matrix, with one row
% per model and a column per parameter (in PEB.Pnames). This induces
% a joint model space over parameters and group effects at the second level
% (encoded by the design matrix, X). Using Bayesian model reduction, this
% joint model space is scored for all combinations of models for the first
% group effect (cconstant terms modelling effects that are common to all
% subjects) and the second group effect (generally modelling between
% subject differences). The particular form of Bayesian model comparison
% and averaging here evaluates the Bayesian model average of the second
% level parameters and uses these as BMA estimates of empirical priors to
% compute subject specific posteriors (returned in BMA.SUB). Using
% partitions of model space one can then computes the posterior probability
% of various combinations of group effects over different parameters. Of
% particular interest are (i) the posterior probabilities over the
% the first two group effects in the design matrix and the posterior
% probability of models with and without each parameter, for the common
% (first) and subject-specific (second) group affects (returned in BMA.P,
% BMA.Pw and BMA.Px respectively. The Bayesian model averages of the second
% level parameters and can be found in BMA.Ep and BMA.Cp.
%
% NB for EEG models the absence of a connection means it is equal to its
% prior mesn, not that is is zero.
%
% see also: spm_dcm_peb.m and spm_dcm_bmr
%__________________________________________________________________________
% Copyright (C) 2005 Wellcome Trust Centre for Neuroimaging

% Karl Friston
% $Id: spm_dcm_peb_bmc.m 6317 2015-01-25 15:15:40Z karl $

% Compute reduced log-evidence
%==========================================================================


% number of parameters and effects
%--------------------------------------------------------------------------
[Np Nx]   = size(PEB.Ep);
if nargin < 2, models = 'B'; end
if ischar(models)
    
    % compare all combinations of field in 'models'
    %----------------------------------------------------------------------
    Pnames = char(PEB.Pnames);
    k      = any(ismember(Pnames,models),2);
    K      = ones(2^sum(k),Np);
    K(:,k) = spm_perm_mtx(sum(k));
    
elseif iscell(models)
    
    % (RFX) BMA – define the model space in terms of a matrix
    %----------------------------------------------------------------------
    Nm    = length(models);
    Np    = length(PEB.Pind);
    K     = ones(Nm,Np);
    for i = 1:Nm
        
        pC = models{i}.M.pC;
        if isstruct(pC)
            k  = find(spm_vec(pC));
        else
            k  = find(diag(pC));
        end
        j      = find(~ismember(PEB.Pind,k));
        K(i,j) = 0;
    end
    
else
    % model space in defined in terms of a matrix
    %----------------------------------------------------------------------
    K     = models;
    k     = find(any(~K));
end
[Nm Np]   = size(K);


% check number of models
%--------------------------------------------------------------------------
i = find(any(~K),1);
if (Nm > 32 && Nx > 1) || (Nm > 1024 && Nx < 2)
    warndlg('please reduce the size of your model space')
    return
end
if isempty(i)
    warndlg('your model space is empty')
    return
end


%-score models with log-evidences
%==========================================================================
fprintf('BMC:     ')

% Get priors and posteriors - of first and second order parameters
%--------------------------------------------------------------------------
qE    = spm_vec(PEB.Ep,PEB.Eh);
qC    = PEB.Cph;
pE    = spm_vec(PEB.M.pE,PEB.M.hE);
pC    = blkdiag(PEB.M.pC,PEB.M.hC);
q     = 1:Nx*Np;
for i = 1:Nm
    
    if Nx > 1
        
        % model comparison over common (constant) and group effects
        %------------------------------------------------------------------
        for j = 1:Nm
            
            % reduced prior
            %--------------------------------------------------------------
            k   = [K(i,:) K(j,:) ones(1,(Nx - 2)*Np)];
            R   = diag(k);
            rE  = spm_vec(R*PEB.M.pE,  PEB.M.hE);
            rC  = blkdiag(R*PEB.M.pC*R,PEB.M.hC);
            
            % Bayesian model reduction (of second level)
            %--------------------------------------------------------------
            [F, sE, sC] = spm_log_evidence_reduce(qE,qC,pE,pC,rE,rC);
            BMR{i,j}.Ep = sE(q);
            BMR{i,j}.Cp = sC(q,q);
            BMR{i,j}.F  = F;
            G(i,j)      = F;
            
            % report progress
            %--------------------------------------------------------------
            fprintf('\b\b\b\b%-3.0f%%',100*((i - 1)*Nm + j)/(Nm*Nm))
            
        end
        
    else
        
        % otherwise, reduced prior over group mean
        %------------------------------------------------------------------
        k   = K(i,:);
        R   = diag(k);
        
        rE  = spm_vec(R*PEB.M.pE,  PEB.M.hE);
        rC  = blkdiag(R*PEB.M.pC*R,PEB.M.hC);
        
        % Bayesian model reduction (of second level)
        %------------------------------------------------------------------
        [F, sE, sC] = spm_log_evidence_reduce(qE,qC,pE,pC,rE,rC);
        BMR{i,j}.Ep = sE(q);
        BMR{i,j}.Cp = sC(q,q);
        BMR{i,j}.F  = F;
        G(i,j)      = F;
    end
end

% family wise inference over models and parameters
%==========================================================================
P     = G;
P(:)  = exp(P(:) - max(P(:)));
P(:)  = P/sum(P(:));

% family wise inference over parameters (present an absent)
%--------------------------------------------------------------------------
k     = find(any(~K));
Nb    = length(k);
Kname = PEB.Pnames(k);
for i = 1:Nb
    Pw(1,i) = mean(sum(P( ~K(:,k(i)),:),2));
    Pw(2,i) = mean(sum(P(~~K(:,k(i)),:),2));
    
    if Nx > 1
        Px(1,i) = mean(sum(P(:, ~K(:,k(i))),1));
        Px(2,i) = mean(sum(P(:,~~K(:,k(i))),1));
    else
        Px(1,i) = 1;
        Px(2,i) = 0;
    end
    
end
Pw    = Pw(2,:)./sum(Pw,1);
Px    = Px(2,:)./sum(Px,1);

% family wise inference over mmodels (commonalities and differences)
%--------------------------------------------------------------------------
P1    = sum(P,2);
P2    = sum(P,1);

%-hierarchical inversion using optimised second level priors
%==========================================================================

% Bayesian model averaging (with an Occam's window of eight)
%--------------------------------------------------------------------------
i     = G(:) > max(G(:) - 8);
BMA   = spm_dcm_bma(BMR(i)');
for i = 1:length(PEB.SUB)
    
    % get empirical prior expectations and reduced 1st level posterior
    %----------------------------------------------------------------------
    rE        = kron(PEB.M.X(i,:),PEB.M.W)*BMA.Ep;
    rC        = PEB.Ce;
    
    pE        = PEB.SUB(i).pE;
    pC        = PEB.SUB(i).pC;
    qE        = PEB.SUB(i).Ep;
    qC        = PEB.SUB(i).Cp;
    
    [F sE sC] = spm_log_evidence_reduce(qE,qC,pE,pC,rE,rC);
    
    % and save
    %----------------------------------------------------------------------
    SUB(i).Ep = sE;
    SUB(i).Cp = sC;
    SUB(i).F  = F;
    
end


% assemble BMA output structure
%--------------------------------------------------------------------------
BMA.Sname = PEB.Snames;
BMA.Pname = PEB.Pnames;
BMA.Pind  = PEB.Pind;
BMA.Kname = Kname;

BMA.SUB   = SUB;
BMA.F     = G;
BMA.P     = P;
BMA.Px    = Px;
BMA.Pw    = Pw;
BMA.M     = PEB.M;
BMA.K     = K;


% Show results
%==========================================================================
spm_figure('Getwin','BMC'); clf

subplot(3,2,1), imagesc(G)
title('Free energy','FontSize',16)
xlabel('Model (differences)','FontSize',12)
ylabel('Model (commonalities)','FontSize',12)
axis square

subplot(3,2,3)
[m i] = max(P1); bar(P1),
text(i - 1/4,m/2,sprintf('%-2.0f%%',m*100),'Color','w','FontSize',8)
title('Commonalities','FontSize',16)
xlabel('Model','FontSize',12)
ylabel('Probability','FontSize',12)
axis([0 (Nm + 1) 0 1]), axis square

subplot(3,2,5), bar(diag(Pw),length(Pw));
title('Commonalities','FontSize',16)
xlabel('Parameter','FontSize',12)
ylabel('Model probability','FontSize',12)
axis([0 (Nb + 1) 0 1]), axis square
legend(Kname)

if Nx < 2, return, end

subplot(3,2,2), imagesc(P)
title('Posterior probabilities','FontSize',16)
xlabel('Model (differences)','FontSize',12)
ylabel('Model (commonalities)','FontSize',12)
axis square

subplot(3,2,4)
[m i] = max(P2); bar(P2),
text(i - 1/4,m/2,sprintf('%-2.0f%%',m*100),'Color','w','FontSize',8)
title('Differences','FontSize',16)
xlabel('Model','FontSize',12)
ylabel('Probability','FontSize',12)
axis([0 (Nm + 1) 0 1]), axis square

subplot(3,2,6), bar(diag(Px),length(Px))
title('Differences','FontSize',16)
xlabel('Parameter','FontSize',12)
ylabel('Model probability','FontSize',12)
axis([0 (Nb + 1) 0 1]), axis square



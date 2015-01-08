function [BMA] = spm_dcm_group_BMA(HDM,models)
% Hierarchical dynamic model comparison and averaging
% FORMAT [BMA] = spm_dcm_group_BMA(HDM,models)
%
% HCM    - GCM of between subject (second level) effects
% ------------------------------------------------------------
%     HDM.Snames - string array of Ns first level model names
%     HDM.Pnames - string array of Np parameters of interest
%
%     HDM.M.X  -   second level (between subject) design matrix
%     HDM.M.Q  -   precision components of second level random effects 
%     HDM.M.pE -   prior expectation of second level parameters
%     HDM.M.pC -   prior covariance  of second level parameters
%     HDM.Ep   -   posterior expectation of second level parameters
%     HDM.Cp   -   posterior covariance  of second level parameters
%
% models - parameter field in HCM.Ep to compare [default: {'B'}]
%          or logical (Nm x Np) matrix of Nm models
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
%     BMA.Ep   - BMA expectation of ssecond level parameters
%     BMA.Cp   - BMA covariances of ssecond level parameters
%     BMA.F    - free energy over model space
%     BMA.P    - posterior probability over models
%     BMA.Px   - posterior probability over parameters (differences)
%     BMA.Pw   - posterior probability over parameters (common)
%     BMA.M    - second level model
%     BMA.K    - model space  
%
%--------------------------------------------------------------------------
% This routine performs Bayesian model comparison & averaging at the second
% level of a hierarchical dynamic model. The model space is defined either
% in terms of fields (e.g. 'A' or 'B') or as a logical matrix, with one row
% per model and one her parameter (in HDM.Pnames). This model space induces
% a joint model space over parameter and group effects at the second level
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
% particular interest are (i) the posterior probabilities with and without 
% the first two group effects in the design matrix and the posterior
% probability of models with and without each parameter, for the common
% (first) and subject-specific (second) group affects (returned in BMA.PX,
% BMA.Pw and BMA.Px respectively. The Bayesian model averages of the second
% level parameters and can be found in BMA.Ep and BMA.Cp.
%
% see also: spm_dcm_group.m
%__________________________________________________________________________
% Copyright (C) 2015 Wellcome Trust Centre for Neuroimaging
 
% Karl Friston
% $Id: spm_dcm_group_BMA.m 6299 2015-01-08 12:56:00Z guillaume $
 
% Compute reduced log-evidence
%==========================================================================


%  mmodel space of within subject effects
%--------------------------------------------------------------------------
[Np,Nx]   = size(HDM.Ep);             % number of parameters and effects
Pname     = char(HDM.Pnames);
if nargin == 2;
    if ischar(models)
        field = models;               % compare all combinations of field         
    else
        K     = models;               % specified model space
        k     = find(any(~K));
        field = Pname(k(1),1);
    end
else
    field  = 'B';
    k      = any(ismember(Pname,field),2); 
    K      = ones(2^sum(k),Np);       
    K(:,k) = spm_perm_mtx(sum(k));
end
[Nm,Np]  = size(K);


% check number of models
%--------------------------------------------------------------------------
if (Nm > 256 && Nx > 1) || (Nm > 1024 && Nx < 2)
    warndlg('please reduce the size of your model space')
    return
end


% prior expectations under the null (for extrinsic EEG connections)
%--------------------------------------------------------------------------
if field == 'A' && isfield(HDM.SUB(1).M.pE,'T')
    rD = -4*kron(spm_speye(Nx,1),k');
else
    rD = 0;
end
    
%-score models with log-evidences
%==========================================================================
fprintf('BMC:     ')

% Get priors and posteriors - of first and second order parameters
%--------------------------------------------------------------------------
qE    = spm_vec(HDM.Ep);
qC    = HDM.Cp;
pE    = spm_vec(HDM.M.pE);
pC    = HDM.M.pC;
for i = 1:Nm
    
    if Nx > 1
        
        % model comparison over common (constant) and group effects
        %------------------------------------------------------------------
        for j = 1:Nm
                        
            % reduced prior
            %--------------------------------------------------------------
            k   = [K(i,:) K(j,:) ones(1,(Nx - 2)*Np)];
            R   = diag(k);
            rE  = R*pE - rD.*(1 - k');
            rC  = R*pC*R;
            
            % Bayesian model reduction (of second level)
            %--------------------------------------------------------------
            [F, sE, sC] = spm_log_evidence_reduce(qE,qC,pE,pC,rE,rC);
            BMR{i,j}.Ep = sE;
            BMR{i,j}.Cp = sC;
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
        rE  = R*pE - rD.*(1 - k');
        rC  = R*pC*R;
        
        % Bayesian model reduction (of second level)
        %------------------------------------------------------------------
        [F, sE, sC] = spm_log_evidence_reduce(qE,qC,pE,pC,rE,rC);
        BMR{i,1}.Ep = sE;
        BMR{i,1}.Cp = sC;
        BMR{i,1}.F  = F;
        G(i,1)      = F;
        
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
Kname = HDM.Pnames(k);
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
bma   = spm_dcm_bma(BMR(i)');

for i = 1:length(HDM.SUB)
    
    % get empirical prior expectations and reduced 1st level posterior
    %----------------------------------------------------------------------
    rE        = kron(HDM.M.X(i,:),HDM.M.W)*bma.mEp;
    rC        = HDM.Ce;
    
    pE        = HDM.SUB(i).pE;
    pC        = HDM.SUB(i).pC;
    qE        = HDM.SUB(i).Ep;
    qC        = HDM.SUB(i).Cp;
    
    [F,sE,sC] = spm_log_evidence_reduce(qE,qC,pE,pC,rE,rC);
        
    % and save
    %----------------------------------------------------------------------
    SUB(i).Ep = sE;
    SUB(i).Cp = sC;
    SUB(i).F  = F;
    
end


% assemble BMA output structure
%--------------------------------------------------------------------------
BMA.Sname = HDM.Snames;
BMA.Pname = HDM.Pnames;
BMA.Pind  = HDM.Pind;
BMA.Kname = Kname;

BMA.SUB   = SUB;
BMA.Ep    = spm_unvec(bma.mEp,HDM.Ep);
BMA.Cp    = spm_unvec(bma.sEp.^2,HDM.Ep);
BMA.F     = G;
BMA.P     = P;
BMA.Px    = Px;
BMA.Pw    = Pw;
BMA.M     = HDM.M;
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
[m,i] = max(P1); bar(P1),
text(i - 1/4,m/2,sprintf('%-2.0f%%',m*100),'Color','w','FontSize',8)
title('Commonalities','FontSize',16)
xlabel('Model','FontSize',12)
ylabel('Probability','FontSize',12)
axis([0 (Nm + 1) 0 1]), axis square

subplot(3,2,5), bar(diag(Pw),length(Pw)); 
title('Commonalities','FontSize',16)
xlabel('Paameter','FontSize',12)
ylabel('Mosterior probability','FontSize',12)
axis([0 (Nb + 1) 0 1]), axis square
legend(Kname)

if Nx < 2, return, end

subplot(3,2,2), imagesc(P)
title('Posterior probabilities','FontSize',16)
xlabel('Model (differences)','FontSize',12)
ylabel('Model (commonalities)','FontSize',12)
axis square

subplot(3,2,4)
[m,i] = max(P2); bar(P2),
text(i - 1/4,m/2,sprintf('%-2.0f%%',m*100),'Color','w','FontSize',8)
title('Differences','FontSize',16)
xlabel('Model','FontSize',12)
ylabel('Probability','FontSize',12)
axis([0 (Nm + 1) 0 1]), axis square

subplot(3,2,6), bar(diag(Px),length(Px))
title('Differences','FontSize',16)
xlabel('Parameter','FontSize',12)
ylabel('Mosterior probability','FontSize',12)
axis([0 (Nb + 1) 0 1]), axis square

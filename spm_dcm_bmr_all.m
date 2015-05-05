function [DCM,BMR,BMA] = spm_dcm_bmr_all(DCM,field)
% Bayesian model reduction of all permutations of model parameters
% FORMAT [RCM,BMR,BMA] = spm_dcm_bmr_all(DCM,field)
%
%  DCM      - DCM structures:
%  DCM.M.pE - prior expectation (with parameters in pE.A, pE.B and pE.C)
%  DCM.M.pC - prior covariance
%  DCM.Ep   - posterior expectation: Bayesian model average
%  DCM.Cp   - posterior covariances; Bayesian model average
%  DCM.Pp   - Model posterior (with and without each parameter)
%
% field     - parameter fields in DCM{i}.Ep to optimise [default: {'A','B'}]
%             'All' will invoke all fields (i.e. random effects)
%             If Ep is not a structure, all parameters will be considered
%
%
% RCM - reduced DCM array
% BMR - (Nsub) summary structure 
%        BMR.name - character/cell array of DCM filenames
%        BMR.F    - their associated free energies
%        BMR.P    - and posterior (model) probabilities
% BMA - Baysian model average (see spm_dcm_bma)
%
%--------------------------------------------------------------------------
% This routine searches over all possible reduced models of a full model
% (DCM) and uses Bayesian model reduction to model average. Reduced
% models mean all permutations of free parameters (parameters with a non-
% zero prior covariance), where models are defined in terms of their prior
% covariance. The full model should be inverted prior to post hoc
% optimization. If there are more than 16 free-parameters, this routine
% will implement a greedy search: This entails searching over all
% permutations of the 8 parameters whose removal (shrinking the prior
% variance to zero) produces the smallest reduction (greatest increase)
% in model evidence. This procedure is repeated until all 8 parameters
% are retained in the best model or there are no more parameters to
% consider.
%
% See also: spm_dcm_post_hoc – this routine is essentially a simplified
% version of spm_dcm_post_hoc
%__________________________________________________________________________
% Copyright (C) 2010-2014 Wellcome Trust Centre for Neuroimaging

% Karl Friston, Peter Zeidman
% $Id: spm_dcm_bmr_all.m 6427 2015-05-05 15:42:35Z karl $


%-Number of parameters to consider before invoking greedy search
%--------------------------------------------------------------------------
nmax = 8;

%-Check fields of parameter stucture
%--------------------------------------------------------------------------
if nargin < 2 || isempty(field)
    field = {'A','B'};
end

%-dela with filenames stucture
%--------------------------------------------------------------------------
if ischar(DCM)
    DCM = load(DCM,'DCM');
    DCM = DCM.DCM;
end

% Get prior covariances
%--------------------------------------------------------------------------
if isstruct(DCM.M.pC), DCM.M.pC = diag(spm_vec(DCM.M.pC)); end

% Get priors and posteriors
%--------------------------------------------------------------------------
qE  = DCM.Ep;
qC  = DCM.Cp;
pE  = DCM.M.pE;
pC  = DCM.M.pC;

% Remove (a priori) null space
%--------------------------------------------------------------------------
U   = spm_svd(pC);
qE  = U'*spm_vec(qE);
pE  = U'*spm_vec(pE);
qC  = U'*qC*U;
pC  = U'*pC*U;


%-Greedy search (GS) - eliminating parameters in a top down fashion
%==========================================================================

% Accumulated reduction vector (C)
%--------------------------------------------------------------------------
q   = diag(DCM.M.pC);
C   = logical(q > mean(q(q < 1024))/1024);
GS  = 1;
while GS
    
    %-Find free coupling parameters
    %----------------------------------------------------------------------
    if isstruct(DCM.Ep)
        k = spm_fieldindices(DCM.Ep,field{:});
    else
        k = 1:spm_length(DCM.Ep);
    end
    k = k(C(k));
    
    % If there are too many find those with the least evidence
    %----------------------------------------------------------------------
    nparam = length(k);
    if nparam > nmax
        
        % Model search over new prior without the i-th parameter
        %------------------------------------------------------------------
        Z     = zeros(1,nparam);
        for i = 1:nparam
            r    = C; r(k(i)) = 0;
            R    = U(r,:)'*U(r,:);
            rE   = R*pE;
            rC   = R*pC*R;
            Z(i) = spm_log_evidence(qE,qC,pE,pC,rE,rC);
        end
        
        % Find parameters with the least evidence
        %------------------------------------------------------------------
        [z,i] = sort(-Z);
        k     = k(i(1:8));
        
        % Flag a greedy search
        %------------------------------------------------------------------
        GS = 1;
        
    elseif isempty(k)
        fprintf('\nThere are no free parameters in this model.\n')
        return
    else
        GS = 0;
    end
    
    % Create model space in terms of free parameter indices
    %----------------------------------------------------------------------
    K     = spm_perm_mtx(length(k));
    
    % Model search over new prior (covariance)
    %----------------------------------------------------------------------
    G     = [];
    for i = 1:size(K,1)
        r    = C; r(k(K(i,:))) = 0;
        R    = U(r,:)'*U(r,:);
        rE   = R*pE;
        rC   = R*pC*R;
        G(i) = spm_log_evidence(qE,qC,pE,pC,rE,rC);
    end
    
    % posterior probability
    %----------------------------------------------------------------------
    p      = spm_softmax(G(:));
    
    %-Get selected model and prune redundant parameters
    %======================================================================
    [z,i]  = max(p);
    C(k(K(i,:))) = 0;
    
    % Continue greedy search if any parameters have been eliminated
    %----------------------------------------------------------------------
    nelim  = full(sum(K(i,:)));
    GS     = GS & nelim;
    
    % Show results
    % ---------------------------------------------------------------------
    spm_figure('Getwin','BMR - all'); clf
    fprintf('%i out of %i free parameters removed \n',nelim,nparam)
    
    subplot(3,2,1)
    if length(K) > 32, plot(G,'k'), else, bar(G,'c'), end
    title('log-posterior','FontSize',16)
    xlabel('model','FontSize',12)
    ylabel('log-probability','FontSize',12)
    axis square
    
    subplot(3,2,2)
    if length(K) > 32, plot(p,'k'), else, bar(p,'r'), end
    title('model posterior','FontSize',16)
    xlabel('model','FontSize',12)
    ylabel('probability','FontSize',12)
    axis square
    drawnow
    
end


%-Inference over families (one family per coupling parameter)
%==========================================================================
for i = 1:length(k)
    Pk(1,i) = mean(p(~K(:,i)));
    Pk(2,i) = mean(p( K(:,i)));
end
Pk    = Pk(1,:)./sum(Pk);
Pp    = full(C);
Pp(k) = Pk;


%-Bayesian model average
%==========================================================================
qE    = DCM.Ep;
qC    = DCM.Cp;
pE    = DCM.M.pE;
pC    = DCM.M.pC;
pE    = spm_vec(pE);
BMA   = {};
Gmax  = max(G);
for i = 1:length(K)
    if G(i) > (Gmax - 8)
        r            = C;
        r(k(K(i,:))) = 0;
        R            = diag(r);
        rE           = R*spm_vec(pE);
        rC           = R*pC*R;
        [F,Ep,Cp]    = spm_log_evidence_reduce(qE,qC,pE,pC,rE,rC);
        BMA{end + 1} = struct('Ep',Ep,'Cp',Cp,'F',F);
    end
end
BMA   = spm_dcm_bma(BMA);
Ep    = BMA.Ep;
Cp    = BMA.Cp;


% Show full and reduced conditional estimates (for Bayesian average)
%--------------------------------------------------------------------------
spm_figure('Getwin','BMR - all');

if isstruct(DCM.Ep)
    i  = spm_find_pC(pC,DCM.Ep,field);
else
    i  = 1:spm_length(DCM.Ep);
end
qE     = spm_vec(qE);
Ep     = spm_vec(Ep);


% BMR summary and plotting
%--------------------------------------------------------------------------
try
    Pnames = spm_fieldindices(DCM.Ep,k);
catch
    Np     = numel(DCM.Pnames);
    Pnames = DCM.Pnames(rem(k - 1,Np) + 1);
end
BMR.name = Pnames;
BMR.F    = G;
BMR.P    = p;
BMR.K    = K;

subplot(3,2,3), spm_plot_ci(qE(i),qC(i,i))
title('MAP (full)','FontSize',16)
axis square, a = axis;

subplot(3,2,4), spm_plot_ci(Ep(i),abs(Cp(i)))
title('MAP (reduced)','FontSize',16), axis square, axis(a)

subplot(3,2,5), imagesc(K')
xlabel('model'), ylabel('parameter'), title('model space','FontSize',16)
set(gca,'YTickLabel',BMR.name);
axis tight, axis square

subplot(3,2,6), bar(diag(Pp(i)),length(i))
xlabel('parameter'), title(' posterior','FontSize',16)
spm_axis tight, axis square
drawnow


%-Save Bayesian parameter average and family-wise model inference
%==========================================================================
if isstruct(DCM.Ep)
    if length(i) < 32
        legend(spm_fieldindices(DCM.Ep,i))
    end
    Pp    = spm_unvec(Pp,DCM.Ep);
    Ep    = spm_unvec(Ep,DCM.Ep);
end

DCM.Pp    = Pp;        % Model posterior over parameters (with and without)
DCM.Ep    = Ep;        % Bayesian model averages
DCM.Cp    = Cp;        % Bayesian model variance
DCM.F     = DCM.F + F; % reduced free energy



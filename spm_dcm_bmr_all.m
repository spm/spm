function [DCM,BMR,BMA] = spm_dcm_bmr_all(DCM,field)
% Bayesian model reduction of all permutations of model parameters
% FORMAT [RCM,BMR,BMA] = spm_dcm_bmr_all(DCM,field)
%
% DCM      - A single estimated DCM (or PEB) structure:
%
%  DCM.M.pE  - prior expectation
%  DCM.M.pC  - prior covariance
%  DCM.Ep    - posterior expectation
%  DCM.Cp    - posterior covariances
%  DCM.beta  - prior expectation of reduced parameters (default: 0)
%  DCM.gamma - prior variance    of reduced parameters (default: 0)
%              NB: beta = 'pE' uses full priors
%
% field      - parameter fields in DCM{i}.Ep to optimise [default: {'A','B'}]
%             'All' will invoke all fields (i.e. random effects)
%             If Ep is not a structure, all parameters will be considered
%
% Returns:
%
% DCM - Bayesian Model Average (BMA) over models in the final iteration of 
%       the search:
%
%       DCM.Ep    - (BMA) posterior expectation
%       DCM.Cp    - (BMA) posterior covariance
%
% BMR -  (Nsub) summary structure reporting the model space from the last
%        iteration of the search:
%
%        BMR.name - character/cell array of parameter names
%        BMR.F    - free energies (relative to full model)
%        BMR.P    - and posterior (model) probabilities
%        BMR.K    - [models x parameters] model space (1 = off, 0 = on)
%        BMR.bma  - cell array of each model's parameters and optimised 
%                   model evidences used to calculate the BMA
%
% BMA - Baysian model average (see spm_dcm_bma)
%
%--------------------------------------------------------------------------
% This routine searches over reduced (nested) models of a full model (DCM) 
% using Bayesian model reduction and performs Bayesian Model Averaging.
% 'Reduced' means some free parameters (parameters with a non-
% zero prior covariance) are switched off by fixing their prior variance 
% to zero. 
%
% If there are fewer than nmax = 8 free parameters, all permutations of 
% switching off parameters will be tested. Otherwise, this routine 
% implements the following greedy search procedure. The nmax parameters 
% are identified which, when switched off individually, produce the least 
% reduction (greatest increase) in model evidence. All permutations of 
% switching off these parameters are then evaluated and the best 
% permutation is retained. This procedure is repeated until all nmax
% parameters are retained or there are no more parameters to consider. 
% Finally, BMA is performed on the models from the last iteration.
% 
% NB: The full model should be estimated prior to running this function. 
%
% See also: spm_dcm_post_hoc - this routine is essentially a simplified
% version of spm_dcm_post_hoc
%__________________________________________________________________________
% Copyright (C) 2010-2014 Wellcome Trust Centre for Neuroimaging

% Karl Friston, Peter Zeidman
% $Id: spm_dcm_bmr_all.m 7476 2018-11-07 15:17:39Z peter $


%-Number of parameters to consider before invoking greedy search
%--------------------------------------------------------------------------
nmax  = 8;

%-specification of null prior covariance
%--------------------------------------------------------------------------
if isfield(DCM,'beta'),  beta  = DCM.beta;  else, beta  = 0; end
if isfield(DCM,'gamma'), gamma = DCM.gamma; else, gamma = 0; end

%-Check fields of parameter stucture
%--------------------------------------------------------------------------
if nargin < 2 || isempty(field)
    field = {'A','B'};
end
if ischar(field)
    field = {field};
end

%-deal with filenames structure
%--------------------------------------------------------------------------
if ischar(DCM)
    DCM = load(DCM,'DCM');
    DCM = DCM.DCM;
end

% Get prior covariances
%--------------------------------------------------------------------------
if isstruct(DCM.M.pC), DCM.M.pC = diag(spm_vec(DCM.M.pC)); end
if spm_length(DCM.M.pE) ~= size(DCM.M.pC,1)
    DCM.M.pC = diag(spm_vec(DCM.M.pC));
end

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
if sum(q < 1024)
    C   = double(q > mean(q(q < 1024))/1024);
else
    C   = double(q > 0);
end
GS  = 1;
while GS
    
    %-Find free coupling parameters
    %----------------------------------------------------------------------
    if isstruct(DCM.Ep)
        k = spm_fieldindices(DCM.Ep,field{:});
    else
        k = 1:spm_length(DCM.Ep);
    end
    k = k(find(C(k))); %#ok<FNDSB>
    
    % If there are too many find those with the least evidence
    %----------------------------------------------------------------------
    nparam = length(k);
    if nparam > nmax
        
        % Model search over new prior without the i-th parameter
        %------------------------------------------------------------------
        Z     = zeros(1,nparam);
        for i = 1:nparam
            
            % Identify parameters to retain r and to remove s
            %--------------------------------------------------------------
            r    = C; r(k(i)) = 0; s = 1 - r;

            % Create reduced prior covariance matrix
            %--------------------------------------------------------------
            R    = U'*diag(r + s*gamma)*U;
            rC   = R*pC*R;
            
            % Create reduced prior means
            %--------------------------------------------------------------
            if isnumeric(beta)
                S    = U'*diag(r)*U;
                rE   = S*pE + U'*s*beta;
            else
                rE   = pE;
            end
            
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
        
        % Identify parameters to retain (r) and to remove (s)
        %------------------------------------------------------------------
        r    = C; r(k(K(i,:))) = 0; s = 1 - r;
        
        % Create reduced prior covariance matrix
        %------------------------------------------------------------------
        R    = U'*diag(r + s*gamma)*U;
        rC   = R*pC*R;
        
        % Create reduced prior means
        %------------------------------------------------------------------
        if isnumeric(beta)
            S    = U'*diag(r)*U;
            rE   = S*pE + U'*s*beta;
        else
            rE   = pE;
        end
        
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
Pp    = C;
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
        s            = 1 - r;
        R            = diag(r + s*gamma);
        rC           = R*pC*R;
        S            = diag(r);
        if isnumeric(beta)
            rE       = S*spm_vec(pE) + s*beta;
        else
            rE       = pE;
        end
      
        [F,Ep,Cp]    = spm_log_evidence_reduce(qE,qC,pE,pC,rE,rC);
        BMA{end + 1} = struct('Ep',Ep,'Cp',Cp,'F',F);
    end
end

BMR.bma = BMA;
BMA     = spm_dcm_bma(BMA);
Ep      = BMA.Ep;
Cp      = BMA.Cp;

if isstruct(Cp) || (spm_length(Cp) == spm_length(Ep))
    Cp = diag(spm_vec(Cp));
end

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

j = i(ismember(i,1:length(spm_vec(Ep))));

% BMR summary and plotting
%--------------------------------------------------------------------------
try
    Pnames     = spm_fieldindices(DCM.Ep,k);
catch
    try
        Np     = numel(DCM.Pnames);
        Pnames = DCM.Pnames(rem(k - 1,Np) + 1);
    catch
        Pnames = 'all parameters';
    end
end
BMR.name = Pnames;
BMR.F    = G;
BMR.P    = p;
BMR.K    = K;
BMR.k    = k;

subplot(3,2,3), spm_plot_ci(qE(i),qC(i,i))
title('MAP (full)','FontSize',16)
axis square, a = axis;

subplot(3,2,4), spm_plot_ci(Ep(j),abs(Cp(j,j)))
title('MAP (reduced)','FontSize',16), axis square, axis(a)

subplot(3,2,5), imagesc(1 - K')
xlabel('model'), ylabel('parameter'), title('model space','FontSize',16)
set(gca,'YTickLabel',BMR.name);
axis tight, axis square

subplot(3,2,6)
Np = length(i);
if Np > 1
    bar(diag(Pp(i)),Np)
else
    bar(Pp)
end
xlabel('parameter'), title(' posterior','FontSize',16)
axis square, drawnow, axis([0 (Np + 1) 0 1])


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



function DCM = spm_dcm_diagnose(DCM,varargin)
% Post hoc diagnosis of DCMs (under the Laplace approximation)
% FORMAT spm_dcm_diagnose(DCM,'field','field',...)
%
% DCM        - DCM stricture (inverted)
% field      - field name(s) of parameters to consider
%
%--------------------------------------------------------------------------
% This routine searches over all possible reduced models of a full model
% (DCM) and uses post hoc model selection to select the best. Reduced
% models mean all permutations of free parameter sets (parameters with non-
% zero prior covariance), where models are defined in terms of their prior
% covariance. The full model should be inverted prior to post hoc
% optimization. If there are more than 8 free-parameter sets, this routine
% will implement a greedy search: This entails searching over all
% permutations of the 8 sets whose removal (shrinking the prior
% variance to zero) produces the smallest reduction (greatest increase)
% in model evidence. This procedure is repeated until all sets
% are retained in the best model or there are no more parameters to
% consider.
%
% A parameter set is specified implicitly by the structure (DCM.Ep). Each
% set corresponds to a column of (the cell arrays or matrix) each field of
% DCM.Ep.
%
% if only one field is specified the log-evidence is computed as a function
% of the scaled prior variance. Redundant parameters have a log-evidence 
% that keeps increasing as the prior variance shrinks.
%
% The outputs of this routine are graphics reporting the model reduction
% (optimisation). Red means weak evidence; blue strong evidence (> 3) and 
% cyan very strong evidence (> 5)
%__________________________________________________________________________
% Copyright (C) 2008-2011 Wellcome Trust Centre for Neuroimaging
 
% Karl Friston
% $Id: spm_dcm_diagnose.m 7679 2019-10-24 15:54:07Z spm $
 
 
 
 
%-Get priors and posteriors
%==========================================================================
qE      = DCM.Ep;
qC      = DCM.Cp;
pE      = DCM.M.pE;
pC      = DCM.M.pC;
 
if isstruct(pC);
    pV  = pC;
    pC  = diag(spm_vec(pC));
end
 
 
 
% find the partition of parameters
%--------------------------------------------------------------------------
str   = {};
ind   = {};
if isstruct(qE)
    
    % get fields (parameters) to consider
    %----------------------------------------------------------------------
    if ~isempty(varargin)
        field = varargin;
    else
        field = fieldnames(qE);
    end
    for i = 1:length(field)
        q = getfield(qE,field{i});
        k = spm_fieldindices(qE,field{i});
        for j = 1:size(q,2)
            u            = 1:length(spm_vec(q(:,j)));
            str{end + 1} = [field{i} '(' num2str(j) ')'];
            ind{end + 1} = k(u); k(u) = [];
        end
    end
    
else
    field = {};
    k     = 1:length(spm_vec(qE));
    for j = 1:size(qE,2)
        u            = 1:length(spm_vec(qE(:,j)));
        str{end + 1} = ['qE(' num2str(j) ')'];
        ind{end + 1} = k(u); k(u) = [];
    end
end
 
% find free parameters
%--------------------------------------------------------------------------
ip    = {};
sp    = {};
for i = 1:length(str)
    j = find(diag(pC(ind{i},ind{i})));
    if ~isempty(j)
        ip{end + 1} = ind{i}(j);
        sp{end + 1} = str{i};
    end
end
 
 
% model search over shrinking priors
% -------------------------------------------------------------------------
n     = length(pC);
if length(field) == 1
    
    r     = (0:64)/32;
    for j = 1:length(r)
        for i = 1:length(ip)
            R      = speye(n,n) - sparse(ip{i},ip{i},1,n,n)*(1 - r(j));
            Z(i,j) = spm_log_evidence(qE,qC,pE,pC,pE,R*pC*R);
            
        end
    end
    
    subplot(2,1,1)
    plot(r,Z)
    title(['log-evidence for ' field{1}],'FontSize',16);
    xlabel('proportion of full prior')
    ylabel('negative free-energy')
    return
    
end
 
 
 
% model search over new prior without the i-th parameter set
% -------------------------------------------------------------------------
for i = 1:length(ip)
    R    = speye(n,n) - sparse(ip{i},ip{i},1,n,n);
    Z(i) = spm_log_evidence(qE,qC,pE,pC,pE,R*pC*R);
end
 
% find parameters with the least evidence
%--------------------------------------------------------------------------
[Z i] = sort(-Z);
ip    = ip(i);
sp    = sp(i);
 
% If there are too many subsets search use those with the least evidence
%--------------------------------------------------------------------------
if length(ip) > 8
    
    % flag a greedy search
    %----------------------------------------------------------------------
    ip     = ip(1:8);
    repeat = 1;
    
else
    repeat = 0;
end
 
% Create model space in terms of free parameter indices
%--------------------------------------------------------------------------
K     = spm_perm_mtx(length(ip));
 
% and get log-evidences
%==========================================================================
for i = 1:size(K,1)
    k    = spm_vec(ip(find(K(i,:))));
    R    = speye(n,n) - sparse(k,k,1,n,n);
    S(i) = spm_log_evidence(qE,qC,pE,pC,pE,R*pC*R);
end
 
 
% Find a good model with the most parameters
% =========================================================================
[i k]  = max(S);
k      = find(K(k,:));
 
nelim  = length(k);
repeat = repeat & nelim;
 
% Continue greedy search if any parameters have been eliminated
%--------------------------------------------------------------------------
if repeat
    
    fprintf('removing %i parameter set: ',nelim)
    fprintf('%s, ',sp{k});fprintf('\n')
    
    % Optimise selected model (prior covariance)
    %----------------------------------------------------------------------
    k         = spm_vec(ip(k));
    R         = speye(n,n) - sparse(k,k,1,n,n);
    rC        = R*pC*R;
    
    % Get posterior of selected model (rC)
    % ---------------------------------------------------------------------
    [F Ep Cp] = spm_log_evidence(qE,qC,pE,pC,pE,rC);
    
    DCM.M.pC  = rC;
    DCM.Ep    = Ep;
    DCM.Cp    = Cp;
    DCM.F     = F;
    DCM       = spm_dcm_diagnose(DCM,field{:});
    return
    
else
    
    % model posterior
    % ---------------------------------------------------------------------
    S     = S - S(end);
    p     = exp(S - max(S));
    p     = p/sum(p);
    
    fprintf('\n%i Irrelevant parameter sets: ',length(k))
    fprintf('%s, ',sp{k});fprintf('\n')
    
    % Organise field names and show results
    % ---------------------------------------------------------------------
    spm_figure('Getwin','Graphics');
    
    
    for i = 1:4
        j = (1:8) + (i - 1)*8;
        try
            tsr{i,1} = sprintf('%s,',sp{j});
        catch
            tsr{i,1} = sprintf('%s,',sp{j(1):end});
        end
    end
       
    subplot(2,1,1)
    bar(Z,'c'),hold on
    bar(Z(Z < 5),'b'),hold on
    bar(Z(Z < 3),'r'),hold off
    title('relative log evidence (with vs. without)','FontSize',16)
    xlabel(tsr,'FontSize',16)
    ylabel('negative free-energy','FontSize',12)
    
    subplot(2,2,3)
    if size(K,1) > 32, plot(S,'k'), else, bar(S,'c'), end
    title('log-posterior','FontSize',16)
    xlabel('model','FontSize',12)
    ylabel('log-probability','FontSize',12)
    axis square
    
    subplot(2,2,4)
    if size(K,1) > 32, plot(p,'k'), else, bar(p,'r'), end
    title('model posterior','FontSize',16)
    xlabel('model','FontSize',12)
    ylabel('probability','FontSize',12)
    axis square
    
end



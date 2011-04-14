function spm_dcm_post_hoc(P)
% Post hoc optimisation of DCMs (under Laplace approximation)
% FORMAT spm_dcm_post_hoc(P)
%
% P         - character/cell array of DCM filenames
%           - or cell array of DCM structures
%
%--------------------------------------------------------------------------
% This routine searches over all possible reduced models of a full model 
% (DCM) and uses post hoc model selection to select the best. Reduced 
% models mean all permutations of free parameters (parameters with a non-
% zero prior covariance), where models are defined in terms of their prior 
% covariance. The full model should be inverted prior to post hoc 
% optimization. If there are more than 16 free-parameters, this routine 
% will implement a greedy search: This entails searching over all 
% permutations of the 8 parameters whose removal (shrinking the prior 
% variance to zero) produces the smallest reduction (greatest increase) 
% in model evidence. This procedure is repeated until all 8 parameters 
% are retained in the best model or there are no more parameters to 
% consider. When several DCMs are optimized together (as in group studies), 
% they are checked to ensure the same free parameters have been specified 
% and the log-evidences are pooled in a fixed effects fashion.
% 
% This application of post hoc optimization assumes the DCMs that are 
% optimized are the same model of different data. Normally, this would be 
% a full model, in the sense of having the maximum number of free 
% parameters, such that the set of reduced models is as large as possible. 
% In contrast spm_dcm_search operates on different DCMs of the same data 
% to identify the best model
% 
% The outputs of this routine are graphics reporting the model reduction 
% (optimisation) and a DCM_opt_??? for every input DCM that contains the 
% reduced conditional parameters estimates (for simplicity, the original 
% kernels and predicted states are retained). The structural and function 
% (spectral embedding) graphs are based on Bayesian parameter averages 
% over multiple DCMs.
%__________________________________________________________________________
% Copyright (C) 2010-2011 Wellcome Trust Centre for Neuroimaging
 
% Karl Friston
% $Id: spm_dcm_post_hoc.m 4307 2011-04-14 16:45:43Z guillaume $
 
% get filenames
%--------------------------------------------------------------------------
try
    P;
catch
    [P, sts] = spm_select([1 Inf],'^DCM.*\.mat$','Select DCM*.mat files');
    if ~sts, return; end
end
 
if ischar(P),   P = cellstr(P); end
if isstruct(P), P = {P}; end
N = numel(P);
 
%-Check models are compatible in terms of their prior variances
%==========================================================================
for j = 1:N
    
    % get prior covariances
    %----------------------------------------------------------------------
    try, load(P{j}); catch, DCM = P{j}; end
    pC   = diag(DCM.M.pC);
    
    
    % and compare it with the first model
    %----------------------------------------------------------------------
    if j == 1
        C = pC;
    else
        if any(~pC - ~C)
            fprintf('Please check model %i for compatibility',j)
            return
        end
    end
end
 
% find free coupling parameters
%--------------------------------------------------------------------------
k     = spm_fieldindices(DCM.Ep,'A','B','D');
k     = k(~~C(k));
n     = length(C);


% If there are too many find those with the least evidence
%--------------------------------------------------------------------------
if length(k) > 16
    
    %-Loop through DCMs and free parameters and get log-evidences
    %----------------------------------------------------------------------
    for j = 1:N
        
        try, load(P{j}); catch, DCM = P{j}; end
        
        % Get priors and posteriors
        % -----------------------------------------------------------------
        qE    = DCM.Ep;
        qC    = DCM.Cp;
        pE    = DCM.M.pE;
        pC    = DCM.M.pC;
        
        % model search over new prior without the i-th parameter
        % -----------------------------------------------------------------
        for i = 1:length(k)
            R      = speye(n,n) - sparse(k(i),k(i),1,n,n);
            Z(i,j) = spm_log_evidence(qE,qC,pE,pC,pE,R*pC*R);
        end
    end
    
    % find parameters with the least evidence
    %----------------------------------------------------------------------
    Z      = sum(Z,2);
    [Z i]  = sort(-Z);
    k      = k(i(1:8));
    
    % flag a greedy search
    %----------------------------------------------------------------------
    repeat = 1;
    
else
    repeat = 0;
end

% Create model space in terms of free parameter indices
%--------------------------------------------------------------------------
K     = spm_perm_mtx(length(k));
 
%-Loop through DCMs and models and get log-evidences
%==========================================================================
for j = 1:N
    
    try, load(P{j}); catch, DCM = P{j}; end
    fprintf('\nsearching (%i): 00%%',j)
    
    % Get priors and posteriors
    % ---------------------------------------------------------------------
    qE    = DCM.Ep;
    qC    = DCM.Cp;
    pE    = DCM.M.pE;
    pC    = DCM.M.pC;
    
    % model search over new prior (covariance)
    % ---------------------------------------------------------------------
    for i = 1:length(K)
        R      = speye(n,n) - sparse(k,k,K(i,:),n,n);
        G(i,j) = spm_log_evidence(qE,qC,pE,pC,pE,R*pC*R);
    
        fprintf('\b\b\b\b')
        fprintf('%-3.0f%%',i*100/length(K))
    end 
end
fprintf('\n')
 
% Pooled model evidence
% =========================================================================
S     = sum(G,2);
S     = S - S(end);
p     = exp(S - max(S));
p     = p/sum(p);
 
% Get selected model (prior covariance)
%--------------------------------------------------------------------------
[q i] = max(p);
R     = speye(n,n) - sparse(k,k,K(i,:),n,n);
rC    = R*pC*R;
 
% Continue greedy search if any parameters have been eliminated
%--------------------------------------------------------------------------
nelim  = full(sum(K(i,:)));
repeat = repeat & nelim;
if repeat
    fprintf('%i paramters eliminated in this step\n',nelim)
end
 
% Show results
% -------------------------------------------------------------------------
spm_figure('Getwin','Graphics');
spm_figure('Clear');
 
subplot(2,2,1)
if length(K) > 32, plot(S,'k'), else, bar(S,'c'), end
title('log-posterior','FontSize',16)
xlabel('model','FontSize',12)
ylabel('log-probability','FontSize',12)
axis square
 
subplot(2,2,2)
if length(K) > 32, plot(p,'k'), else, bar(p,'r'), end
title('model posterior','FontSize',16)
xlabel('model','FontSize',12)
ylabel('probability','FontSize',12)
axis square
 
% conditional estimates of selected model
% =========================================================================
Eq    = 0;
Pq    = 0;
for j = 1:N
    
    % Get priors and posteriors
    % ---------------------------------------------------------------------
    try, load(P{j}); catch, DCM = P{j}; end
    
    qE   = DCM.Ep;
    qC   = DCM.Cp;
    pE   = DCM.M.pE;
    pC   = DCM.M.pC;
    
    % Get posterior of selected model - rC
    % ---------------------------------------------------------------------
    [F Ep Cp] = spm_log_evidence(qE,qC,pE,pC,pE,rC);
    
    % Bayesian parameter average
    %----------------------------------------------------------------------
    Pp   = spm_inv(Cp);
    Eq   = Eq + Pp*spm_vec(Ep);
    Pq   = Pq + Pp;
 
 
    % Put reduced conditional estimates in DCM
    % =====================================================================
 
    % Bayesian inference and variance
    %----------------------------------------------------------------------
    sw   = warning('off','SPM:negativeVariance');
    Pp   = spm_unvec(1 - spm_Ncdf(DCM.T,abs(spm_vec(Ep)),diag(Cp)),Ep);
    Vp   = spm_unvec(diag(Cp),Ep);
    warning(sw);
    
    
    % Store parameter estimates
    %----------------------------------------------------------------------
    DCM.M.pC = rC;
    DCM.Ep   = Ep;
    DCM.Cp   = Cp;
    DCM.Pp   = Pp;
    DCM.Vp   = Vp;   
    
    % Save approximations to model evidence: negative free energy, AIC, BIC
    %----------------------------------------------------------------------
    evidence = spm_dcm_evidence(DCM);
    DCM.F    = F;
    DCM.AIC  = evidence.aic_overall;
    DCM.BIC  = evidence.bic_overall;
    
    
    %-Save optimised DCM
    %======================================================================
    try
        [pth, name] = fileparts(P{j});
        if ~strncmp(name,'DCM_opt_',8)
            name = ['DCM_opt_' name(5:end) '.mat'];
        end
        P{j} = fullfile(pth,name);
    catch
        P{j} = fullfile(pwd,['DCM_opt_' date '.mat']);
    end
    if spm_check_version('matlab','7') >= 0
        save(P{j},'-V6','DCM','F','Ep','Cp');
    else
        save(P{j},'DCM','F','Ep','Cp');
    end
    
end
 
% Show full and reduced conditional estimates (for last DCM)
%--------------------------------------------------------------------------
spm_figure('Getwin','Graphics');
 
i   = spm_fieldindices(DCM.Ep,'A','B','C','D');
qE  = spm_vec(qE);
Ep  = spm_vec(Ep);
 
subplot(2,2,3)
spm_plot_ci(qE(i),qC(i,i)), hold on
title('MAP connections (full)','FontSize',16)
axis square
a   = axis;
 
subplot(2,2,4)
spm_plot_ci(Ep(i),Cp(i,i)), hold off
title('MAP connections (reduced)','FontSize',16)
axis square
axis(a)
drawnow
 
% repeat optimisation (greedy search) if necessary
%--------------------------------------------------------------------------
if repeat
    spm_dcm_post_hoc(P);
    return
end
 
 
% Show structural and functional graphs
%--------------------------------------------------------------------------
spm_figure('Getwin','Graph'); clf
 
% Bayesian parameter average
%--------------------------------------------------------------------------
Cq  = spm_inv(Pq);
Eq  = Cq*Eq;
Eq  = spm_unvec(Eq,pE);
A   = Eq.A + sum(Eq.B,3) + sum(Eq.D,3);
 
spm_dcm_graph(DCM.xY,A)

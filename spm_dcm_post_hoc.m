function spm_dcm_post_hoc(P,fun)
% Post hoc optimisation of DCMs (under Laplace approximation)
% FORMAT spm_dcm_post_hoc(P,[fun])
%
% P         - character/cell array of DCM filenames
%           - or cell array of DCM structures
%
% fun       - optional family definition function: k = fun(A,B,C,D)
%             k = 1,2,...,K for K families or proper subsets of a partition
%             of model space, as function of the adjacency matrices: e.g.,
%
%             fun = inline('any(spm_vec(B(:,:,2))) + 1','A','B','C','D')
%
%             return 1 if there are no bilinear parameters for the 2nd
%             bilinear effect and 2 if there are. fun should be an inline
%             function or script. NB: Model posteriors over families with
%             and without free parameters (in A,B,C and D) are evaluated
%             automatically and saved in DCM_bpa (DCM.Pp)
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
% to identify the best model, after inverting the full(est) model
%
% The outputs of this routine are graphics reporting the model reduction
% (optimisation) and a DCM_opt_??? for every specified DCM that contains
% reduced conditional parameters estimates (for simplicity, the original
% kernels and predicted states are retained). The structural and function
% (spectral embedding) graphs are based on Bayesian parameter averages
% over multiple DCMs, which are stored in DCM_bpa.mat. This DCM also
% contains the posterior probability of models partitioned according to
% whether a particular parameter exists or not:
%
% DCM.Pp     -  Model posterior over parameters (with and without)
% DCM.Ep     -  Bayesian parameter average under selected model
% DCM.Pf     -  Model posteriors over user specified families
% DCM.fun    -  User-specified family definition function
% DCM.files  -  List of DCM files used for Bayesian averaging
%
%__________________________________________________________________________
% Copyright (C) 2010-2011 Wellcome Trust Centre for Neuroimaging

% Karl Friston
% $Id: spm_dcm_post_hoc.m 4711 2012-04-10 13:20:39Z karl $

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

% check family definition functions
%--------------------------------------------------------------------------
if nargin < 2; fun = {}; Pf = 0; end

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
            fprintf('Please check model %i for compatibility ',j)
            return
        end
    end
end

% find free coupling parameters
%--------------------------------------------------------------------------
k     = spm_fieldindices(DCM.Ep,'A','B','C','D');
k     = k(~~C(k));
n     = length(C);


% If there are too many find those with the least evidence
%--------------------------------------------------------------------------
if length(k) > 16
    
    %-Loop through DCMs and free parameters and get log-evidences
    %----------------------------------------------------------------------
    for j = 1:N
        
        try, load(P{j}); catch, DCM = P{j}; end
        fprintf('\ninitial search (%i): 00%%',j)
        
        % Get priors and posteriors
        % -----------------------------------------------------------------
        qE    = DCM.Ep;
        qC    = DCM.Cp;
        pE    = DCM.M.pE;
        pC    = DCM.M.pC;
        
        % Remove (a priori) null space
        % -----------------------------------------------------------------
        U     = spm_svd(pC);
        qE    = U'*spm_vec(qE);
        pE    = U'*spm_vec(pE);
        qC    = U'*qC*U;
        pC    = U'*pC*U;
        
        % model search over new prior without the i-th parameter
        % -----------------------------------------------------------------
        for i = 1:length(k)
            R      = speye(n,n) - sparse(k(i),k(i),1,n,n);
            R      = U'*R*U;
            Z(i,j) = spm_log_evidence(qE,qC,pE,pC,pE,R*pC*R);
            
            fprintf('\b\b\b\b')
            fprintf('%-3.0f%%',i*100/length(k))
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
    
elseif isempty(k)
    
    sprintf('\n There are no free parameters in this model\n')
    return
    
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
    
    % Remove (a priori) null space
    % -----------------------------------------------------------------
    U     = spm_svd(pC);
    qE    = U'*spm_vec(qE);
    pE    = U'*spm_vec(pE);
    qC    = U'*qC*U;
    pC    = U'*pC*U;
    
    % model search over new prior (covariance)
    % ---------------------------------------------------------------------
    for i = 1:length(K)
        R      = speye(n,n) - sparse(k,k,K(i,:),n,n);
        R      = U'*R*U;
        G(i,j) = spm_log_evidence(qE,qC,pE,pC,pE,R*pC*R);
        
        fprintf('\b\b\b\b')
        fprintf('%-3.0f%%',i*100/length(K))
    end
end
fprintf('\n')

% Pooled model evidence
%==========================================================================
S     = sum(G,2);
S     = S - S(end);
p     = exp(S - max(S));
p     = p/sum(p);

% Inference over families (one family per coupling parameter)
%==========================================================================
pE    = DCM.M.pE;
pC    = DCM.M.pC;
for i = 1:length(k)
    Pk(1,i) = mean(p(find(~K(:,i))));
    Pk(2,i) = mean(p(find( K(:,i))));
end
C     = double(~~C);
Pk    = Pk(1,:)./sum(Pk);
Pn    = C;
Pn(k) = Pk;
Pk    = spm_unvec(Pn,pE);

% Inference over families (specified by fun)
%==========================================================================
if ~isempty(fun)
    for i = 1:length(K)
        Pn     = C;
        Pn(find(K(i,:))) = 0;
        Pn     = spm_unvec(Pn,pE);
        Kf(i)  = fun(Pn.A,Pn.B,Pn.C,Pn.D);
    end
    for i = 1:max(Kf)
        Pf(i) = mean(p(find(Kf == i)));
    end
    Pf    = Pf/sum(Pf);
end



% Get selected model (prior covariance)
%==========================================================================
[q i] = max(p);
R     = speye(n,n) - sparse(k,k,K(i,:),n,n);
rC    = R*pC*R;

% record pruned parameters
%--------------------------------------------------------------------------
R     = spm_unvec(diag(R),pE);
R.a   = DCM.a & R.A;
R.b   = DCM.b & R.B;
R.c   = DCM.c & R.C;
R.d   = DCM.d & R.D;


% Continue greedy search if any parameters have been eliminated
%--------------------------------------------------------------------------
nelim  = full(sum(K(i,:)));
repeat = repeat & nelim;
if repeat
    fprintf('%i out of %i free parameters removed \n',nelim,sum(C))
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
drawnow

% conditional estimates of selected model for each data set
%==========================================================================
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
    
    % Bayesian parameter average (for selected model)
    %----------------------------------------------------------------------
    Pp   = spm_inv(Cp);
    Eq   = Eq + Pp*spm_vec(Ep);
    Pq   = Pq + Pp;
    
    
    % Put reduced conditional estimates in DCM
    %======================================================================
    
    % Bayesian inference and variance
    %----------------------------------------------------------------------
    T    = full(spm_vec(pE)) + DCM.T;
    sw   = warning('off','SPM:negativeVariance');
    Pp   = spm_unvec(1 - spm_Ncdf(T,abs(spm_vec(Ep)),diag(Cp)),Ep);
    Vp   = spm_unvec(diag(Cp),Ep);
    warning(sw);
    
    % Store parameter estimates
    %----------------------------------------------------------------------
    DCM.M.pC    = rC;
    DCM.Ep      = Ep;
    DCM.Cp      = Cp;
    DCM.Pp      = Pp;
    DCM.Vp      = Vp;
    
    % and in DEM format
    %----------------------------------------------------------------------
    DCM.qP.P{1} = Ep;
    DCM.qP.C    = Cp;
    DCM.qP.V{1} = spm_unvec(diag(Cp),Ep);
    
    % and in prior constraints fields
    %----------------------------------------------------------------------
    DCM.a       = R.a;
    DCM.b       = R.b;
    DCM.c       = R.c;
    DCM.d       = R.d;
    
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
    save(P{j},'DCM','F','Ep','Cp', spm_get_defaults('mat.format'));
    
end


% repeat optimisation (greedy search) if necessary
%==========================================================================
if repeat
    save('tmp','P'),  clear, load tmp;
    spm_dcm_post_hoc(P);
    return
else
    try, delete('tmp.mat'), end
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


% Show structural and functional graphs
%--------------------------------------------------------------------------
spm_figure('Getwin','Parameter posterior'); clf

% Bayesian parameter average
%--------------------------------------------------------------------------
Cq  = spm_inv(Pq);
Eq  = Cq*Eq;
Eq  = spm_unvec(Eq,pE);

spm_dcm_graph(DCM.xY,Eq.A)

% Show coupling matrices
%--------------------------------------------------------------------------
spm_figure('Getwin','Bayesian parameter average (selected model)'); clf
spm_dcm_fmri_image(Eq)

spm_figure('Getwin','Model posterior (over parameters)'); clf
spm_dcm_fmri_image(Pk)

if ~isempty(fun)
    spm_figure('Getwin','Model posterior (over families)'); clf
    subplot(2,1,1)
    bar(Pf)
    xlabel('familiy')
    title('Model posterior (over families)','FontSize',16)
    axis square
end

%-Save Bayesian parameter average and family-wise model inference
%==========================================================================

% Get original (first) DCM
% -------------------------------------------------------------------------
try, load(P{1}); catch, DCM = P{1}; end

DCM.Pp    = Pk;       % Model posterior over parameters (with and without)
DCM.Ep    = Eq;       % Bayesian parameter average under selected model
DCM.fun   = fun;      % user-specified family definition function
DCM.files = P;        % list of DCM files used for Bayesian averaging
DCM.Pf    = Pf;       % Model posteriors over user specified families

% and save as DCM_bpa
% -------------------------------------------------------------------------
try
    [pth, name] = fileparts(P{1});
    name = ['DCM_bpa.mat'];
    name = fullfile(pth,name);
catch
    name = fullfile(pwd,['DCM_bpa.mat']);
end
save(name,'DCM', spm_get_defaults('mat.format'));

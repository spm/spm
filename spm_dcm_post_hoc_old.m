function DCM = spm_dcm_post_hoc_old(P,fun,varargin)
% Post hoc optimisation of DCMs (under the Laplace approximation)
% FORMAT DCM = spm_dcm_post_hoc(P,fun,field,field,...)
%
% P         - character/cell array of DCM filenames
%           - or cell array of DCM structures; where
%  DCM.M.pE - prior expectation (with parameters in pE.A, pE.B and pE.C)
%  DCM.M.pC - prior covariance
%  DCM.Ep   - posterior expectations
%  DCM.Cp   - posterior covariance
%
% fun       - optional family definition function: k = fun(A,B,C)
%             k = 1,2,...,K for K families or proper subsets of a partition
%             of model space - a function of the adjacency matrices: e.g.,
%
%             fun = @(A,B,C) any(spm_vec(B(:,:,2))) + 1;
%
%             returns 1 if there are no bilinear parameters for the 2nd
%             bilinear effect and 2 if there are. fun should be an inline
%             function or script. NB: Model posteriors over families with
%             and without free parameters (in A,B,C and D) are evaluated
%             automatically and saved in DCM_BPA (DCM.Pp)
%
% field     - the field nsmes of the parameters in the structure pE and Ep 
%             that are to be inlcudied in Baysian model reduction.
%             The default is the cell array 'A','B','C'
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
% kernels and predicted states are retained). The structural and functional
% (spectral embedding) graphs are based on Bayesian parameter averages
% over multiple DCMs, which are stored in DCM_BPA.mat. This DCM also
% contains the posterior probability of models partitioned according to
% whether a particular parameter exists or not:
%
% DCM.Pp     -  Model posterior (with and without each parameter)
% DCM.Ep     -  Bayesian parameter average under selected model
% DCM.Cp     -  Bayesian parameter covariance under selected model
% DCM.Pf     -  Model posteriors over user specified families
% DCM.fun    -  User-specified family definition function
% DCM.files  -  List of DCM files used for Bayesian averaging

% See alos: spm_dcm_search.m
%
%__________________________________________________________________________
% Copyright (C) 2010-2012 Wellcome Trust Centre for Neuroimaging

% Karl Friston
% $Id: spm_dcm_post_hoc_old.m 6724 2016-02-19 19:13:07Z karl $


% number of parameters to consider before invoking greedy search
%--------------------------------------------------------------------------
nmax = 8;

% Get filenames
%--------------------------------------------------------------------------
if nargin < 1
    [P, sts] = spm_select([1 Inf],'^DCM.*\.mat$','Select DCM*.mat files');
    if ~sts, return; end
end

if ischar(P),   P = cellstr(P); end
if isstruct(P), P = {P};        end
N   = numel(P);
TOL = exp(-16);

% Check family definition functions
%--------------------------------------------------------------------------
if nargin < 2; fun = {}; Pf = 0; end
if ~numel(fun); Pf = 0;          end

% Check fields of parameter stucture
%--------------------------------------------------------------------------
if nargin < 3
    field = {'A','B','C'};
else
    field = varargin;
end

%-Check models are compatible in terms of their prior variances
%==========================================================================
for j = 1:N
    
    % Get prior covariances
    %----------------------------------------------------------------------
    try, load(P{j});          catch, DCM = P{j};              end
    try, pC = diag(DCM.M.pC); catch, pC  = spm_vec(DCM.M.pC); end
    
    % and compare it with the first model
    %----------------------------------------------------------------------
    if j == 1
        C = pC;
    else
        if any(xor(pC,C))
            fprintf('Please check model %i for compatibility.\n',j)
            return
        end
    end
end


%-Greedy search (GS) - eliminating parameters in a top down fashion
%==========================================================================

% Accumulated reduction vector (C)
%--------------------------------------------------------------------------
C  = logical(C);
GS = 1;
while GS
    
    % Find free coupling parameters
    %----------------------------------------------------------------------
    k = spm_fieldindices(DCM.Ep,field{:});
    k = k(C(k));
    nparam = length(k);
    
    
    % If there are too many find those with the least evidence
    %----------------------------------------------------------------------
    if nparam > nmax
        
        % Loop through DCMs and free parameters and get log-evidences
        %------------------------------------------------------------------
        Z     = [];
        for j = 1:N
            
            try, load(P{j}); catch, DCM = P{j}; end
            if isstruct(DCM.M.pC), DCM.M.pC = diag(spm_vec(DCM.M.pC)); end
            
            fprintf('\ninitial search (%i): 00%%',j)
            
            % Get priors and posteriors
            %--------------------------------------------------------------
            qE    = DCM.Ep;
            qC    = DCM.Cp;
            pE    = DCM.M.pE;
            pC    = DCM.M.pC;
            
            % Remove (a priori) null space
            %--------------------------------------------------------------
            U     = spm_svd(pC);
            qE    = U'*spm_vec(qE);
            pE    = U'*spm_vec(pE);
            qC    = U'*qC*U;
            pC    = U'*pC*U;
            
            % Model search over new prior without the i-th parameter
            %--------------------------------------------------------------
            for i = 1:length(k)
                r      = C; r(k(i)) = 0;
                R      = U(r,:)'*U(r,:);
                Z(i,j) = spm_log_evidence(qE,qC,pE,pC,pE,R*pC*R);
                
                fprintf('\b\b\b\b')
                fprintf('%-3.0f%%',i*100/length(k))
            end
        end
        
        % Find parameters with the least evidence
        %------------------------------------------------------------------
        Z          = sum(Z,2);
        [unused,i] = sort(-Z);
        k          = k(i(1:8));
        
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
    
    %-Loop through DCMs and models and get log-evidences
    %======================================================================
    G     = [];
    for j = 1:N
        
        try, load(P{j}); catch, DCM = P{j}; end
        if isstruct(DCM.M.pC), DCM.M.pC = diag(spm_vec(DCM.M.pC)); end
        fprintf('\nsearching (%i): 00%%',j)
        
        % Get priors and posteriors
        %------------------------------------------------------------------
        qE    = DCM.Ep;
        qC    = DCM.Cp;
        pE    = DCM.M.pE;
        pC    = DCM.M.pC;
        
        % Remove (a priori) null space
        %------------------------------------------------------------------
        U     = spm_svd(pC);
        qE    = U'*spm_vec(qE);
        pE    = U'*spm_vec(pE);
        qC    = U'*qC*U;
        pC    = U'*pC*U;
        
        % Model search over new prior (covariance)
        %------------------------------------------------------------------
        for i = 1:length(K)
            r      = C; r(k(K(i,:))) = 0;
            R      = U(r,:)'*U(r,:);
            G(i,j) = spm_log_evidence(qE,qC,pE,pC,pE,R*pC*R);
            
            fprintf('\b\b\b\b')
            fprintf('%-3.0f%%',i*100/length(K))
        end
    end
    fprintf('\n')
    
    % Pooled model evidence
    %----------------------------------------------------------------------
    S     = sum(G,2);
    S     = S - S(end);
    p     = exp(S - max(S));
    p     = p/sum(p);
    
    %-Get selected model and prune redundant parameters
    %======================================================================
    [unused, i]  = max(p);
    C(k(K(i,:))) = 0;
    
    % Continue greedy search if any parameters have been eliminated
    %----------------------------------------------------------------------
    nelim  = full(sum(K(i,:)));
    GS     = GS & nelim;
    
    % Show results
    % ---------------------------------------------------------------------
    spm_figure('Getwin','Graphics'); clf    
    fprintf('%i out of %i free parameters removed \n',nelim,nparam)
    
    subplot(3,2,1)
    if length(K) > 32, plot(S,'k'), else, bar(S,'c'), end
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
pE    = DCM.M.pE;
pC    = DCM.M.pC;
for i = 1:length(k)
    Pk(1,i) = mean(p(~K(:,i)));
    Pk(2,i) = mean(p( K(:,i)));
end
Pk    = Pk(1,:)./sum(Pk);
Pn    = full(C);
Pn(k) = Pk;
Pk    = spm_unvec(Pn,pE);


%-Inference over families (specified by fun)
%==========================================================================
if ~isempty(fun)
    for i = 1:length(K)
        Pn     = full(C);
        Pn(K(i,:)) = 0;
        Pn     = spm_unvec(Pn,pE);
        try
            Kf(i)     = fun(Pn.A,Pn.B,Pn.C);
        catch
            try
                Kf(i) = fun(Pn.A,Pn.B);
            catch
                Kf(i) = fun(Pn.A);
            end
        end
            
    end
    for i = 1:max(Kf)
        Pf(i) = mean(p(Kf == i));
    end
    Pf(isnan(Pf)) = 0;
    Pf    = Pf/sum(Pf);
end


%-Conditional estimates of selected model for each data set
%==========================================================================

% Selected model (reduced prior covariance);
%--------------------------------------------------------------------------
rC      = diag(C)*pC*diag(C);

% Record pruned parameters
%--------------------------------------------------------------------------
R       = spm_unvec(full(C),pE);
try
    R.a = DCM.a & R.A;
    R.b = DCM.b & R.B;
    R.c = DCM.c & R.C;
    R.d = DCM.d & R.D;
end

EQ      = 0;
PQ      = 0;
Eq      = 0;
Pq      = 0;
for j = 1:N
    
    % Get priors and posteriors
    %----------------------------------------------------------------------
    try, load(P{j}); catch, DCM = P{j}; end
    
    qE  = DCM.Ep;
    qC  = DCM.Cp;
    pE  = DCM.M.pE;
    pC  = DCM.M.pC;
    
    % Get posterior of selected model - rC
    %----------------------------------------------------------------------
    [F,Ep,Cp] = spm_log_evidence_reduce(qE,qC,pE,pC,pE,rC);
    
    % Bayesian parameter average (for full and reduced selected model)
    %----------------------------------------------------------------------
    PP  = spm_inv(qC,TOL);
    Pp  = spm_inv(Cp,TOL);
    EQ  = EQ + PP*spm_vec(qE);
    Eq  = Eq + Pp*spm_vec(Ep);
    PQ  = PQ + PP;
    Pq  = Pq + Pp;
    
    
    %-Put reduced conditional estimates in DCM
    %======================================================================
    
    % Bayesian inference and variance
    %----------------------------------------------------------------------
    try
        T = full(spm_vec(pE)) + DCM.T;
    catch
        T = full(spm_vec(pE));
    end
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
    try
        DCM.a   = R.a;
        DCM.b   = R.b;
        DCM.c   = R.c;
        DCM.d   = R.d;
    end  
    
    % approximations to model evidence: negative free energy, AIC, BIC
    %------------------------------------------------------------------
    DCM.F    = F;
    try
        evidence = spm_dcm_evidence(DCM);
        DCM.AIC  = evidence.aic_overall;
        DCM.BIC  = evidence.bic_overall;
    end
    
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


%-Bayesian parameter average
%==========================================================================
if isstruct(pC), pC = diag(spm_vec(pC)); end

CQ  = spm_inv(PQ,TOL);
Cq  = spm_inv(Pq,TOL);
EQ  = CQ*EQ;
Eq  = Cq*Eq;
EQ  = spm_unvec(EQ,pE);
Eq  = spm_unvec(Eq,pE);

% Show full and reduced conditional estimates (for Bayesian average)
%--------------------------------------------------------------------------
spm_figure('Getwin','Graphics');

i   = spm_fieldindices(DCM.Ep,'A','B','C');
pE  = spm_vec(pE);
EP  = spm_vec(EQ);
Ep  = spm_vec(Eq);

subplot(3,2,3)
spm_plot_ci(EP(i),CQ(i,i))
title('MAP connections (full)','FontSize',16)
axis square
a   = axis;

subplot(3,2,4)
spm_plot_ci(Ep(i),abs(Cq(i,i)))
title('MAP connections (reduced)','FontSize',16)
axis square
axis(a)

subplot(3,2,5)
bar(EP(i) - pE(i))
xlabel('parameter')
title('MAP minus prior','FontSize',16)
spm_axis tight
axis square

subplot(3,2,6)
bar(Ep(i) - EP(i))
xlabel('parameter')
title('differences in MAP','FontSize',16)
spm_axis tight
axis square
drawnow


% Show structural and functional graphs
%--------------------------------------------------------------------------
if ~nargout
    spm_figure('Getwin','Graph analysis'); clf
    try
        spm_dcm_graph(DCM.xY,Eq.A);
    catch
        try
            spm_dcm_graph(DCM,Eq.A);
        end
    end
end


% Show coupling matrices
%--------------------------------------------------------------------------
if ~nargout && numel(field) == 3;
    
    spm_figure('Getwin','Bayesian parameter average (selected model)'); clf
    spm_dcm_fmri_image(Eq)
    
    spm_figure('Getwin','Model posterior (over parameters)'); clf
    spm_dcm_fmri_image(Pk)

end

% illustrate family-wise inference
%--------------------------------------------------------------------------
if ~isempty(fun) && ~nargout
    
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
DCM.Cp    = Cq;       % Bayesian parameter covariance under selected model
DCM.fun   = fun;      % user-specified family definition function
DCM.files = P;        % list of DCM files used for Bayesian averaging
DCM.Pf    = Pf;       % Model posteriors over user specified families

% and save as DCM_BPA
% -------------------------------------------------------------------------
try
    pth  = fileparts(P{1});
    name = 'DCM_BPA.mat';
    name = fullfile(pth,name);
catch
    name = fullfile(pwd,'DCM_BPA.mat');
end

save(name,'DCM', spm_get_defaults('mat.format'));


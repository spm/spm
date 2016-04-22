function [PEB,P]   = spm_dcm_peb_full(P,M,field)
% Hierarchical (PEB) inversion of DCMs using BMR and VL
% FORMAT [PEB,DCM] = spm_dcm_peb_full(DCM,M,field)
% FORMAT [PEB,DCM] = spm_dcm_peb_full(DCM,X,field)
%
% DCM    - {N [x M]} structure array of DCMs from N subjects
% ------------------------------------------------------------
%     DCM{i}.M.pE - prior expectation of parameters
%     DCM{i}.M.pC - prior covariances of parameters
%     DCM{i}.Ep   - posterior expectations
%     DCM{i}.Cp   - posterior covariance
%
% M.X    - second level design matrix, where X(:,1) = ones(N,1) [default]
% M.pE   - second level prior expectation of parameters
% M.pC   - second level prior covariances of parameters
% M.hE   - second level prior expectation of log precisions
% M.hC   - second level prior covariances of log precisions
% 
% field  - parameter fields in DCM{i}.Ep to optimise [default: {'A','B'}]
%          'All' will invoke all fields. this argument effectively allows 
%          one to specify which parameters constitute random effects.     
% 
% PEB    - hierarchical dynamic model
% -------------------------------------------------------------
%     PEB.Snames - string array of first level model names
%     PEB.Pnames - string array of parameters of interest
%     PEB.Pind   - indices of parameters in spm_vec(DCM{i}.Ep) 
% 
%     PEB.M.X  -   second level (between subject) design matrix
%     PEB.M.W  -   second level (within  subject) design matrix
%     PEB.M.Q  -   precision [components] of second level random effects 
%     PEB.M.pE -   prior expectation of second level parameters
%     PEB.M.pC -   prior covariance  of second level parameters
%     PEB.M.hE -   prior expectation of second level log-precisions
%     PEB.M.hC -   prior covariance  of second level log-precisions
%     PEB.Ep   -   posterior expectation of second level parameters
%     PEB.Eh   -   posterior expectation of second level log-precisions
%     PEB.Cp   -   posterior covariance  of second level parameters
%     PEB.Ch   -   posterior covariance  of second level log-precisions
%     PEB.Ce   -   expected covariance of second level random effects
%     PEB.F    -   free energy of second level model
%
% DCM    - 1st level (reduced) DCM structures with emprical priors
%
%          If DCM is an an (N x M} array, hierarchicial inversion will be
%          applied to each model (i.e., each row) - and PEB will be a 
%          {1 x M} cell array.
%
%--------------------------------------------------------------------------
% This routine inverts a hierarchical DCM using variational Laplace and
% Bayesian model reduction. In essence, it optimises the empirical priors
% over the parameters of a set of first level DCMs, using  second level or
% between subject constraints specified in the design matrix X. This scheme
% is efficient in the sense that it does not require inversion of the first
% level DCMs - it just requires the prior and posterior densities from each
% first level DCMs to compute empirical priors under the implicit
% hierarchical model. The output of this scheme (PEB) can be re-entered
% recursively to invert deep hierarchical models. Furthermore, Bayesian
% model comparison (BMC) can be specified in terms of the empirical
% priors to perform BMC at the group level. Alternatively, subject-specific
% (first level) posterior expectations can be used for classical inference
% in the usual way. Note that these (summary statistics) and  optimal in
% the sense that they have been estimated under empirical (hierarchical) 
% priors.
%
% If called with a single DCM, there are no between subject effects and the
% design matrix is assumed to model mixtures of parameters at the first
% level.
%
% If called with a cell array, each column is assumed to contain 1st level
% DCMs inverted under the same model.
%__________________________________________________________________________
% Copyright (C) 2015-2016 Wellcome Trust Centre for Neuroimaging
 
% Karl Friston
% $Id: spm_dcm_peb_full.m 6778 2016-04-22 11:51:29Z guillaume $
 

% get filenames and set up
%==========================================================================
if ~nargin
    [P, sts] = spm_select([2 Inf],'^DCM.*\.mat$','Select DCM*.mat files');
    if ~sts, return; end
end
if ischar(P),   P = cellstr(P);  end
if isstruct(P), P = {P};         end


% parameter fields
%--------------------------------------------------------------------------
try, load(P{1}); catch, DCM = P{1}; end
if nargin < 2
    M.X   = ones(length(P),1);
end
if isempty(M)
    M.X   = ones(length(P),1);
end
if nargin < 3;
    field = {'A','B'};
end
if strcmpi(field,'all');
    field = fieldnames(DCM.M.pE);
end
if ischar(field)
    field = {field};
end

% repeat for each model (column) if P is an array
%==========================================================================
if size(P,2) > 1
    for i = 1:size(P,2)
        [p,q]   = spm_dcm_peb(P(:,i),M,field);
        PEB(i)  = p;
        P(:,i)  = q;
    end
    return
end


% second level model (ensure first is a constant or main effect)
%--------------------------------------------------------------------------
if ~isstruct(M)
    M = struct('X',M);
end

% use priors from the first level if necessary
%--------------------------------------------------------------------------
if ~isfield(M,'pE')
    M.pE = spm_vec(DCM.M.pE);
end
if ~isfield(M,'pC')
    if isstruct(DCM.M.pC)
        M.pC = diag(spm_vec(DCM.M.pC));
    else
        M.pC = DCM.M.pC;
    end
end
if isstruct(M.pE), M.pE =      spm_vec(M.pE) ; end
if isstruct(M.pC), M.pC = diag(spm_vec(M.pC)); end


% get (first level) densities (summary statistics)
%==========================================================================
q     = spm_find_pC(DCM.M.pC,DCM.M.pE,field);
Pstr  = spm_fieldindices(DCM.M.pE,q);
for i = 1:length(P)
    
    % get first(within subject) level DCM
    %----------------------------------------------------------------------
    try, load(P{i}); catch, DCM = P{i}; end
    
    % posterior densities over all parameters
    %----------------------------------------------------------------------
    if isstruct(DCM.M.pC)
        pC{i} = diag(spm_vec(DCM.M.pC));
    else
        pC{i} = DCM.M.pC;
    end
    pE{i} = spm_vec(DCM.M.pE); 
    qE{i} = spm_vec(DCM.Ep);
    qC{i} = DCM.Cp;
    
    % and free energy of model with full priors
    %----------------------------------------------------------------------
    iF(i) = DCM.F;
    
end

% hierarchical model design and defaults
%==========================================================================
Ns     = numel(P);                    % number of subjects
Np     = length(pE{1});               % number of RFX parameters
Nq     = length(q);                   % number of RFX parameters
Q      = {};                          % precision components


% lower bound on prior precision
%--------------------------------------------------------------------------
pP = spm_inv(M.pC);


% precision components for empirical covariance
%--------------------------------------------------------------------------
if Ns > 1
    OPTION = 'single';
else
    OPTION = 'none';
end
switch OPTION
    
    case {'single'}
        % one between subject precision component
        %------------------------------------------------------------------
        Q{1}          = sparse(Np,Np);
        Q{1}(q,q)     = pP(q,q);
        
    case {'fields'}
        % between subject precision components (one for each field)
        %------------------------------------------------------------------
        for i = 1:length(field)
            j         = spm_fieldindices(DCM.M.pE,field{i});
            j         = find(ismember(j,q));
            Q{i}      = sparse(Np,Np);
            Q{i}(j,j) = pP(j,j);
        end
        
    case {'all'}
        % between subject precision components (one for each parameter)
        %------------------------------------------------------------------
        for i = 1:Nq
            j         = q(i);
            Q{i}(j,j) = sparse(i,i,pP(j,j),Np,Np);
        end
        
    otherwise
end


% priors for empirical expectations
%--------------------------------------------------------------------------
if Ns > 1;
    
    % between-subject design matrices and prior expectations
    %======================================================================
    X     = M.X;
    W     = sparse(q,1:Nq,1,Np,Nq);
    bE    = M.pE(q);
    bC    = M.pC(q,q);

else
    
    % within subject design
    %======================================================================
    if nargin > 1
        W      = sparse(Np,size(M.X,2));
        W(q,:) = M.X;
        X      = 1;
    else
        W = M.pE(q);
        X = 1;
    end
    Nw      = size(W,2);
    try, bE = M.pE; catch, bE = zeros(Nw,1); end
    try, bC = M.bC; catch, bC = eye(Nw,Nw);  end

end


% check for user-specified priors on log precision of second level effects
%--------------------------------------------------------------------------
% Y     = spm_cat(qE)';
% bE    = pinv(X)*Y;
% gE    = -log(trace(pP*cov(Y - X*bE))/Np);

gE    = 0;
gC    = 1;
try, gE = M.hE; end
try, gC = M.hC; end


% prior expectations and precisions of second level parameters
%--------------------------------------------------------------------------
Nx    = size(X,2);                   % number of between subject effects
Nw    = size(W,2);                   % number of within  subject effects
Ng    = length(Q);                   % number of precision components
Nb    = Nw*Nx;                       % number of second level parameters

bX    = speye(Nx,Nx);
bE    = kron(spm_speye(Nx,1),bE);    % prior expectation of group effects
gE    = zeros(Ng,1) + gE;            % prior expectation of log precisions
bC    = kron(bX,bC);                 % prior covariance of group effects
gC    = eye(Ng,Ng)*gC;               % prior covariance of log precisions
bP    = spm_inv(bC);
gP    = spm_inv(gC);


% initialise parameters 
%--------------------------------------------------------------------------
b     = bE;
g     = gE;
p     = [b; g];
ipC   = spm_cat({bP [];
                [] gP});

% variational Laplace
%--------------------------------------------------------------------------
t     = -2;                           % Fisher scoring parameter
for n = 1:32

    % compute prior covariance
    %----------------------------------------------------------------------
    rP    = 0;
    for i = 1:Ng
        rP = rP + exp(g(i))*Q{i};
    end
    rC    = spm_inv(rP);
    
    % update model parameters
    %======================================================================
    
    % Gradient and curvature with respect to free energy
    %----------------------------------------------------------------------
    F     = 0;
    dFdb  = -bP*(b - bE);
    dFdbb = -bP;
    dFdg  = -gP*(g - gE);
    dFdgg = -gP;
    dFdbg = zeros(Nb,Ng);
    for i = 1:Ns
        
        % get empirical prior expectations and reduced 1st level posterior
        %------------------------------------------------------------------
        Xi         = kron(X(i,:),W);
        rE         = Xi*b;
        [Fi,sE,sC] = spm_log_evidence_reduce(qE{i},qC{i},pE{i},pC{i},rE,rC);
        
        % supplement gradients and curvatures
        %------------------------------------------------------------------
        F     = F  + Fi + iF(i);
        dE    = sE - rE;
        dFdb  = dFdb  + Xi'*rP*dE;
        dFdbb = dFdbb + Xi'*(rP*sC*rP - rP)*Xi;
        for j = 1:Ng
            dFdgj      = exp(g(j))*(trace((rC - sC)*Q{j}) - dE'*Q{j}*dE)/2;
            dFdg(j)    = dFdg(j) + dFdgj;
            dFdgg(j,j) = dFdgg(j,j) + dFdgj;
            
            dFdbgj     = exp(g(j))*(Xi - sC*rP*Xi)'*Q{j}*dE;
            dFdbg(:,j) = dFdbg(:,j) + dFdbgj;
            
            for k = 1:Ng
                dFdggj = exp(g(j) + g(k))*(trace((rC*Q{k}*rC - sC*Q{k}*sC)*Q{j})/2 - dE'*Q{k}*sC*Q{j}*dE);
                dFdgg(j,k) = dFdgg(j,k) - dFdggj;
            end
        end
    end

    % Free-energy and update parameters
    %======================================================================
    dFdp  = [dFdb; dFdg];
    dFdpp = spm_cat({dFdbb  dFdbg;
                     dFdbg' dFdgg});
    Cp    = spm_inv(-dFdpp);
    Cb    = spm_inv(-dFdbb);
    Cg    = spm_inv(-dFdgg);
    
    % second level complexity component of free energy
    %----------------------------------------------------------------------
    Fc    = b'*bP*b/2 + g'*gP*g/2 - spm_logdet(ipC*Cp)/2;
    F     = F - Fc;
    
    % best free energy so far
    %----------------------------------------------------------------------
    if n == 1, F0 = F; end
    
    % convergence and regularisation
    %======================================================================
    
    % if F is increasing save current expansion point
    %----------------------------------------------------------------------
    if F >= F0
        
        dF = F - F0;
        F0 = F;
        save('tmp.mat','b','g','F0','dFdb','dFdbb','dFdg','dFdgg');
        
        % decrease regularisation
        %------------------------------------------------------------------
        t  = min(t + 1/4,2);
        
    else
        
        % otherwise, retrieve expansion point and increase regularisation
        %------------------------------------------------------------------
        t  = t - 1;
        load('tmp.mat');
        
    end
    

    % Fisher scoring
    %----------------------------------------------------------------------
    dp      = spm_dx(dFdpp,dFdp,{t});
    [db,dg] = spm_unvec(dp,b,g);
    p       = p + dp;
    b       = b + db;
    g       = g + dg;
    
    % Convergence
    %======================================================================
    fprintf('VL Iteration %-8d: F = %-3.2f dF: %2.4f  [%+2.2f]\n',n,full(F),full(dF),t); 
    if t < -4 || (dF < 1e-4 && n > 4) , break, end
     
end


% assemble output structure
%==========================================================================
for i = 1:Ns
    
    % get first(within subject) level DCM
    %----------------------------------------------------------------------
    try
        load(P{i});
        Sstr{i} = P{i};
    catch
        DCM     = P{i};
        try
            Sstr{i} = DCM.name;
        catch
            Sstr{i} = sprintf('Subject %i',i);
        end
    end
    
    
    % evaluate reduced first level parameters if required
    %----------------------------------------------------------------------
    if nargout > 1
        
        % First level BMR (supplemented with second level complexity)
        %------------------------------------------------------------------
        [Fi,sE,sC] = spm_log_evidence_reduce(qE{i},qC{i},pE{i},pC{i},rE,rC);

        DCM.M.pE = rE;
        DCM.M.pC = rC;
        DCM.Ep   = sE;
        DCM.Cp   = sC;
        DCM.F    = Fi + iF(i) - Fc/Ns;
        
        P{i}     = DCM;
    end
    
end


% second level results
%--------------------------------------------------------------------------
PEB.Snames = Sstr';
PEB.Pnames = Pstr';
PEB.Pind   = q;

PEB.M.X  = X;
PEB.M.W  = W;
PEB.M.pE = bE;
PEB.M.pC = bC;
PEB.M.hE = gE;
PEB.M.hC = gC;
PEB.Ep   = reshape(b,size(W,2),size(X,2));
PEB.Eh   = g;
PEB.Cp   = Cb;
PEB.Ch   = Cg;
PEB.Cph  = Cp;
PEB.Ce   = rC;
PEB.F    = F;

spm_unlink('tmp.mat');

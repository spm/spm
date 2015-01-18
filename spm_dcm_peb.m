function [PEB,P]   = spm_dcm_peb(P,X,field)
% Hierarchical (PEB) inversion of DCMs using BMR and VL
% FORMAT [PEB,DCM] = spm_dcm_peb(DCM,X,field)
%
% DCM    - {N [x M]} structure array of DCMs from N subjects
% ------------------------------------------------------------
%     DCM{i}.M.pE - prior expectation of parameters
%     DCM{i}.M.pC - prior covariances of parameters
%     DCM{i}.M.eE - prior covariances of log precisions
%     DCM{i}.M.eC - prior covariances of log precisions
%     DCM{i}.Ep   - posterior expectations
%     DCM{i}.Cp   - posterior covariance
%
% X      - second level design matrix, where X(:,1) = ones(N,1) [default]
% field  - parameter fields in DCM{i}.Ep to optimise [default: {'A','B'}]
%          'All' will invoke all fields
% 
% PEB    - hierarchical dynamic model
% -------------------------------------------------------------
%     PEB.Snames - string array of first level model names
%     PEB.Pnames - string array of parameters of interest
%     PEB.Pind   - indices of parameters in spm_vec(DCM{i}.Ep) 
% 
%     PEB.SUB  -   first level (within subject) models
%         SUB(i).M  - model and full first level priors
%         SUB(i).pE - empirical prior expectation of first level parameters
%         SUB(i).pC - empirical prior covariance  of first level parameters
%         SUB(i).Ep - empirical posterior expectations
%         SUB(i).Cp - empirical posterior covariance
%         SUB(i).F  - empirical (reduced) free energy
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
% level DCMs – it just requires the prior and posterior densities from each
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
% Copyright (C) 2005 Wellcome Trust Centre for Neuroimaging
 
% Karl Friston
% $Id: spm_dcm_peb.m 6306 2015-01-18 20:50:38Z karl $
 

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
if nargin < 3;
    field  = {'A','B'};
end
if strcmpi(field,'all');
    field = fieldnames(DCM.M.pE);
end

% repeat for each column if P is an array
%--------------------------------------------------------------------------
if isvector(P)
    
    % between subject effects (ensure first is a constant or main effect)
    %----------------------------------------------------------------------
    X(:,1) = ones(length(P),1);
    
else
    
    % between subject effects
    %----------------------------------------------------------------------
    X(:,1) = ones(size(P,1),1);
    
    % loop over models in each column
    %----------------------------------------------------------------------
    for i = 1:size(P,2)
        [p,q]   = spm_dcm_peb(P(:,i),X,field);
        PEB(i)  = p;
        P(:,i)  = q;
    end
    
    return
end



% get (first level) densities (summary statistics)
%==========================================================================
j     = spm_fieldindices(DCM.M.pE,field{:});
q     = j(find(diag(DCM.Cp(j,j))));
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

    
    % select parameters in field
    %----------------------------------------------------------------------
    pE{i} = pE{i}(q); 
    pC{i} = pC{i}(q,q); 
    qE{i} = qE{i}(q); 
    qC{i} = qC{i}(q,q); 
    
    % and free energy of model with full priors
    %----------------------------------------------------------------------
    iF(i) = DCM.F;
    
    if i == 1
            try, gE = DCM.M.eE; catch gE = 0;   end
            try, gC = DCM.M.eC; catch gC = 1/2; end
    end
    
end

% hierarchical dynamic model design and defaults
%==========================================================================
pP     = spm_inv(pC{1});              % lower bound on prior precision
Ns     = numel(P);                    % number of subjects
Np     = length(qC{1});               % number of parameters
Q      = {};                          % precision components


% precision components for empirical covariance
%--------------------------------------------------------------------------
OPTION = 'all';
switch OPTION
    
    case{'single'}
        % one between subject precision component
        %------------------------------------------------------------------
        Q = {pP};
        
    case{'fields'}
        % between subject precision components (one for each field)
        %------------------------------------------------------------------
        for i = 1:length(field)
            j    = spm_fieldindices(DCM.M.pE,field{i});
            j    = find(ismember(q,j));
            Q{i} = sparse(j,j,1,Np,Np);
        end
        
    case{'all'}
        % between subject precision components (one for each parameter)
        %------------------------------------------------------------------
        for i = 1:Np
            p    = pP(i,i);
            if p < exp(16);
                Q{end + 1} = sparse(i,i,pP(i,i),Np,Np);
            end
        end
        
    otherwise
end


% priors for empirical expectations
%--------------------------------------------------------------------------
if Ns > 1;
    
    % between subject design matrices and prior expectations
    %----------------------------------------------------------------------
    W     = speye(Np,Np);
    bE    = pE{1};
    bC    = pC{1};

    
else
    
    % within subject design
    %----------------------------------------------------------------------
    if nargin > 1
        W = X;
        X = 1;
        Q = {};
    else
        W = pE{1};
        X = 1;
    end
    Nw = size(W,2);
    bE = zeros(Nw,1);
    bC = eye(Nw,Nw);

end


% prior expectations and precisions of second level parameters
%--------------------------------------------------------------------------
Nx    = size(X,2);                   % number of between subject effects
Nw    = size(W,2);                   % number of within  subject effects
Ng    = length(Q);                   % number of precision components
Nb    = Nw*Nx;                       % number of second level parameters
bE    = kron(spm_speye(Nx,1),bE);    % prior expectation of group effects
gE    = zeros(Ng,1) + gE;            % prior expectation of log precisions
bC    = kron(eye(Nx,Nx),bC);         % prior covariance of group effects
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
    rP    = pP;
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
        [Fi sE sC] = spm_log_evidence_reduce(qE{i},qC{i},pE{i},pC{i},rE,rC);
        
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

    % Free-energy
    %======================================================================
    dFdp  = [dFdb; dFdg];
    dFdpp = spm_cat({dFdbb  dFdbg;
                     dFdbg' dFdgg});
    Cp    = spm_inv(-dFdpp);
    Cb    = spm_inv(-dFdbb);
    Cg    = spm_inv(-dFdgg);
    F     = F - b'*bP*b/2 - g'*gP*g/2 + spm_logdet(ipC*Cp)/2;
    
    % best free energy so far
    %----------------------------------------------------------------------
    if n == 1, F0 = F; end
    
    % update parameters
    %======================================================================
    
    % if F is increasing save current expansion point
    %----------------------------------------------------------------------
    if F >= F0
        
        dF = F - F0;
        F0 = F;
        save tmp b g dFdb dFdbb dFdg dFdgg
        
        % decrease regularisation
        %------------------------------------------------------------------
        t  = min(t + 1/4,0);
        
    else
        
        % otherwise, retrieve expansion point and increase regularisation
        %------------------------------------------------------------------
        t  = t - 1;
        F  = F0;
        dF = F - F0;
        load tmp
        
    end
    

    % Fisher scoring
    %----------------------------------------------------------------------
    dp      = spm_dx(dFdpp,dFdp,{t});
    [db dg] = spm_unvec(dp,b,g);
    p       = p + dp;
    b       = b + db;
    g       = g + tanh(dg);
    
    % Convergence
    %======================================================================
    fprintf('VL Iteration %-8d: F = %-3.2f dF: %2.4f  [%+2.2f]\n',n,full(F),full(dF),t); 
    if t < -4 || (dF < 1e-2 && ~t && n > 4) , break, end
    
    
end


% assemble output structure
%==========================================================================

% place new (emprical) priors in structures
%--------------------------------------------------------------------------
r    = spm_vec(DCM.M.pE);
r(q) = rE; 
RE   = spm_unvec(r,DCM.M.pE);
if isstruct(DCM.M.pC)
    r      = spm_vec(DCM.M.pC);
    r(q)   = diag(rC); 
    RC     = spm_unvec(r,DCM.M.pE);
else
    r      = spm_inv(DCM.M.pC);
    r(q,q) = r(q,q) + rP - pP;
    RC     = spm_inv(r);
end

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
    
    % get first level posteriors under expected empirical priors
    %----------------------------------------------------------------------
    [Fi sE sC] = spm_log_evidence_reduce(qE{i},qC{i},pE{i},pC{i},rE,rC);
    
    % and save
    %----------------------------------------------------------------------
    SUB(i).M  = DCM.M;
    SUB(i).pE = rE;
    SUB(i).pC = rC;
    SUB(i).Ep = sE;
    SUB(i).Cp = sC;
    SUB(i).F  = Fi + iF(i);
    
    
    % repeat for all parameters
    %----------------------------------------------------------------------
    if nargout > 1
        
        % posterior densities over all parameters
        %------------------------------------------------------------------
        pE{i} = DCM.M.pE;
        pC{i} = DCM.M.pC;
        qE{i} = DCM.Ep;
        qC{i} = DCM.Cp;
        
        % First level BMR
        %------------------------------------------------------------------
        [Fi sE sC] = spm_log_evidence_reduce(qE{i},qC{i},pE{i},pC{i},RE,RC);

        DCM.M.pE = RE;
        DCM.M.pC = RC;
        DCM.Ep   = sE;
        DCM.Cp   = sC;
        DCM.F    = Fi + iF(i);
        
        P{i}     = DCM;
    end
    
end


% second level hierarchical dynamic model
%--------------------------------------------------------------------------
PEB.Snames = Sstr';
PEB.Pnames = Pstr';
PEB.Pind   = q;

PEB.SUB  = SUB;
PEB.M.X  = X;
PEB.M.W  = W;
PEB.M.Q  = Q;
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

try, delete tmp.mat, end




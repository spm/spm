function [PEB,P]   = spm_dcm_peb(P,M,field)
% Hierarchical (PEB) inversion of DCMs using BMR and VL
% FORMAT [PEB,DCM] = spm_dcm_peb(DCM,M,field)
% FORMAT [PEB,DCM] = spm_dcm_peb(DCM,X,field)
%
% DCM    - {N [x M]} structure array of DCMs from N subjects
% -------------------------------------------------------------------------
%     DCM{i}.M.pE	- prior expectation of parameters
%     DCM{i}.M.pC	- prior covariances of parameters
%     DCM{i}.Ep		- posterior expectations
%     DCM{i}.Cp		- posterior covariance
%     DCM{i}.F		- free energy
%
% M.X	   - 2nd-level design matrix: X(:,1) = ones(N,1) [default]
% M.bE	   - 3rd-level prior expectation [default: DCM{1}.M.pE]
% M.bC	   - 3rd-level prior covariance  [default: DCM{1}.M.pC/M.alpha]
% M.pC     - 2nd-level prior covariance  [default: DCM{1}.M.pC/M.beta]
% M.hE	   - 2nd-level prior expectation of log precisions [default: 0]
% M.hC	   - 2nd-level prior covariances of log precisions [default: 1/16]
%
% M.Q      - covariance components: {'single','fields','all','none'}
% M.alpha  - optional scaling to specify M.bC [default = 1]
% M.beta   - optional scaling to specify M.pC [default = 16]
%
% NB: the prior covariance of 2nd-level random effects is:
%            exp(M.hE)*DCM{1}.M.pC/M.beta [default DCM{1}.M.pC/16]
%
% M.Xnames - cell array of names for second level parameters [default: {}]
% 
% field    - parameter fields in DCM{i}.Ep to optimise [default: {'A','B'}]
%          	 'All' will invoke all fields. This argument effectively allows 
%          	 one to specify the parameters that constitute random effects.     
% 
% PEB      - hierarchical dynamic model
% -------------------------------------------------------------------------
%     PEB.Snames 	- string array of first level model names
%     PEB.Pnames 	- string array of parameters of interest
%     PEB.Pind   	- indices of parameters in spm_vec(DCM{i}.Ep) 
%     PEB.Xnames - names of second level parameters
% 
%     PEB.M.X  		- second level (between-subject) design matrix
%     PEB.M.W  		- second level (within-subject) design matrix
%     PEB.M.Q  		- precision [components] of second level random effects 
%     PEB.M.pE 		- prior expectation of second level parameters
%     PEB.M.pC 		- prior covariance of second level parameters
%     PEB.M.hE 		- prior expectation of second level log-precisions
%     PEB.M.hC 		- prior covariance of second level log-precisions
%     PEB.Ep		- posterior expectation of second level parameters
%     PEB.Eh   		- posterior expectation of second level log-precisions
%     PEB.Cp  		- posterior covariance of second level parameters
%     PEB.Ch   		- posterior covariance of second level log-precisions
%     PEB.Ce   		- expected  covariance of second level random effects
%     PEB.F    		- free energy of second level model
%
% DCM    		- 1st level (reduced) DCM structures with empirical priors
%
%          If DCM is an an (N x M} array, hierarchical inversion will be
%          applied to each model (i.e., each row) - and PEB will be a 
%          {1 x M} cell array.
% 
% This routine inverts a hierarchical DCM using variational Laplace and
% Bayesian model reduction. In essence, it optimises the empirical priors
% over the parameters of a set of first level DCMs, using second level or
% between subject constraints specified in the design matrix X. This scheme
% is efficient in the sense that it does not require inversion of the first
% level DCMs - it just requires the prior and posterior densities from each
% first level DCM to compute empirical priors under the implicit
% hierarchical model. The output of this scheme (PEB) can be re-entered
% recursively to invert deep hierarchical models. Furthermore, Bayesian
% model comparison (BMC) can be specified in terms of the empirical priors
% to perform BMC at the group level. Alternatively, subject-specific (first
% level) posterior expectations can be used for classical inference in the
% usual way. Note that these (summary statistics) are optimal in the sense
% that they have been estimated under empirical (hierarchical)  priors. 
%
% If called with a single DCM, there are no between-subject effects and the
% design matrix is assumed to model mixtures of parameters at the first
% level. If called with a cell array, each column is assumed to contain 1st
% level DCMs inverted under the same model.
%
%__________________________________________________________________________
% Copyright (C) 2015-2016 Wellcome Trust Centre for Neuroimaging
 
% Karl Friston
% $Id: spm_dcm_peb.m 6778 2016-04-22 11:51:29Z guillaume $
 

% get filenames and set up
%==========================================================================
if ~nargin
    [P, sts] = spm_select([2 Inf],'^DCM.*\.mat$','Select DCM*.mat files');
    if ~sts, return; end
end
if ischar(P),   P = cellstr(P);  end
if isstruct(P), P = {P};         end

% check for DEM structures
%--------------------------------------------------------------------------
try
    DEM = P;
    P   = spm_dem2dcm(P);
end


% check parameter fields and design matrices
%--------------------------------------------------------------------------
try, load(P{1}); catch, DCM = P{1}; end

if nargin < 2,      M.X = ones(length(P),1); end
if isnumeric(M),    M   = struct('X',M);     end
if ~isfield(M,'X'), M.X = ones(length(P),1); end
if isempty(M.X),    M.X = ones(length(P),1); end

if nargin < 3; field = {'A','B'};  end
if strcmpi(field,'all');  field = fieldnames(DCM.M.pE);end
if ischar(field), field = {field}; end

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


% get (first level) densities (summary statistics)
%==========================================================================
Ns    = numel(P);                               % number of subjects
if isfield(M,'bC') && Ns > 1
    q = spm_find_pC(M.bC,M.bE,field);           % parameter indices
else
    q = spm_find_pC(DCM,field);                 % parameter indices
end
Pstr  = spm_fieldindices(DCM.M.pE,q);           % field names 
Np    = numel(q);                               % number of parameters
if Np == 1
    Pstr = {Pstr}; 
end

for i = 1:Ns
    
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
    pC{i} = pC{i};
    pE{i} = spm_vec(DCM.M.pE);
    qE{i} = spm_vec(DCM.Ep);
    qC{i} = DCM.Cp;
    
    % deal with rank deficient priors
    %----------------------------------------------------------------------
    if i == 1
        PE = pE{i}(q);
        PC = pC{i}(q,q);
        U  = spm_svd(PC);
        Ne = numel(pE{i});
    else
        if numel(pE{i}) ~= Ne
            error('Please ensure all DCMs have the same parameterisation');
        end
    end

    % select parameters in field
    %----------------------------------------------------------------------
    pE{i} = U'*pE{i}(q); 
    pC{i} = U'*pC{i}(q,q)*U; 
    qE{i} = U'*qE{i}(q); 
    qC{i} = U'*qC{i}(q,q)*U;
    
    % shrink posterior to accommodate inefficient inversions
    %----------------------------------------------------------------------
    if Ns > 1
       qC{i} = spm_inv(spm_inv(qC{i}) + spm_inv(pC{i})/16);
    end
    
    % and free energy of model with full priors
    %----------------------------------------------------------------------
    iF(i) = DCM.F;
    
end

% hierarchical model design and defaults
%==========================================================================

% second level model
%--------------------------------------------------------------------------
if  isfield(M,'beta'), alpha  = M.alpha; else, alpha = 1;         end
if  isfield(M,'beta'), beta   = M.beta;  else, beta  = 16;        end
if  isfield(M,'Q'),    OPTION = M.Q;     else, OPTION = 'single'; end
if ~isfield(M,'W'),    M.W    = speye(Np,Np);                     end


% design matrices
%--------------------------------------------------------------------------
if Ns > 1;
    
    % between-subject design matrices and prior expectations
    %======================================================================
    X   = M.X;
    W   = U'*M.W*U;
    
else
    
    % within subject design
    %======================================================================
    OPTION = 'none';
    U      = 1;
    X      = 1;
    W      = M.W;

end

% variable names
%--------------------------------------------------------------------------
if isfield(M,'Xnames')
    Xnames = M.Xnames;
else
    Nx = size(X,2);
    Xnames = cell(1,Nx);
    for i = 1:Nx
        Xnames{i} = sprintf('Covariate %d',i);
    end
end

% get priors (from DCM if necessary) and ensure correct sizes
%--------------------------------------------------------------------------
if isfield(M,'bE')
    M.bE = spm_vec(M.bE);
    if size(M.bE,1) > Np && Ns > 1, M.bE = M.bE(q);   end
else
    M.bE = PE;
end
if isfield(M,'bC')
    if isstruct(M.bC),    M.bC = diag(spm_vec(M.bC)); end
    if size(M.bC,1) > Np && Ns > 1, M.bC = M.bC(q,q); end
else
    M.bC = PC/alpha;
end
if isfield(M,'pC')
    if isstruct(M.pC),    M.pC = diag(spm_vec(M.pC)); end
    if size(M.pC,1) > Np && Ns > 1, M.pC = M.pC(q,q); end
else
    M.pC = PC/beta;
end


% prior precision (pP) and components (Q) for empirical covariance
%--------------------------------------------------------------------------
pQ    = spm_inv(U'*M.pC*U);
rP    = pQ;
Q     = {};
switch OPTION
    
    case{'single'}
        % one between subject precision component
        %------------------------------------------------------------------
        Q = {pQ};
        
    case{'fields'}
        % between subject precision components (one for each field)
        %------------------------------------------------------------------
        for i = 1:length(field)
            j    = spm_fieldindices(DCM.M.pE,field{i});
            j    = find(ismember(q,j));
            Q{i} = sparse(Np,Np);
            Q{i}(j,j) = pQ(j,j);
            Q{i} = U'*Q{i}*U;
        end
        
    case{'all'}
        % between subject precision components (one for each parameter)
        %------------------------------------------------------------------
        for i = 1:Np
            Q{i} = sparse(i,i,pQ(i,i),Np,Np);
            Q{i} = U'*Q{i}*U;
        end
        
    otherwise
end


% number of parameters and effects
%--------------------------------------------------------------------------
Nx    = size(X,2);                  % number of between subject effects
Nw    = size(W,2);                  % number of within  subject effects
Ng    = numel(Q);                   % number of precision components
Nb    = Nw*Nx;                      % number of second level parameters

% check for user-specified priors on log precision of second level effects
%--------------------------------------------------------------------------
gE    = 0;
gC    = 1/16;
try, gE = M.hE; end
try, gC = M.hC; end
try, Xc = M.bX; catch
    
    % adjust (second level) priors for the norm of explanatory variables
    %----------------------------------------------------------------------
    Xc  = diag(size(X,1)./sum(X.^2));
    
end

% prior expectations and precisions of second level parameters
%--------------------------------------------------------------------------
Xe    = spm_speye(Nx,1);
bE    = kron(Xe,U'*M.bE);            % prior expectation of group effects
gE    = zeros(Ng,1) + gE;            % prior expectation of log precision
bC    = kron(Xc,U'*M.bC*U);          % prior covariance of group effects
gC    = eye(Ng,Ng)*gC;               % prior covariance of log precision
bP    = spm_inv(bC);
gP    = spm_inv(gC);

% initialise parameters 
%--------------------------------------------------------------------------
b     = bE;
g     = gE;
ipC   = spm_cat({bP [];
                [] gP});
            
% variational Laplace
%--------------------------------------------------------------------------
t     = -4;                         % Fisher scoring parameter
for n = 1:64

    % compute prior precision (with a lower bound of pQ/4096)
    %----------------------------------------------------------------------
    if Ng > 0
        rP   = pQ/4096;        
        for i = 1:Ng
            rP = rP + exp(g(i))*Q{i};
        end
    end
    rC      = spm_inv(rP);
    
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
        rE{i}      = Xi*b;
        [Fi,sE,sC] = spm_log_evidence_reduce(qE{i},qC{i},pE{i},pC{i},rE{i},rC);
        
        % supplement gradients and curvatures
        %------------------------------------------------------------------
        F     = F  + Fi + iF(i);
        dE    = sE - rE{i};
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
    Fb    = b'*bP*b;
    Fg    = g'*gP*g;
    Fc    = Fb/2 + Fg/2 - spm_logdet(ipC*Cp)/2;
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
    b       = b + db;
    g       = g + dg;
    
    % Convergence
    %======================================================================
    if ~isfield(M,'noplot')
        fprintf('VL Iteration %-8d: F = %-3.2f dF: %2.4f  [%+2.2f]\n',n,full(F),full(dF),t); 
    end
    if t < -4 || (dF < 1e-4 && n > 4)
        fprintf('VL Iteration %-8d: F = %-3.2f dF: %2.4f  [%+2.2f]\n',n,full(F),full(dF),t); 
        break
    end
     
end


% assemble output structure
%==========================================================================
for i = 1:Ns
    
    % get first(within subject) level DCM
    %----------------------------------------------------------------------
    try
        load(P{i},'DCM');
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
        
        % posterior densities over all parameters
        %------------------------------------------------------------------
        if isstruct(DCM.M.pC)
            pC{i} = diag(spm_vec(DCM.M.pC));
        else
            pC{i} = DCM.M.pC;
        end
        pE{i} = DCM.M.pE;
        qE{i} = DCM.Ep;
        qC{i} = DCM.Cp;
        
        % augment empirical priors
        %------------------------------------------------------------------
        RP       = spm_inv(pC{i});
        RP(q,q)  = U*rP*U';
        RC       = spm_inv(RP);
        
        RE       = spm_vec(pE{i});
        RE(q)    = U*rE{i};
        RE       = spm_unvec(RE,pE{i});
       
        % First level BMR (supplemented with second level complexity)
        %------------------------------------------------------------------
        [Fi,sE,sC] = spm_log_evidence_reduce(qE{i},qC{i},pE{i},pC{i},RE,RC);

        DCM.M.pE = RE;
        DCM.M.pC = RC;
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
PEB.Xnames = Xnames;

Ub       = kron(eye(Nx,Nx),U);
PEB.M.X  = X;
PEB.M.W  = W;
PEB.M.U  = U;
PEB.M.pE = kron(Xe,M.bE);
PEB.M.pC = kron(Xc,M.bC);
PEB.M.hE = gE;
PEB.M.hC = gC;
PEB.Ep   = U*reshape(b,Nw,Nx);
PEB.Eh   = g;
PEB.Ch   = Cg;
PEB.Cp   = Ub*Cb*Ub';
PEB.Ce   = U*rC*U';
PEB.F    = F;

spm_unlink('tmp.mat');

% check for DEM structures
%--------------------------------------------------------------------------
try, P   = spm_dem2dcm(DEM,P); end


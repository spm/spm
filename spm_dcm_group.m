function [HDM] = spm_dcm_group(P,X,field)
% hhierarchical inversion of DCMs using Bayesian model reduction and VL
% FORMAT [HDM] = spm_dcm_group(P,X,field)
%
% DCM    - {N} structure array of DCMs from N subjects
% ------------------------------------------------------------
%     DCM{i}.M.pE - prior expectations of M parameters
%     DCM{i}.M.pC - prior covariance
%     DCM{i}.Ep   - posterior expectations
%     DCM{i}.Cp   - posterior covariance
%
% X      - second level design matrix, where X(:,1) = ones(N,1) [default]
% field  - parameter fields in DCM{i}.Ep to optimise [default: {'A','B'}]
%          'All' will invoke all fields
% 
% HDM    - hierarchical dynamic model
% -------------------------------------------------------------
%     HDM.Snames - string array of first level model names
%     HDM.Pnames - string array of parameters of interest
%     HDM.Pind   - indices of parameters in spm_vec(DCM{i}.Ep) 
% 
%     HDM.SUB  -   first level (within subject) models
%         SUB(i).M  - model and full first level priors
%         SUB(i).pE - empirical prior expectation of first level parameters
%         SUB(i).pC - empirical prior covariance  of first level parameters
%         SUB(i).Ep - empirical posterior expectations
%         SUB(i).Cp - empirical posterior covariance
%         SUB(i).F  - empirical (reduced) free energy
%
%     HCM.M.X  -   second level (between subject) design matrix
%     HCM.M.W  -   second level (within  subject) design matrix
%     HDM.M.Q  -   precision [components] of second level random effects 
%     HDM.M.pE -   prior expectation of second level parameters
%     HDM.M.pC -   prior covariance  of second level parameters
%     HDM.M.hE -   prior expectation of second level log-precisions
%     HDM.M.hC -   prior covariance  of second level log-precisions
%     HDM.Ep   -   posterior expectation of second level parameters
%     HDM.Eh   -   posterior expectation of second level log-precisions
%     HDM.Cp   -   posterior covariance  of second level parameters
%     HDM.Ch   -   posterior covariance  of second level log-precisions
%     HDM.Ce   -   expected covariance of second level random effects
%     HDM.F    -   free energy of second level model
%
%--------------------------------------------------------------------------
% This routine inverts a hierarchical DCM using variational Laplace and
% Bayesian model reduction. In essence, it optimises  the empirical priors
% over the parameters of a set of first level DCMs, using  second level or
% between subject constraints specified in the design matrix X.this scheme
% is efficient in the sense that it does not require inversion of the first
% level DCMs – it just requires the prior and posterior densities from each
% first level DCMs  to compute empirical priors under the implicit
% hierarchical model. The output of this scheme (HDM) can be re-entered
% recursively to invert deep hierarchical models. Furthermore, Bayesian
% model comparison (BMC) can be specified in terms of the empirical
% priors to perform BMC at the group level. Alternatively, subject-specific
% (first level) posterior expectations can be used for classical inference
% in the usual way. Note that these (summary statistics) and now optimal in
% the sense that they have been estimated under empirical (hierarchical) 
% priors.
%
% If called with a single DCM, there are no between subject effects and the
% design matrix is assumed to model mixtures of parameters at the first
% level.
%__________________________________________________________________________
% Copyright (C) 2005 Wellcome Trust Centre for Neuroimaging
 
% Karl Friston
% $Id: spm_log_evidence_reduce.m 5394 2013-04-07 14:51:28Z karl $
 
% Compute reduced log-evidence
%==========================================================================

% get filenames and set up
%--------------------------------------------------------------------------
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
if strcmp(lower(field),'all');
    field = fieldnames(DCM.M.pE);
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
    pE{i} = spm_vec(DCM.M.pE); 
    pC{i} = spm_vec(DCM.M.pC);
    qE{i} = spm_vec(DCM.Ep);
    qC{i} = DCM.Cp;
    if isvector(pC{i})
        pC{i} = diag(pC{i});
    end
    
    % select parameters in field
    %----------------------------------------------------------------------
    pE{i} = pE{i}(q); 
    pC{i} = pC{i}(q,q); 
    qE{i} = qE{i}(q); 
    qC{i} = qC{i}(q,q); 
    
    % and free energy of full model
    %----------------------------------------------------------------------
    iF(i) = DCM.F;
    
    
end

% hierarchical dynamic model design and defaults
%==========================================================================
pP     = spm_inv(pC{1});              % lower bound on prior precision
Ns     = numel(P);                    % number of subjects
Np     = length(qC{1});               % number of parameters
X(:,1) = ones(Ns,1);                  % between subject effects

OPTION = 'all';

if Ns > 1;
    
    % between subject design matrices and prior expectations
    %----------------------------------------------------------------------
    W     = speye(Np,Np);
    bE    = pE{1};
    bC    = pC{1};
    
    switch OPTION
        case{'single'}       
            % one between subject precision component
            %--------------------------------------------------------------
            Q = {pP};
            
        case{'fields'}
            % between subject precision components (one for each field)
            %--------------------------------------------------------------
            for i = 1:length(field)
                j    = spm_fieldindices(DCM.M.pE,field{i});
                j    = find(ismember(q,j));
                Q{i} = sparse(j,j,1,Np,Np);
            end
            
        case{'all'}
            % between subject precision components (one for each parameter)
            %--------------------------------------------------------------
            for i = 1:Np
                Q{i} = sparse(i,i,pP(i,i),Np,Np);
            end
            
        otherwise
    end
    
else
    
    % within subject design
    %----------------------------------------------------------------------
    if nargin > 1
        W = X;
        X = 1;
    else
        W = pE{1};
    end
    Nb = size(W,2);
    bE = zeros(Nb,1);
    bC = eye(Nb,Nb);
    Q  = {};
    
end

% prior expectations and precisions of second level parameters
%--------------------------------------------------------------------------
Nx    = size(X,2);                   % number of between subject effects
Ng    = length(Q);                   % number of precision components
Nb    = Np*Nx;                       % number of second level parameters
bE    = kron(spm_speye(Nx,1),bE);    % prior expectation of group effects
gE    = zeros(Ng,1) + 0;             % prior expectation of log precisions
bC    = kron(eye(Nx,Nx),bC);         % prior covariance of group effects
gC    = eye(Ng,Ng)/2;                % prior covariance of log precisions
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
t     = 0;                           % Fisher scoring parameter
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
        F     = F  + Fi;
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
    Cp    = -spm_inv(dFdpp);
    Cb    = -spm_inv(dFdbb);
    Cg    = -spm_inv(dFdgg);
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
        t  = min(t + 1/4,4);
        
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
    dp    = spm_dx(dFdpp,dFdp,{t});
    b     = b + dp(1:Nb);
    g     = g + dp(Nb + 1:end);
    
    % Convergence (1% change in log-evidence)
    %======================================================================
    fprintf('VL Iteration %-8d: F = %-3.2f dF: %2.4f  [%+2.2f]\n',n,full(F),full(dF),t); 
    if dF < 1e-1 && n > 4, break, end
    
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
    
    % get reduced first level posteriors
    %----------------------------------------------------------------------
    [Fi sE sC] = spm_log_evidence_reduce(qE{i},qC{i},pE{i},pC{i},rE,rC);
    
    % and save
    %----------------------------------------------------------------------
    SUB(i).M  = DCM.M;
    SUB(i).pE = rE;
    SUB(i).pC = rC;
    SUB(i).Ep = sE;
    SUB(i).Cp = sC;
    SUB(i).F  = Fi;
    
end


% second level hierarchical dynamic model
%--------------------------------------------------------------------------
HDM.Snames = Sstr';
HDM.Pnames = Pstr';
HDM.Pind   = q;

HDM.SUB  = SUB;
HDM.M.X  = X;
HDM.M.W  = W;
HDM.M.Q  = Q;
HDM.M.pE = bE;
HDM.M.pC = bC;
HDM.M.hE = gE;
HDM.M.hC = gC;
HDM.Ep   = reshape(b,size(W,2),size(X,2));
HDM.Eh   = g;
HDM.Cp   = Cb;
HDM.Ch   = Cg;
HDM.Ce   = rC;
HDM.F    = F;


try, delete tmp.mat, end




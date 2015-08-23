function [DCM,PEB,M,HCM] = spm_dcm_peb_fit(GCM,M,field)
% Bayesian group inversion using empirical Bayes
% FORMAT [DCM,PEB,M] = spm_dcm_peb_fit(DCM,M,field)
%
% DCM    - {N [x M]} structure array of DCMs from N subjects
% ------------------------------------------------------------
%     DCM{i}.M.pE - prior expectation of parameters
%     DCM{i}.M.pC - prior covariances of parameters
%     DCM{i}.Ep   - posterior expectations
%     DCM{i}.Cp   - posterior covariance
%     DCM{i}.FEB  - free energy over empirical Bayes iterations
%
% M.X    - second level design matrix, where X(:,1) = ones(N,1) [default]
% M.pE   - second level prior expectation of parameters
% M.pC   - second level prior covariances of parameters
% M.hE   - second level prior expectation of log precisions
% M.hC   - second level prior covariances of log precisions
%
% field  - parameter fields in DCM{i}.Ep to optimise [default: {'A','B'}]
%          'All' will invoke all fields (i.e. random effects)
%
%
% DCM  - DCM structures inverted with emprical priors
% PEB  - second level model structure
% F    - ssecond level free energy over iterations
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
%--------------------------------------------------------------------------
% This routine performs hierarchical empirical Bayesian inversion of a
% group DCM study. It uses Bayesian model reduction to place second
% (between subject) level constraints on the coordinate descent implicit
% in the inversion of DCMs at the first (within subject) level. In other
% words, at each iteration (or small number of iterations) of the within
% subject inversion, the priors are updated using empirical priors from
% the second level. The free energy of this hierarchical model comprises
% the complexity of group effects plus the sum of free energies from each
% subject - evaluated under the empirical priors  provided by the second
% level.
%
% If called with a cell array, each column is assumed to contain the same
% model of a different subject or dataset, while each row contains
% different models of the same dataset. Bayesian model reduction will be
% applied automatically, after inversion of the full model, which is
% assumed to occupy the first column.
%
% The posterior densities of subject or session specific DCMs are adjusted
% so that they correspond to what would have been obtained under the
% original priors. Effectively, this group inversion is used to suppress
% local minima, prior to inference on group means.
%
% see also: spm_dcm_fit.m; spm_dcm_peb.m; spm_dcm_bmr.m
%__________________________________________________________________________
% Copyright (C) 2005 Wellcome Trust Centre for Neuroimaging

% Karl Friston
% $Id: spm_dcm_peb_fit.m 6532 2015-08-23 13:59:19Z karl $


% set up
%==========================================================================

% Number of subjects (data) and models (of those data)
%--------------------------------------------------------------------------
[Ns,Nm]   = size(GCM);
for i = 1:Ns
    for j = 1:Nm
        if ischar(GCM{i,j})
            model    = load(GCM{i,j});
            GCM{i,j} = model.DCM;
        end
    end
end

% get priors
%--------------------------------------------------------------------------
DCM       = GCM(:,1);
[~,rC,rE] = spm_find_pC(DCM{1});

% priors and parameter fields
%--------------------------------------------------------------------------
if nargin < 2;
    M.X   = ones(Ns,1);
end
if nargin < 3;
    field = 'all'; 
end


% enforce fixed priors at second level
%--------------------------------------------------------------------------
if ~isfield(M,'Q');  M.Q    = 'single'; end
if ~isfield(M,'bE'); M.bE   = rE;       end
if ~isfield(M,'bC'); M.bC   = rC;       end


% reinvert (full) model with initialization; recursively
%==========================================================================
for i = 1:Ns
    DCM{i,1}.M.Nmax = 64;
    DCM{i,1}.M.hE   = 6;
    DCM{i,1}.M.hC   = 1/32;
    try, dipfit{i}  = DCM{i,1}.M.dipfit; end
end
for k = 1:4
    
    % replace spatial model if necessary
    %----------------------------------------------------------------------
    for i = 1:Ns
        try, DCM{i,1}.M.dipfit = dipfit{i}; end
    end
    
    
    % re-initialise and invert the full (first) model
    %----------------------------------------------------------------------
    try, DCM  = spm_dcm_fit(DCM); catch, break;  end
    
     
    % empirical Bayes - over subjects
    %----------------------------------------------------------------------
    [PEB,DCM] = spm_dcm_peb(DCM,M,field);
    
    
    % convergence
    %----------------------------------------------------------------------
    E(k) = PEB.Eh;
    F(k) = PEB.F;
    H(k) = spm_logdet(PEB.Cp);
    
    disp('log precision      : ');disp(E);
    disp('free energy        : ');disp(F);
    disp('conditional entropy: ');disp(H);
    
    if k > 1
        if H(k) > H(k - 1)
            DCM = tmpDCM;
            PEB = tmpPEB;
            break
        else
            tmpDCM = DCM;
            tmpPEB = PEB;
        end
    else
        tmpDCM = DCM;
        tmpPEB = PEB;
    end
    
    % save
    %----------------------------------------------------------------------
    HCM(:,k) = DCM;
    
end

% save second level free energy
%--------------------------------------------------------------------------
for i = 1:Ns
    DCM{i}.EEB = E;
    DCM{i}.FEB = F;
    DCM{i}.HEB = H;
end

% restore original priors
%==========================================================================
DCM = spm_dcm_reduce(DCM,rE,rC);
HCM = spm_dcm_reduce(HCM,rE,rC);
   
% Bayesian model reduction if necessary
%==========================================================================
if Nm > 1
    GCM(:,1) = DCM;
    DCM      = spm_dcm_bmr(GCM,'none');
end

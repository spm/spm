function [DCM,PEB,F] = spm_dcm_peb_fit(GCM,M,field)
% Bayesian ggroup inversion  using empirical Bayes
% FORMAT [DCM,PEB,F] = spm_dcm_peb_fit(DCM,M,field)
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
% subject – evaluated under the empirical priors  provided by the second
% level.
%
% If called with a cell array, each column is assumed to contain the same
% model of a different subject or dataset, while each row contains
% different models of the same dataset. Bayesian model reduction will be
% applied automatically, after inversion of the full model, which is
% assumed to occupy the first column.
%
% see also: spm_dcm_fit.m; spm_dcm_peb.m; spm_dcm_bmr.m
%__________________________________________________________________________
% Copyright (C) 2005 Wellcome Trust Centre for Neuroimaging

% Karl Friston
% $Id: spm_dcm_fit.m 6305 2015-01-17 12:40:51Z karl $


% set up
%==========================================================================

% Number of subjects (data) and models (of those data)
%--------------------------------------------------------------------------
[Ns,Nm]   = size(GCM);
DCM       = GCM(:,1);

% parameter fields
%--------------------------------------------------------------------------
if nargin < 2;
    M.X   = ones(Ns,1);
end
if nargin < 3;
    field = fieldnames(DCM{1,1}.M.pE);
end


% reinvert (full) model with initialization; recursively
%==========================================================================
for i = 1:Ns
    DCM{i,1}.M.Nmax = 8;
    try, dipfit{i}  = DCM{i,1}.M.dipfit; end
end
for k = 1:64
    
    % replace spatial model if necessary
    %----------------------------------------------------------------------
    for i = 1:Ns
        try, DCM{i,1}.M.dipfit = dipfit{i}; end
    end
    
    
    % re-initialise and invert the full (first) model
    %----------------------------------------------------------------------
    try, DCM  = spm_dcm_fit(DCM); catch, break;  end
     
    % empirical Bayes – over subjects
    %----------------------------------------------------------------------
    [PEB,DCM] = spm_dcm_peb(DCM,M,field);
    
    % get intial parameters
    %----------------------------------------------------------------------
    for i = 1:Ns
        % DCM{i,1}.M.pC = DCM{i,1}.M.pC*2;
        DCM{i,1}.M.P  = DCM{i,1}.Ep;
    end
    
    % convergence
    %----------------------------------------------------------------------
    if k == 1
        F(k) = PEB.F;
        dF   = 0;
    else
        dF   = PEB.F - F(k - 1);
        F(k) = PEB.F;
    end
    if dF < 1e-2 && k > 32; break, end
    
end


% Bayesian model reduction if necessary
%==========================================================================
if Nm > 1
    
    % place posteriors in full model (first column)
    %----------------------------------------------------------------------
    GCM(:,1) = DCM;
    
    % place emprical priros in reduced models
    %----------------------------------------------------------------------
    pE    = GCM{1,1}.M.pE;
    pC    = GCM{1,1}.M.pC;
    Np    = size(pC,1);
    for j = 1:Nm
        k     = spm_find_pC(GCM{1,j});
        R     = sparse(k,k,1,Np,Np);
        for i = 1:Ns
            GCM{i,j}.M.pE = spm_unvec(R*spm_vec(pE),pE);
            GCM{i,j}.M.pC = R*pC*R;
        end
    end
    
    % and reduce
    %----------------------------------------------------------------------
    DCM   = spm_dcm_bmr(GCM);
    
end



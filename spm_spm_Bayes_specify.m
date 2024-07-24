function sP = spm_spm_Bayes_specify(SPM)
% Specification of a PEB model for a voxel with empirical priors
%
% SPM - standard SPM structure (see spm_spm.m)
%
% sP(i).P{1}.X - 1st level design matrix
% sP(i).P{1}.C - 1st level prior covariance (see spm_est_non_sphericity.m)
% sP(i).P{2}.X - 2nd level expected values (zeros)
% sP(i).P{2}.C - 2nd level prior covariance (empirical prior)
% sP(i).u      - indices of scans
% sP(i).v      - indices of regressors
%
% ...for each separable partition (e.g. session) of the design i
%__________________________________________________________________________
%
% Creates a structure for a 2-level hierarchical regression model, 
% compatible with spm_PEB.m. The spatial covariance of the betas over 
% voxels is used as an empirical prior for voxel-wise estimation.
%__________________________________________________________________________

% Will Penny
% Copyright (C) 2002-2024 Wellcome Centre for Human Neuroimaging

xX = SPM.xX;

% Number of scans
nScan = size(xX.X,1);

% Number of hyperparameters
try
    nHp       = length(SPM.nscan);
catch
    nHp       = nScan;
    SPM.nscan = nScan;
end

% Number of separable partitions
s = nHp;

% get row u{i} and column v{i}/v0{i} indices for separable designs
%--------------------------------------------------------------------------
if isfield(SPM,'Sess')
    for i = 1:s
         u{i} = SPM.Sess(i).row;
         v{i} = SPM.Sess(i).col;
        v0{i} = xX.iB(i);
    end
else
     u{1} = [1:nScan];
     v{1} = [xX.iH xX.iC]; % Indices of conditions, covariates
    v0{1} = [xX.iB xX.iG]; % Indices of blocks, nuisance variables
end

% cycle over separarable partitions
%--------------------------------------------------------------------------
for i = 1:s

    % Get design X and confounds X0
    %----------------------------------------------------------------------
    fprintf('%-30s\n',sprintf('  ReML Session %i',i));                  %-#
    X     = xX.X(u{i}, v{i});
    X0    = xX.X(u{i},v0{i});
    [m,n] = size(X);

    % add confound in 'filter'
    %----------------------------------------------------------------------
    if isstruct(xX.K)
        X0 = full([X0 xX.K(i).X0]);
    end

    % orthogonalize X w.r.t. X0
    %----------------------------------------------------------------------
    X = X - X0*(pinv(X0)*X);

    % covariance components induced by parameter variations {Q}
    %----------------------------------------------------------------------
    for j = 1:n
        Q{j} = X*sparse(j,j,1,n,n)*X';
    end

    % covariance components induced by error non-sphericity {V}
    %----------------------------------------------------------------------
    Q{n + 1} = SPM.xVi.V(u{i},u{i});

    % ReML covariance component estimation
    %----------------------------------------------------------------------
    [C,h]   = spm_reml(SPM.xVi.CY,X0,Q);

    % check for negative variance components
    %----------------------------------------------------------------------
    h       = abs(h);

    % 2-level model for this partition using prior variances sP(i)
    % treat confounds as fixed (i.e. infinite prior variance)
    %----------------------------------------------------------------------
    n0      = size(X0,2);
    Cb      = blkdiag(diag(h(1:n)),speye(n0,n0)*1e8);
    P{1}.X  = [X X0];
    P{1}.C  = {SPM.xVi.V};
    P{2}.X  = sparse(size(P{1}.X,2),1);
    P{2}.C  = Cb;

    sP(i).P = P;
    sP(i).u = u{:};
    sP(i).v = [v{:} v0{:}];
end

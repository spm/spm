function sP = spm_spm_Bayes_Cb(SPM)

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
    X     = X - X0*(pinv(X0)*X);

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

function varargout = spm_bias_lib(action,varargin)
% Library of functions for Bias correction
%
% Bias correction is performed by optimising a GMM fit to the data.
%
%--------------------------------------------------------------------------
% Basis functions
% ---------------
%
% nbcmp = spm_bias_lib('fwhm2nbcomp', lattice, voxel_size, fwhm)
% basis = spm_bias_lib('dcbasis',     lattice, nb_components)
% prec  = spm_bias_lib('regulariser', mode, lattice, nb_components, voxel_size)
% field = spm_bias_lib('reconstruct', basis, coefficients, ['mult']/'add')
% [TODO] coeff = spm_bias_lib('rescale',     coeff, centre)
%
%--------------------------------------------------------------------------
% Optimisation
% ------------
%
% [g,H]     = spm_bias_lib('derivatives', p, obs, basis, resp, cluster, codes, binvar)
% [ll,bias] = spm_bias_lib('objective',   obs, resp, bias, mean, prec, codes, binvar)
% [TODO] ll = spm_bias_lib('prior',       coeff, precision)
%
%--------------------------------------------------------------------------
% Visualisation
% -------------
%
% spm_bias_lib('plot', 'LB',  lb)
% spm_bias_lib('plot', 'Bias', X, B)
%__________________________________________________________________________
% Copyright (C) 2018-2020 Wellcome Centre for Human Neuroimaging

% $Id: spm_bias_lib.m 7856 2020-05-19 22:54:41Z spm $


switch lower(action)
    case 'fwhm2nbcomp'
        [varargout{1:nargout}] = fwhm2nbcomp(varargin{:});
    case 'dcbasis'
        [varargout{1:nargout}] = dcbasis(varargin{:});
    case 'dfbasis'
        [varargout{1:nargout}] = dfbasis(varargin{:});
    case 'reconstruct'
        [varargout{1:nargout}] = reconstruct(varargin{:});
    case 'regulariser'
        [varargout{1:nargout}] = regulariser(varargin{:});
    case 'derivatives'
        [varargout{1:nargout}] = derivatives(varargin{:});
    case 'objective'
        [varargout{1:nargout}] = objective(varargin{:});
    case 'plot'
        [varargout{1:nargout}] = biasplot(varargin{:});
    otherwise
        error('Unknown function %s.', action);
end

% =========================================================================
function nbcmp = fwhm2nbcomp(lattice, vs, fwhm)
% FORMAT nbcmp = spm_bias_lib('fwhm2nbcomp', lattice, voxel_size, fwhm)
%
% lattice      - Dimensions of the lattice [dx dy ...]
% voxel_size   - Voxel size of the lattice [vx vy ...]
% fwhm         - Full-width half-max of the highest frequency basis (mm)
%
% The number of components is chosen so that the full-width half-max of
% the highest frequency basis function is smaller than fwhm. The bias
% field cannot model effects whose spatial frequency is higher than this
% value.
%
% If only one value is provided for voxel_size or fwhm, the same value is
% used along all dimensions.

% -------------------------------------------------------------------------
% Preprocess input arguments
ndim = numel(lattice);
vs = reshape(vs, 1, []);
if numel(vs) < ndim
    vs = spm_padarray(vs, [0 ndim-numel(vs)], 'replicate', 'post');
end
fwhm = reshape(fwhm, 1, []);
if numel(fwhm) < ndim
    fwhm = spm_padarray(fwhm, [0 ndim-numel(fwhm)], 'replicate', 'post');
end

% -------------------------------------------------------------------------
% Compute number of components per direction
nbcmp = ceil(2 * vs .* lattice ./ fwhm);
nbcmp = max(nbcmp, 1);


% =========================================================================
function varargout = dcbasis(lattice, nb_component)
% FORMAT [Bx,By,Bz,...] = spm_bias_lib('dcbasis', lattice, nb_component)
%
% lattice      - Dimensions of the lattice [dx dy ...]
% nb_component - Number of basis functions along each dimension [nx ny ...]
%
% Bx - Smooth basis along the x dimension [dx*nx]
% By - Smooth basis along the y dimension [dy*ny]
% ...
%
% There are as many basis objects as elements in `lattice`

ndim = numel(lattice);

% -------------------------------------------------------------------------
% Preprocess input arguments
nb_component = reshape(nb_component, 1, []);
if numel(nb_component) < ndim
    nb_component = spm_padarray(nb_component, [0 ndim-numel(nb_component)], 'replicate', 'post');
end

% -------------------------------------------------------------------------
% Compute each basis
varargout = cell(1,min(ndim, nargout));
for d=1:min(ndim, nargout)
    varargout{d} = spm_dctmtx(lattice(d),nb_component(d));
end

% =========================================================================
function varargout = dfbasis(lattice, nb_component)
% FORMAT [Bx,By,Bz,...] = spm_bias_lib('dfbasis', lattice, nb_component)
%
% lattice      - Dimensions of the lattice [dx dy ...]
% nb_component - Number of basis functions along each dimension [nx ny ...]
%
% Bx - Smooth basis along the x dimension [dx*nx]
% By - Smooth basis along the y dimension [dy*ny]
% ...
%
% There are as many basis objects as elements in `lattice`

ndim = numel(lattice);

% -------------------------------------------------------------------------
% Preprocess input arguments
nb_component = reshape(nb_component, 1, []);
if numel(nb_component) < ndim
    nb_component = spm_padarray(nb_component, [0 ndim-numel(nb_component)], 'replicate', 'post');
end

% -------------------------------------------------------------------------
% Compute each basis
varargout = cell(1,min(ndim, nargout));
for d=1:min(ndim, nargout)
    varargout{d} = spm_dftmtx(lattice(d),nb_component(d));
end

% =========================================================================
function L = regulariser(mode, lattice, nb_component, vs, bnd)
% FORMAT L = regulariser(param, lattice, nb_component, voxel_size)
% FORMAT L = regulariser(mode,  lattice, nb_component, voxel_size)
%
% param        - Parameters for absolute, membrane and bending energies
% mode         - Name of a single energy ('absolute'/'membrane'/'bending')
% lattice      - Dimensions of the lattice [dx dy ...]
% nb_component - Number of basis functions along each dimension [nx ny ...]
% voxel_size   - Voxel size of the lattice [vx vy ...]
%
% L            - Precision matrix [(nx*ny*...)^2]
%
% If numerical parameters are provided, a weighted combination of the
% three types of regularisation is returned.
% If an energy name is provided, the matrix that allows to compute it is
% returned (without weighting: the regularisation parameter should be
% multiplied with this matrix)
%
% If only one value is provided for nb_component or voxel_size, the
% same value is used along all dimensions.

if nargin < 5
    bnd = 'neumann';
end

% -------------------------------------------------------------------------
% Special case: mixture of regularisers
if ~ischar(mode)
    param = mode;
    L = 0;
    for i=1:numel(param)
        if param(i) ~= 0
            switch i
                case 1
                    mode = 'absolute';
                case 2
                    mode = 'membrane';
                case 3
                    mode = 'bending';
                case 4
                    mode = 'linearelastic1';
                case 5
                    mode = 'linearelastic2';
            end
            L1 = param(i) * regulariser(mode, lattice, nb_component, vs, bnd);
            if numel(L) == 1 || size(L,1) == size(L1,1)
                L = L + L1;
            else
                L0   = L;
                nprm = size(L,1);
                ndim = size(L1,1)/nprm;
                L = zeros(ndim*nprm);
                for d=1:ndim
                    L(((d-1)*nprm+1):d*nprm,((d-1)*nprm+1):d*nprm) = L0;
                end
                clear L0
                L = L + L1;
            end
        end
    end
    return
end


% -------------------------------------------------------------------------
% Preprocess input arguments
ndim = numel(lattice);
nb_component = reshape(nb_component, 1, []);
if numel(nb_component) < ndim
    nb_component = spm_padarray(nb_component, [0 ndim-numel(nb_component)], 'replicate', 'post');
end
if nargin < 4
    vs = 1;
end
vs = reshape(vs, 1, []);
if numel(vs) < ndim
    vs = spm_padarray(vs, [0 ndim-numel(vs)], 'replicate', 'post');
end

% -------------------------------------------------------------------------
% Mode-specific options
switch lower(mode)
    case {'absolute' 'abs' 'a'}
        maxdiff = 0;
    case {'membrane' 'mem' 'm' ...
          'linear-elastic1' 'linearelastic1' 'le1' ...
          'linear-elastic2' 'linearelastic2' 'le2'}
        maxdiff = 1;
    case {'bending' 'ben' 'b'}
        maxdiff = 2;
    otherwise
        error('Unknown mode %s, should be ''absolute'', ''membrane'' or ''bending''.', mode);
end

% -------------------------------------------------------------------------
% Compute each basis + square it
switch lower(bnd)
    case {0, 'circulant', 'circ', 'c'}
        mtxfun = @spm_dftmtx;
    case {1, 'neumann', 'neu', 'n'}
        mtxfun = @spm_dctmtx;
    case {2, 'dirichlet', 'dir', 'd'}
        mtxfun = @spm_dstmtx;
    otherwise
        error('Unknown boundary condition');
end

basis = cell(ndim, maxdiff + 1);
nbprm = 1;
for d=1:ndim
    for diff=0:maxdiff
        switch diff
            case 0
                basis{d,diff+1} = mtxfun(lattice(d),nb_component(d));
                nbprm = nbprm * size(basis{d,diff+1}, 2);
            case 1
                basis{d,diff+1} = mtxfun(lattice(d),nb_component(d),'diff') / vs(d);
            case 2
                basis{d,diff+1} = mtxfun(lattice(d),nb_component(d),'diff2') / vs(d)^2;
        end
        if any(strcmpi(mode, {'absolute' 'abs' 'a' 'membrane' 'mem' 'm' 'bending' 'ben' 'b'}))
            basis{d,diff+1} = basis{d,diff+1}.' * basis{d,diff+1};
        end
    end
end

% -------------------------------------------------------------------------
% Compute precision matrix
switch lower(mode)
    case {'absolute' 'abs' 'a'}
        L = 1;
        for d=1:ndim
            L = spm_krutil(basis{d,1}, L);
        end
    case {'membrane' 'mem' 'm'}
        L = 0;
        for dd=1:ndim               % Which dimension to differentiate
            L1 = 1;
            for d=1:ndim            % Kronecker loop
                if d == dd
                    L1 = spm_krutil(basis{d,2}, L1);
                else
                    L1 = spm_krutil(basis{d,1}, L1);
                end
            end
            L = L + L1;
        end
    case {'bending' 'ben' 'b'}
        L = 0;
        for dd1=1:ndim              % First dimension to differentiate
            L1 = 1;
            for d=1:ndim            % Kronecker loop
                if d == dd1
                    L1 = spm_krutil(basis{d,3}, L1);
                else
                    L1 = spm_krutil(basis{d,1}, L1);
                end
            end
            L = L + L1;
            for dd2=dd1+1:ndim      % Second dimension to differentiate
                L1 = 1;
                for d=1:ndim        % Kronecker loop
                    if d == dd1 || d == dd2
                        L1 = spm_krutil(basis{d,2}, L1);
                    else
                        L1 = spm_krutil(basis{d,1}, L1);
                    end
                end
                L = L + 2 * L1;
            end
        end
    case {'linear-elastic1' 'linearelastic1' 'le1'}
        L = zeros(nbprm,ndim,nbprm,ndim);
        for h1=1:ndim               % First Hessian dimension
            for dd=1:ndim          % First dimension to differentiate
                if dd == h1
                    coeff = 1;
                else
                    coeff = 0.5;
                end
                L1 = 1;
                for d=1:ndim        % Kronecker loop
                    if d == dd
                        L1 = spm_krutil(basis{d,2}.' * basis{d,2}, L1);
                    else
                        L1 = spm_krutil(basis{d,1}.' * basis{d,1}, L1);
                    end
                end
                L(:,h1,:,h1) = L(:,h1,:,h1) + coeff * reshape(L1, [nbprm 1 nbprm]);
            end

            for h2=h1+1:ndim        % Second Hessian dimension
                L1 = 1;
                for d=1:ndim        % Kronecker loop
                    if d == h1
                        L1 = spm_krutil(basis{d,1}.' * basis{d,2}, L1);
                    elseif d == h2
                        L1 = spm_krutil(basis{d,2}.' * basis{d,1}, L1);
                    else
                        L1 = spm_krutil(basis{d,1}.' * basis{d,1}, L1);
                    end
                end
                L(:,h1,:,h2) = L(:,h1,:,h2) + 0.5 * reshape(L1,  [nbprm 1 nbprm]);
                L(:,h2,:,h1) = L(:,h2,:,h1) + 0.5 * reshape(L1', [nbprm 1 nbprm]);
            end

        end
        L = reshape(L, nbprm*ndim, nbprm*ndim);
    case {'linear-elastic2' 'linearelastic2' 'le2'}
        L = zeros(nbprm,ndim,nbprm,ndim);
        for h1=1:ndim               % First Hessian dimension
            L1 = 1;
            for d=1:ndim        % Kronecker loop
                if d == h1
                    L1 = spm_krutil(basis{d,2}.' * basis{d,2}, L1);
                else
                    L1 = spm_krutil(basis{d,1}.' * basis{d,1}, L1);
                end
            end
            L(:,h1,:,h1) = L(:,h1,:,h1) + 0.5 * reshape(L1, [nbprm 1 nbprm]);

            for h2=h1+1:ndim        % Second Hessian dimension
                L1 = 1;
                for d=1:ndim        % Kronecker loop
                    if d == h1
                        L1 = spm_krutil(basis{d,2}.' * basis{d,1}, L1);
                    elseif d == h2
                        L1 = spm_krutil(basis{d,1}.' * basis{d,2}, L1);
                    else
                        L1 = spm_krutil(basis{d,1}.' * basis{d,1}, L1);
                    end
                end
                L(:,h1,:,h2) = L(:,h1,:,h2) + 0.5 * reshape(L1,  [nbprm 1 nbprm]);
                L(:,h2,:,h1) = L(:,h2,:,h1) + 0.5 * reshape(L1', [nbprm 1 nbprm]);
            end

        end
        L = reshape(L, nbprm*ndim, nbprm*ndim);
end

% =========================================================================
function field = reconstruct(basis, coeff, mode)
% FORMAT field = spm_bias_lib('reconstruct', {Bx,By,...}, coefficients)

if nargin < 3
    mode = 'mult';
end

% -------------------------------------------------------------------------
% Get number of components per basis
ndim = numel(basis);
ncomp = zeros(1,ndim);
lat   = zeros(1,ndim);
for i=1:ndim
    lat(i)   = size(basis{i}, 1);
    ncomp(i) = size(basis{i}, 2);
end
P = numel(coeff)/prod(ncomp); % Number of bias fields
ncomp(ndim+1) = P;

% -------------------------------------------------------------------------
% Coefficients provided
if ~isempty(coeff)
    field = coeff;
    for d=1:ndim
        field    = reshape(field, ncomp(1), []); % Coeffs in matrix form
        field    = basis{d} * field;             % Basis x Coeff
        ncomp(1) = size(field, 1);               % Update size (nbcoeffs -> nbvoxels)
        field    = reshape(field, ncomp);        % Coeffs in ND-array form
        ncomp    = circshift(ncomp, -1);         % Shift dimensions
        field    = shiftdim(field, 1);           % Shift dimensions
    end
    if P > 1
        ncomp    = circshift(ncomp, -1);         % Shift dimensions
        field    = shiftdim(field, 1);           % Shift dimensions
    end
    field = reshape(field, ncomp);               % Final reshape
    switch lower(mode)
        case 'add'
        case 'mult'
            field = exp(field);
        otherwise
            error('Unknown bias mode %s. Should be ''add'' or ''mult''.', mode)
    end

% -------------------------------------------------------------------------
% No coefficients provided
else
    switch lower(mode)
        case 'add'
            field  = zeros(lat, 'single');
        case 'mult'
            field  = ones(lat, 'single');
        otherwise
            error('Unknown bias mode %s. Should be ''add'' or ''mult''.', mode)
    end
end

% =========================================================================
function [g,H] = derivatives(p, X, B, Z, cluster, codes, binvar)
% FORMAT [g,H] = derivatives(p, obs, basis, resp, cluster, codes, binvar)
%
% p        -     Channel to process. If empty: compute joint gradient/Hessian
% obs      - NxP Bias corrected image
% basis    -     Smooth basis functions {Bx,By,...}
% resp     - NxK Cluster responsibilites
% cluster  -     Either {MU,A} or {MU,V,n} -> Gaussian mixture parameters
% codes    - Nx1 List of cdes encoding missing configuations: C or {C,L}
% binvar   - 1xP Binning uncertainty
%
% g - (J1xJ2x...)xP                Gradient w.r.t. bias coefficients
% H - (J1xJ2x...)xPx(J1xJ2x...)xP  Hessian  w.r.t. bias coefficients
%
% Compute gradient and Hessian of the conditional term


MU = [];
A  = [];
V  = [];
n  = [];
C  = [];
L  = [];
if nargin < 7
    binvar = 0;
end


%--------------------------------------------------------------------------
% Read input arguments
if ~iscell(cluster)
    MU = cluster;
else
    if numel(cluster) >= 1
        MU = cluster{1};
        if numel(cluster) >= 2
            A = cluster{2};
            if numel(cluster) >= 3
                V = A;
                n = cluster{3};
            end
        end
    end
end
if nargin >= 6
    if ~iscell(codes)
        C = codes;
    else
        if numel(codes) >= 1
            C = codes{1};
            if numel(codes) >= 2
                L = codes{2};
            end
        end
    end
    if isempty(L)
        L = unique(C);
    end
end

% -------------------------------------------------------------------------
% Dimensions
N  = size(X,1);
P  = size(MU,1);
K  = size(MU,2);
if isempty(L), L = 2^P - 1;    end % None missing
ndim  = numel(B);
ncomp = zeros(1,ndim);
lat   = zeros(1,ndim);
for i=1:ndim
    lat(i)   = size(B{i}, 1);
    ncomp(i) = size(B{i}, 2);
end

% -------------------------------------------------------------------------
% Initialise arrays to store statistics for gradient and Hessian
if ~isempty(p)
    g = zeros(N,1);    % <- 0.5 * Spp * x_p^2 - x_p * [Sp*(mu-0.5*x)]
    H = zeros(N,1);    % <- 1.5 * Spp * x_p^2 - x_p * [Sp*(mu-0.5*x)]
else
    g = zeros(N,P);    % [p]   <- 0.5 * Spp * x_p^2 - x_p * [Sp*(mu-0.5*x)]
    H = zeros(N,P,P);  % [p,p] <- 1.5 * Spp * x_p^2 - x_p * [Sp*(mu-0.5*x)]
                       % [p,q] <- 0.5 * Spq * x_p * x_q
end

% -------------------------------------------------------------------------
% For each combination of missing voxels
for i=1:numel(L)

    % ---------------------------------------------------------------------
    % Get mask of missing modalities (with this particular code)
    c        = L(i);
    observed = spm_gmm_lib('code2bin', c, P);
    missing  = ~observed;
    if ~isempty(p) && missing(p), continue; end
    if isempty(C), msk = ones(N, 1, 'logical');
    else,          msk = (C == c);
    end
    Pm      = sum(missing);
    Po      = P-Pm;
    Nc      = sum(msk);
    if Nc == 0, continue; end

    % ---------------------------------------------------------------------
    % Convert channel indices to observed indices
    if isempty(p)
        list_p = 1:Po;
    else
        list_p = 1:P;
        list_p = list_p(observed);
        [~,list_p] = find(list_p == p);
        if isempty(list_p)
            continue;
        end
    end

    X1 = X(msk,observed);

    for k=1:K

        Z1 = Z(msk,k);

        % -----------------------------------------------------------------
        % Compute expected precision (see GMM+missing data)
        if sum(n) > 0
            Ao = V(observed,observed,k) - V(observed,missing,k)*(V(missing,missing,k)\V(missing,observed,k));
            Ao = (n(k)-Pm) * Ao;
        else
            Ao = A(observed,observed,k) - A(observed,missing,k)*(A(missing,missing,k)\A(missing,observed,k));
        end
        MUo = MU(observed,k);

        % -----------------------------------------------------------------
        % Compute statistics
        sk1 = zeros(Nc,numel(list_p));
        sk2 = zeros(Nc,numel(list_p),numel(list_p));
        for iq=1:numel(list_p)
            q = list_p(iq);
            sk1(:,iq)   = X1(:,q) .* (bsxfun(@minus, X1, MUo.') * Ao(q,:).');
            sk2(:,iq,iq) = Ao(q,q) * X1(:,q).^2;
            if numel(binvar) > 1
                sk1(:,iq)    = sk1(:,iq)    + Ao(q,q) * binvar(msk,q);
                sk2(:,iq,iq) = sk2(:,iq,iq) + Ao(q,q) * binvar(msk,q);
            end
            if numel(list_p) > 1
                for qq=q+1:Po
                    sk2(:,q,qq) = Ao(q,qq) * X1(:,q) .* X1(:,qq);
                    sk2(:,qq,q) = sk2(:,q,qq);
                end
            end
        end
        sk1 = bsxfun(@times, sk1, Z1);
        sk2 = bsxfun(@times, sk2, Z1);

        % -----------------------------------------------------------------
        % Accumulate
        if isempty(p)
            g(msk,observed)          = g(msk,observed)          + sk1;
            H(msk,observed,observed) = H(msk,observed,observed) + sk2;
        else
            g(msk) = g(msk) + sk1;
            H(msk) = H(msk) + sk2;
        end
        clear sk1 sk2

    end
    % ---------------------------------------------------------------------
    % Normalisation term
    if isempty(p)
        g(msk,observed) = g(msk,observed) - 1;
    else
        g(msk) = g(msk) - 1;
    end

end

% -------------------------------------------------------------------------
% Multiply with basis functions
if isempty(p)
    % ---------------------------------------------------------------------
    % Gradient
    dimG = [lat P];
    g    = reshape(g, dimG);
    for d=1:ndim
        g       = reshape(g, dimG(1), []);  % Stats in matrix form
        g       = B{d}.' * g;               % Basis x Stats
        dimG(1) = size(g, 1);               % Update size (nbvoxels -> nbcoeffs)
        g       = reshape(g, dimG);         % Coeffs in ND-array form
        dimG    = circshift(dimG, -1);      % Shift dimensions
        g       = shiftdim(g, 1);           % Shift dimensions
    end
    if P > 1
        g       = shiftdim(g, 1);
    end
    g = reshape(g, [ncomp P]);
    % ---------------------------------------------------------------------
    % Hessian
    BB = ones(1, 'like', B{1});
    for d=1:ndim
        BB = spm_krutil(B{d}, BB);
    end
    BB = reshape(BB, [lat ncomp]);
    H  = bsxfun(@times, BB, reshape(H, [lat ones(1,numel(ncomp)) P P]));
    BB = reshape(BB, prod(lat), prod(ncomp));
    H  = BB' * reshape(H, prod(lat), prod(ncomp)*P*P); clear BB
    H  = reshape(H, [prod(ncomp) prod(ncomp) P P]);
    H  = reshape(permute(H, [1 3 2 4]), [ncomp P ncomp P]);
else
    % ---------------------------------------------------------------------
    % Gradient
    dimG = lat;
    g    = reshape(g, dimG);
    for d=1:ndim
        g       = reshape(g, dimG(1), []);  % Stats in matrix form
        g       = B{d}.' * g;               % Basis x Stats
        dimG(1) = size(g, 1);               % Update size (nbvoxels -> nbcoeffs)
        g       = reshape(g, dimG);         % Coeffs in ND-array form
        dimG    = circshift(dimG, -1);      % Shift dimensions
        g       = shiftdim(g, 1);           % Shift dimensions
    end
    g = reshape(g, ncomp);
    % ---------------------------------------------------------------------
    % Hessian
    BB = ones(1, 'like', B{1});
    for d=1:ndim
        BB = spm_krutil(B{d}, BB);
    end
    BB = reshape(BB, [lat ncomp]);
    H  = bsxfun(@times, BB, reshape(H, lat));
    BB = reshape(BB, [], prod(ncomp));
    H  = BB' * reshape(H, [], prod(ncomp)); clear BB
    H  = reshape(H, [ncomp ncomp]);
end


% =========================================================================
function lb = objective(X, Z, B, mean, prec, codes, binwidth)
% FORMAT lb = spm_bias_lib('objective', obs, resp, bias, mean, prec, codes, binvar)
%
% MANDATORY
% ---------
% obs    - NxP      observations (non-corrected)
% resp   - NxK      responsibilities
% bias   - NxP      bias field (exponentiated)
% mean   - PxK      GMM mean:      MU or {MU,b}
% prec   - PxPxK    GMM precision: A  or {V,n}
%
% OPTIONAL
% --------
% codes    - Nx1   image of missing codes (and code list): C or {C,L}
% binwidth - 1xP   bin width
%
% Compute the conditional data term of the objective function:
% sum_n { sum_k log zk * N( Bx | MUk, Ak ) }

if nargin < 7, binwidth = 0; end

C  = [];
L  = [];

%--------------------------------------------------------------------------
% Read input arguments
if nargin >= 6
    if ~iscell(codes)
        C = codes;
    else
        if numel(codes) >= 1
            C = codes{1};
            if numel(codes) >= 2
                L = codes{2};
            end
        end
    end
    if isempty(L)
        L = unique(C);
    end
end

%--------------------------------------------------------------------------
% Normalisation term
if sum(binwidth) > 0
    binwidth = reshape(binwidth, 1, []);
    binvar   = (bsxfun(@times, binwidth, B).^2)/12;
end

%--------------------------------------------------------------------------
% Compute GMM likelihood from bias-corrected data
X = X .* B;
[lSS0,lSS1,lSS2] = spm_gmm_lib('SuffStat', 'base', X, Z, 1, {C,L});
SS2b = 0;
if sum(binwidth) > 0
    SS2b = spm_gmm_lib('SuffStat', 'bin', binvar, Z, 1, {C,L});
end
lb = spm_gmm_lib('MarginalSum', lSS0, lSS1, lSS2, mean, prec, L, SS2b);

% =========================================================================
function varargout = biasplot(action,varargin)
% Custom visualisation tools for Gaussian Mixture modelling
%
% spm_bias_lib('plot', 'lb', lb, [figname])
% > Plot lower bound
% spm_bias_lib('plot', 'bias', X, B, lat, [figname])
% > Plot bias field

switch lower(action)
    case {'lowerbound','lb'}
        [varargout{1:nargout}] = plot_lowerbound(varargin{:});
    case {'bias'}
        [varargout{1:nargout}] = plot_bias(varargin{:});
    otherwise
        error('Unknown function %s.', action);
end

% =========================================================================
function plot_lowerbound(lb, figname)

% -------------------------------------------------------------------------
% Get figure (create if it does not exist)
if nargin < 2
    figname = '(SPM) Plot Bias Lower Bound';
end
f = findobj('Type', 'Figure', 'Name', figname);
if isempty(f)
    f = figure('Name', figname, 'NumberTitle', 'off');
end
set(0, 'CurrentFigure', f);
clf(f);

% -------------------------------------------------------------------------
% Plots
subplot(1, 3, 1);
plot(lb.sum)
title('Lower Bound')
subplot(1, 3, 2);
plot(sum(lb.X,1) + sum(lb.XB,1));
title('Conditional')
subplot(1, 3, 3);
plot(lb.B)
title('Regularisation)')
drawnow


% =========================================================================
function plot_bias(X, B, lat, figname)

% -------------------------------------------------------------------------
% Get figure (create if it does not exist)
if nargin < 4
    figname = '(SPM) Plot bias field';
end
f = findobj('Type', 'Figure', 'Name', figname);
if isempty(f)
    f = figure('Name', figname, 'NumberTitle', 'off');
end
set(0, 'CurrentFigure', f);
clf(f);

% -------------------------------------------------------------------------
% Lattice
if nargin < 4
    lat = [size(X) 1];
    lat = lat(1:3);
end
clat = num2cell(lat);
X = reshape(X, clat{:}, []);
B = reshape(B, clat{:}, []);
P = size(X,4);
z = ceil(size(X,3)/2);

% -------------------------------------------------------------------------
% Choose type
nrow = P;
ncol = 3;

% -------------------------------------------------------------------------
% Plots
handles = cell(nrow,ncol);
for p=1:P
    subplot(nrow, ncol, sub2ind([ncol nrow], 1, p));
    tmp = X(:,:,z,p);
    minval = min(tmp(:));
    maxval = max(tmp(:));
    handles{p,1} = imagesc(tmp(end:-1:1,end:-1:1)');
    colormap(handles{p,1}.Parent, 'gray')
    colorbar
    axis off
    box on
    title(sprintf('Original %d', p));

    subplot(nrow, ncol, sub2ind([ncol nrow], 2, p));
    tmp = tmp .* B(:,:,z,p);
    minval = min(minval, min(tmp(:)));
    maxval = max(maxval, max(tmp(:)));
    handles{p,2} = imagesc(tmp(end:-1:1,end:-1:1)');
    colormap(handles{p,2}.Parent, 'gray')
    colorbar
    axis off
    box on
    title(sprintf('Corrected %d', p));

    caxis(handles{p,1}.Parent, [minval maxval]);
    caxis(handles{p,2}.Parent, [minval maxval]);

    subplot(nrow, ncol, sub2ind([ncol nrow], 3, p));
    tmp = B(:,:,z,p);
    handles{p,3} = imagesc(tmp(end:-1:1,end:-1:1)');
    colorbar
    axis off
    box on
    title(sprintf('Bias field %d', p));
end
drawnow

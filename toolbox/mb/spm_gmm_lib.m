function varargout = spm_gmm_lib(varargin)
%__________________________________________________________________________
%
% Library of functions for Gaussian Mixture modelling
%
%--------------------------------------------------------------------------
% Main function
% -------------
% * Loop                >> Core fit (assumes well formatted input)
%--------------------------------------------------------------------------
% Main helpers
% ------------
% * InferMissing        >> Infer missing values
% * Marginal            >> Compute log marginal distribution: E[log p(X|...)]
% * Responsibility      >> Compute posterior responsibilities: q(Z|...)
% * SuffStat            >> Compute sufficient statistics 0th/1st/2nd order
% * Normalisation       >> Compute normalisation term of the log marginal
%--------------------------------------------------------------------------
% Update parameters
% -----------------
% * UpdateClusters      >> Update posterior cluster parameters: q(MU,A)
% * UpdateProportions   >> Update posterior proportion parameters: q(PI)
% * updateHyperPars     >> Update prior cluster parameters: MU0,b0,V0,n0
%--------------------------------------------------------------------------
% Lower bound
%------------
% * MarginalSum         >> Compute the sum of the log marginal from suff stat
% * KL                  >> Useful KL divergences
%--------------------------------------------------------------------------
% "Missing code" image
%---------------------
% * obs2cell            >> Convert obs matrix to missing-friendly cell
% * cell2obs            >> Convert missing-friendly cell to obs matrix
%--------------------------------------------------------------------------
% Plot
%-----
% * Plot                >> Plotting utilitites
%--------------------------------------------------------------------------
% Extras
%-------
% * Extras              >> Extra utilitites
%__________________________________________________________________________

%--------------------------------------------------------------------------
% Convention
% ----------
% N  - Number of observations
% Nm - Number of observations that correspond to a given "missing code"
% P  - Dimension of observations
% Po - Number of non-missing observations in a given "missing code"
% Pm - Number of     missing observations in a given "missing code"
% K  - Number of clusters
% M  - Number of unique "missing codes"
%
% X   -  N  x P        Observations                    (+ inferred missing)
% Z   -  N  x K        Clusters' responsibility
% PI  - [N] x K        Clusters' proportions
% a   -  1  x K        Proportions concentration       (Dirichlet prior)
% W   - [N] x 1        Observations weights            (Histogram GMM)
% MU  -  P  x K        Clusters' [expected] mean
% A   -  P  x P  x K   Clusters' [expected] precision matrix
% b   -  1  x K        Mean degrees of freedom         (Gauss prior)
% V   -  P  x P  x K   Precision scale matrix          (Wishart prior)
% n   -  1  x K        Precision degrees of freedom    (Wishart prior)
% U   - [N] x P        Uncertainty about observations  (Histogram GMM)
% SS0 -  1  x K        Zeroth order sufficient statistics
% SS1 -  P  x K        First  order sufficient statistics
% SS2 -  P  x P  x K   Second order sufficient statistics
% C   -  N  x 1        Code image: one code / missing combination
% mask-  M  x P        Mask of observed channels per 'missing code'
%--------------------------------------------------------------------------

% $Id$

if nargin == 0
    help spm_gmm_lib
    error('Not enough argument. Type ''help spm_gmm_lib'' for help.');
end
id = varargin{1};
varargin = varargin(2:end);
switch lower(id)
    case 'loop'
        [varargout{1:nargout}] = loop(varargin{:});
%--------------------------------------------------------------------------
% Main helpers
%--------------------------------------------------------------------------
    case 'infermissing'
        [varargout{1:nargout}] = infermissing(varargin{:});
        % X = spm_gmm_lib('InferMissing', X, Z, {MU,A}, {code_im,code_list})
        % > Infer missing values (this should only be done once, at the end)
    case 'marginal'
        [varargout{1:nargout}] = marginal(varargin{:});
        % logp = spm_gmm_lib('Marginal', X, {MU,A},   const, mask, U)
        % logp = spm_gmm_lib('Marginal', X, {MU,V,n}, const, mask, U)
        % > Observation's marginal log probability within each cluster
    case 'responsibility'
        [varargout{1:nargout}] = responsibility(varargin{:});
        % Z = spm_gmm_lib('Responsibility', logpX, logPI)
        % > Compute & normalise responsibilities (safe softmax)
    case 'suffstat'
        [varargout{1:nargout}] = suffstat(varargin{:});
        % FORMAT [SS0,SS1,SS2] = spm_gmm_lib('SuffStat', X, Z, W, (mask))
        % FORMAT [SS0,SS1,SS2] = spm_gmm_lib('SuffStat', 'infer', SS0, SS1, SS2, {MU,A}, mask)
        % FORMAT         [SS2] = spm_gmm_lib('SuffStat', 'uncertainty', U, Z, W, (mask))
        % > Compute sufficient statistics (0th, 1st, 2nd order)
        %   default: no missing data            => E[z], E[z]*x, E[z]*xx'
        %            with mising data           => E[z], E[z]*g, E[z]*gg'
        %   'infer': base to full statistics    => E[z], E[z*x], E[z*xx']
        %   'uncertainty': observed uncertainty => Tr(S\cov[gg'])
    case {'normalisation' 'normalization'}
        [varargout{1:nargout}] = normalisation(varargin{:});
        % const = spm_gmm_lib('Normalisation',  MU,     A,    (mask))
        % const = spm_gmm_lib('Normalisation', {MU,b}, {V,n}, (mask))
        % > Normalisation term of a Gaussian log-distribution
        %   If mask is provided -> marginal distributions
%--------------------------------------------------------------------------
% Update parameters
%--------------------------------------------------------------------------
    case 'updateclusters'
        [varargout{1:nargout}] = updateclusters(varargin{:});
        % [MU,A]       = spm_gmm_lib('UpdateClusters', SS0, SS1, SS2)
        % [MU,A,b,V,n] = spm_gmm_lib('UpdateClusters', SS0, SS1, SS2, {MU0,b0,V0,n0})
        % > Update GMM parameters (ML or Bayesian posterior)
    case 'updateproportions'
        [varargout{1:nargout}] = updateproportions(varargin{:});
        % [PI,logPI,a] = spm_gmm_lib('UpdateProportions', SS0, a0)
        % > Update cluster proportions
    case 'updatehyperpars'
        [varargout{1:nargout}] = updatehyperpars(varargin{:});
        % [GaussPrior,extras] = spm_gmm_lib('updatehyperpars',cluster,GaussPrior,varargin)
        % > Update VB-GMM hyper-parameters (MU,b,V,n)
%--------------------------------------------------------------------------
% Lower bound
%--------------------------------------------------------------------------
    case 'marginalsum'
        [varargout{1:nargout}] = marginalsum(varargin{:});
        % [lb,const] = spm_gmm_lib('MarginalSum', SS0, SS1, SS2,  MU,     A,    mask, SS2u)
        % [lb,const] = spm_gmm_lib('MarginalSum', SS0, SS1, SS2, {MU,b}, {V,n}, mask, SS2u)
        % > Compute conditional datasum: E[ln p(g|MU,A,Z)]
        %   Also returns the result of spm_gmm_lib('const')
    case 'kl'
        [varargout{1:nargout}] = kl(varargin{:});
        % [klMU,klA] = spm_gmm_lib('KL', 'GaussWishart', {MU,b}, {V,n}, {MU0,b0}, {V0,n0})
        % > KL divergence between two Gauss-Wishart distributions
        %
        % klP = spm_gmm_lib('KL', 'Dirichlet', a, a0)
        % > KL divergence between two Dirichlet distributions
        %
        % klZ = spm_gmm_lib('KL', 'Categorical', Z, W, logPI)
        % > KL divergence between two Categorical distributions
%--------------------------------------------------------------------------
% "Missing code" image
%--------------------------------------------------------------------------
    case 'obs2cell'
        [varargout{1:nargout}] = obs2cell(varargin{:});
        % [X,code_im,mask] = spm_gmm_lib('obs2cell', X)
        % > Transform a matrix of observations into a cell of matrices for
        %   each missing pattern.
    case 'cell2obs'
        [varargout{1:nargout}] = cell2obs(varargin{:});
        % X = spm_gmm_lib('cell2obs', X, code_im, mask)
        % > Create a matrix of observations from a cell of matrices.
%--------------------------------------------------------------------------
% Visualisation
%--------------------------------------------------------------------------
    case 'plot'
        [varargout{1:nargout}] = gmmplot(varargin{:});
        % spm_gmm_lib('Plot', 'LB', lb)
        % > Plot lower bound
        %
        % spm_gmm_lib('Plot', 'GMM', X, W, mask, {MU,A}, PI)
        % > Plot mixture fit
        %
        % spm_gmm_lib('plot', 'cat', dm, Z, Template, (wintitle))
        % > Plot (categorical) responsibilities and template (if available)
        %
        % spm_gmm_lib('plot', 'gaussprior', GaussPrior, (wintitle))
        % > Plot VB-GMM hyper-parameters
        %
        % c = spm_gmm_lib('plot', 'cat2rgb', f, pal)
        % > Generate an RGB volume from a categorical (e.g. responsibilities) volume.
%--------------------------------------------------------------------------
% Extras
%--------------------------------------------------------------------------
    case 'extras'
        [varargout{1:nargout}] = gmm_extras(varargin{:});
        % gmm = spm_gmm_lib('extras', 'more_gmms', gmm, part)
        % > A crude heuristic to replace a single Gaussian by a bunch of Gaussians.
    otherwise
        help spm_gmm_lib
        error('Unknown function %s. Type ''help spm_gmm_lib'' for help.', id)
end

%--------------------------------------------------------------------------
% Main loop
%--------------------------------------------------------------------------

function [Z,cluster,prop,lb,mg_w] = loop(X, weights, cluster, props, varargin)
%__________________________________________________________________________
%
% Fit a [Bayesian] Gaussian mixture model to observed [weighted] data.
%
% This function is the core of the fitting process. However, it needs all
% inputs to be well formatted and initialised and is, thus, not the usual
% entry point. To fit a GMM without having to bother with these issues, use
% spm_gmm instead.
%
% FORMAT [resp,cluster,prop,lb,mg_w,ss] = spm_gmm_lib('loop',obs,weights,cluster,prop,...)
%
% MANDATORY
% ---------
%
% obs <- X
%   X - {NoxPo} observations
%
% weights <- W
%   W - {Nox1}  weights [1]
%
% cluster <- {MU,A}, {{MU,b},A}, {MU,{V,n}}, or {{MU,b},{V,n}}
%   MU - PxK   means
%   b  - 1xK   mean d.f. [0=ML]
%   A  - PxPxK precision matrices
%   V  - PxPxK scale matrices
%   n  - 1xK   precision d.f. [0=ML]
%
% prop <- LogPi or {('LogProp', LogPi), ('Prop', Pi), ('Dir', a)}
%   LogPi - {NoxK} Fixed voxel-wise log-proportions
%           1xK    Pre-computed log(Pi) or E[log(Pi)]
%   Pi    - {NoxK} Fixed voxel-wise proportions
%           1xK    Pre-computed Pi or E[Pi]
%   a     - 1xK    Posterior concentration parameter (Dirichlet)
%
% KEYWORD
% -------
%
% LowerBound     - Pre-computed lower bound structure with fields:
%                   sum, last, X, Z, P, MU, A
% GaussPrior     - {MU0,b0,V0,n0} [{}=ML]
% PropPrior      - a0 [0=ML]
% Missing        - Infer missing data [true]
% Missing        - MxP Mask of observed channels per code
% IterMax        - Max number of EM iterations [1024]
% Tolerance      - Gain tolerance to stop the EM algorithm [1e-4]
% SubIterMax     - Max number of sub-EM iterations (Missing == true) [1024]
% SubTolerance   - Sub-EM gain tolerance (Missing == true) [1e-4]
% ObsUncertainty - 1xP     Uncertainty (= variance) about the observations
%                  {NoxPo} Voxel-wise uncertainty
% Verbose        - Verbosity level: [0]= quiet
%                                    1 = write (lower bound)
%                                    2 = plot (lower bound)
%                                    3 = plot more (gmm fit)
% Labels         - {NoxK} Log of voxel-wise labels (from, e.g., manual
%                  segmentations) [[]]
% MultGaussPi    - {[1xKmg],[1xKmg]} For using multiple Gaussians per class in
%                  proportions (Pi). Defined by a cell array with: the first
%                  element being a vector of length Kmg that maps indices of
%                  Gaussians to classes in Pi; the second element being a
%                  vector of the same length with mixing proportions [{}]
%
% OUTPUT
% ------
% resp       - Responsibilities
% cluster    - Structure with fields: MU, b, A, V, n
% prop       - Structure with fields: LogProp, Prop, Dir
% lb         - Structure with fields: sum, last, X, Z, P, MU, A
% mg_w       - Vector with weights
%
%__________________________________________________________________________
% Copyright (C) 2018 Wellcome Centre for Human Neuroimaging

lb0 = struct('sum', NaN, ...
             'X', [], 'Z', [], 'P', [], 'MU', [], 'A', []);

% -------------------------------------------------------------------------
% Parse inputs
p = inputParser;
p.FunctionName = 'spm_gmm_loop';
p.addParameter('LowerBound',     lb0,   @isstruct);
p.addParameter('Resp',           [],    @(X) isnumeric(X) || iscell(X));
p.addParameter('GaussPrior',     {},    @iscell);
p.addParameter('PropPrior',      0,     @isnumeric);
p.addParameter('Missing',        {},    @(X) islogical(X) || iscell(X));
p.addParameter('IterMax',        1024,  @(X) isscalar(X) && isnumeric(X));
p.addParameter('Tolerance',      1e-4,  @(X) isscalar(X) && isnumeric(X));
p.addParameter('SubIterMax',     1024,  @(X) isscalar(X) && isnumeric(X));
p.addParameter('SubTolerance',   1e-4,  @(X) isscalar(X) && isnumeric(X));
p.addParameter('ObsUncertainty', 0,     @(X) isnumeric(X) || iscell(X));
p.addParameter('Verbose',        0,     @(X) isnumeric(X) || islogical(X));
p.addParameter('Labels',         [],    @(X) isnumeric(X) || iscell(X));
p.addParameter('MultGaussPi',    {},    @iscell);
p.parse(varargin{:});
lb              = p.Results.LowerBound;
prop_prior      = p.Results.PropPrior;
obs_channels    = p.Results.Missing;
obs_uncertainty = p.Results.ObsUncertainty;
gauss_prior     = p.Results.GaussPrior;
iter_max        = p.Results.IterMax;
tolerance       = p.Results.Tolerance;
subiter_max     = p.Results.SubIterMax;
subtolerance    = p.Results.SubTolerance;
verbose         = p.Results.Verbose;
labels          = p.Results.Labels;
mult_gauss      = p.Results.MultGaussPi;

% -------------------------------------------------------------------------
% Unfold inputs
b     = 0;  % Mean degrees of freedom (posterior)
n     = 0;  % Precision degrees of freedom (posterior)
V     = []; % Scale matrix (posterior)
if ~iscell(cluster) || numel(cluster) < 2
    error('At least one mean and one precision matrix are needed.');
else
    if ~iscell(cluster{1})
        MU = cluster{1};
    else
        MU = cluster{1}{1};
        if numel(cluster{1}) >= 2
            b = cluster{1}{2};
        end
    end
    if ~iscell(cluster{2})
        A = cluster{2};
    else
        A = cluster{2}{1};
        if numel(cluster{2}) >= 2
            n = cluster{2}{2};
            if sum(n) > 0
                V = A;
                A = bsxfun(@times, V, reshape(n, 1, 1, []));
            end
        end
    end
end
if sum(b) > 0
    mean = {MU,b};
else
    mean = MU;
end
if sum(n) > 0
    prec = {V,n};
else
    prec = A;
end

% --- Gauss-Wishart prior
b0    = 0;  % Mean degrees of freedom (prior)
n0    = 0;  % Precision degrees of freedom (prior)
MU0   = []; % Mean (prior)
V0    = []; % Scale matrix (prior)
if numel(gauss_prior) >= 1
    MU0 = gauss_prior{1};
    if numel(gauss_prior) >= 2
        b0 = gauss_prior{2};
        if numel(gauss_prior) >= 3
            V0 = gauss_prior{3};
            if numel(gauss_prior) >= 4
                n0 = gauss_prior{4};
            end
        end
    end
end
mean0 = {MU0 b0};
prec0 = {V0 n0};

% --- Proportions
log_prop       = []; % [expected] log-proportions
prop           = []; % [expected] proportions
prop_posterior = []; % Dirichlet posterior
if ~iscell(props)
    log_prop = props;
else
    i = 1;
    while i < numel(props)
        switch lower(props{i})
            case 'logprop'
                i = i + 1;
                log_prop = props{i};
            case 'prop'
                i = i + 1;
                prop = props{i};
            case 'dir'
                i = i + 1;
                prop_posterior = props{i};
            otherwise
                log_prop = props{i};
        end
        i = i + 1;
    end
end

% -------------------------------------------------------------------------
% For multiple Gaussians per class in Pi
if isempty(mult_gauss)
    mg_ix = 1:size(MU,2);
    mg_w  = ones([1 size(MU,2)]);
else
    mg_ix = mult_gauss{1};
    mg_w  = mult_gauss{2};
end

% -------------------------------------------------------------------------
% Compute log-prop if needed
if isempty(log_prop)
    if sum(prop_posterior) > 0
        if isempty(prop)
            prop = prop_posterior ./ sum(prop_posterior);
        end
        log_prop = psi(prop_posterior) - psi(sum(prop_posterior));
    elseif ~isempty(prop)
        log_prop = log(bsxfun(@rdivide, prop+eps, sum(prop+eps, 2)));
    else
        error('At least one of Prop, LogProp or Dir must be provided.');
    end
end

% -------------------------------------------------------------------------
% logpX (needed to initialise Z)
norm_term = normalisation(mean, prec, obs_channels);
logpX     = marginal(X, [{MU} prec], norm_term, obs_channels, obs_uncertainty);

% -------------------------------------------------------------------------
% EM loop
for em=1:iter_max

    % ---------------------------------------------------------------------
    % Compute responsibilities
    Z = responsibility(logpX, log_prop, labels, log(mg_w));
    clear logpX

    % ---------------------------------------------------------------------
    % Compute sufficient statistics (bin uncertainty part)
    if iscell(obs_uncertainty) ||  sum(obs_uncertainty) > 0
    	SS2u = suffstat_uncertainty(obs_uncertainty, Z, weights, obs_channels);
    else
        SS2u = 0;
    end

    if ~isempty(obs_channels)
    % ---------------------------------------------------------------------
    % sub-EM algorithm to update Mean/Precision with missing data
    % . Responsibilities (E[z]) are kept fixed
    % . Missing values (E[z*h], E[z*hh']) are updated
    % . Cluster parameters (MU,b,A/V,n) are updated

        % -----------------------------------------------------------------
        % Compute fast sufficient statistics:
        % > sum{E[z]}, sum{E[z]*g}, sum{E[z]*gg'}
        %   for each configuration of missing data
        [SS0m,SS1m,SS2m] = suffstat_missing(X, Z, weights, obs_channels);

        % -----------------------------------------------------------------
        % Initialise objective function
        [LMU,LA]       = kl_gausswishart(mean, prec, mean0, prec0);
        [LX,norm_term] = marginalsum(SS0m, SS1m, SS2m, mean, prec, obs_channels, SS2u);
        LB             = NaN(1,subiter_max);
        LB(1)          = LMU + LA + LX;
        for i=1:subiter_max

            % -------------------------------------------------------------
            % Save previous value
            Ap      = A;
            Vp      = V;
            np      = n;

            if numel(SS0m)==1, subsubiter_max = 1; else subsubiter_max = 4; end

            for ii=1:subsubiter_max
                % -------------------------------------------------------------
                % Infer missing suffstat
                % sum{E[z]}, sum{E[z*x]}, sum{E[z*xx']}
                [SS0,SS1,SS2] = suffstat_infer(SS0m, SS1m, SS2m, {MU,A}, obs_channels);
                SS2 = SS2 + SS2u;

                % -------------------------------------------------------------
                % Update GMM
                [MU,A,b,V,n] = updateclusters(SS0, SS1, SS2, [mean0 prec0]);
                for k=1:size(MU,2)
                    [~,cholp] = chol(A(:,:,k));
                    if cholp ~= 0
                        warning('A not positive definite - reverting to previous version')
                        A(:,:,k) = Ap(:,:,k);
                        if sum(n) > 0
                            V(:,:,k) = Vp(:,:,k);
                            n(k)     = np(k);
                        end
                    end
                end
                mean = {MU,b};
                if ~sum(n), prec = {A};
                else        prec = {V,n};   end
            end

            % -------------------------------------------------------------
            % Marginal / Objective function
            [LMU,LA]       = kl_gausswishart(mean, prec, mean0, prec0);
            [LX,norm_term] = marginalsum(SS0m, SS1m, SS2m, mean, prec, obs_channels, SS2u);
            LB(i+1)        = LMU+LA+LX;
            subgain        = (LB(i+1)-LB(i))/(max(LB(2:i+1))-min(LB(2:i+1)));
            subgain1       = (LB(i+1)-LB(i))/abs(LB(i+1));

            % -------------------------------------------------------------
            % Print stuff
            if numel(verbose) > 1 && verbose(2) > 0
                switch sign(subgain)
                    case 1,     incr = '(+)';
                    case -1,    incr = '(-)';
                    case 0,     incr = '(=)';
                    otherwise,  incr = '';
                end
                fprintf('%-5s | %4d | lb = %-12.6g | gain = %-10.4g | %3s\n', 'sub', i, LB(i+1), subgain, incr);
            end
            if numel(SS0m)==1 || subgain < subtolerance || subgain1 < eps('single')
                break
            end
        end

    else
        % ---------------------------------------------------------------------
        % Classical M-step

        % -----------------------------------------------------------------
        % Compute sufficient statistics
        [SS0,SS1,SS2] = suffstat_classic(X, Z, weights);
        SS2 = SS2 + SS2u;

        % -------------------------------------------------------------
        % Update GMM
        [MU,A,b,V,n] = updateclusters(SS0, SS1, SS2, [mean0 prec0]);
        for k=1:size(MU,2)
            [~,cholp] = chol(A(:,:,k));
            if cholp ~= 0
                warning('A not positive definite');
            end
        end
        mean = {MU,b};
        if ~sum(n), prec =  {A};
        else        prec = {V,n};   end
        norm_term = normalisation(mean, prec);
    end

    % ---------------------------------------------------------------------
    % Update Proportions
    if size(prop,1) == 1
        [prop,log_prop,prop_posterior] = updateproportions(SS0, prop_prior);
    end

    % ---------------------------------------------------------------------
    % Update weight for multiple Gaussians per prop class
    for k=1:size(MU,2)
        tmp     = SS0(mg_ix == mg_ix(k));
        mg_w(k) = (SS0(k) + eps*eps)/sum(tmp + eps*eps);
    end

    % ---------------------------------------------------------------------
    % Plot GMM
    if verbose(1) >= 3
        plot_gmm(X, weights, obs_channels, {MU,A}, prop);
    end

    % ---------------------------------------------------------------------
    % Marginal / Objective function
    logpX = marginal(X, [{MU} prec], norm_term, obs_channels, obs_uncertainty);

    % ---------------------------------------------------------------------
    % Compute lower bound
    lb.P(end+1) = kl_dirichlet(prop_posterior, prop_prior);
    lb.Z(end+1) = kl_categorical(Z, weights, log_prop, labels, log(mg_w));
    if isempty(obs_channels)
        [lb.MU(end+1),lb.A(end+1)] = kl_gausswishart({MU,b}, prec, {MU0,b0}, {V0,n0});
        lb.X(end+1) = sum(sum(bsxfun(@times, logpX, bsxfun(@times, Z, weights)),2),'double');
    else
        lb.MU(end+1) = LMU;
        lb.A(end+1)  = LA;
        lb.X(end+1)  = LX;
    end
    if isfield(lb, 'XB')
        if isempty(lb.XB)
            lb.XB = 0;
        else
            lb.XB(:,end+1) = lb.XB(:,end); % < Add bias normalisation term
        end
    end

    % ---------------------------------------------------------------------
    % Check convergence
    [lb,gain] = check_convergence(lb, em, verbose(1));
    if gain < tolerance
        break;
    end

end

% ---------------------------------------------------------------------
% Compute final responsibilities
Z = responsibility(logpX, log_prop, labels, log(mg_w));
clear logpX

% ---------------------------------------------------------------------
% Recompute parts of lower bound that depends on responsibilities
lb.Z(end+1) = kl_categorical(Z, weights, log_prop, labels, log(mg_w));
if isempty(obs_channels)
    lb.X(end+1) = sum(sum(bsxfun(@times, logpX, bsxfun(@times, Z, weights)),2),'double');
else
    [SS0m,SS1m,SS2m] = suffstat_missing(X, Z, weights, obs_channels);
    LX               = marginalsum(SS0m, SS1m, SS2m, mean, prec, obs_channels, SS2u);
    lb.X(end+1)      = LX;
end
lb = check_convergence(lb, em, verbose(1)); % Make sure lb.sum is correct

% -------------------------------------------------------------------------
% Format output
cluster    = struct('MU', MU, 'b', b, 'A', A, 'V', V, 'n', n);
prop       = struct('LogProp', log_prop, 'Prop', prop, 'Dir', prop_posterior);

% =========================================================================
function [lb,gain] = check_convergence(lb, em, verbose)
% FORMAT [lb,gain] = check_convergence(lb, em, verbose)
% lb      - Lower bound structure with fields X, Z, P, MU, A, sum, last
% em      - EM iteration
% verbose - Verbosity level (>= 0)
%
% Compute lower bound (by summing its parts) and its gain
% + print info

fields = fieldnames(lb);
lb.sum(end+1) = 0;
for i=1:numel(fields)
    field = fields{i};
    if ~any(strcmpi(field, {'sum' 'last'})) && ~isempty(lb.(field)) && ~isnan(lb.(field)(end))
        lb.sum(end) = lb.sum(end) + sum(lb.(field)(:,end));
    end
end
gain = (lb.sum(end) - lb.sum(end-1))/(max(lb.sum(:))-min(lb.sum(:)));
if verbose >= 1
    if verbose >= 2
        plot_lowerbound(lb);
    end
    switch sign(gain)
        case  1,    incr = '(+)';
        case -1,    incr = '(-)';
        case  0,    incr = '(=)';
        otherwise,  incr = '';
    end
    fprintf('%-5s | %4d | lb = %-12.6g | gain = %-10.4g | %3s\n', 'gmm', em, lb.sum(end), gain, incr);
end
gain = abs(gain);

%--------------------------------------------------------------------------
% Update functions
%--------------------------------------------------------------------------

% =========================================================================
function X = infermissing(X, Z, cluster, codes, sample)
% FORMAT X = spm_gmm_lib('missing', X, Z, {MU,A}, {C,L})
% X  - NxP   observations
% Z  - NxK   responsibilities
% MU - PxK   (expected) means
% A  - PxPxK (expected) precision matrices
% C  - Nx1   "missing value" code image
% L  -       list of existing codes
%
% X - NxP    observations with inferred values
%
% Compute the mean expected value of missing voxels.

if nargin < 5
    sample = false;
end

MU = [];
A  = [];
C  = [];
L  = [];

%--------------------------------------------------------------------------
% Read input arguments
if ~iscell(cluster)
    MU = cluster;
else
    if numel(cluster) >= 1
        MU = cluster{1};
        if numel(cluster) >= 2
            A = cluster{2};
        end
    end
end
if nargin >= 4
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
% Dimensions
P = size(X, 2);
K = size(Z, 2);
if isempty(L)
    L = 2^P - 1; % None missing
end

% -------------------------------------------------------------------------
% For each missing combination
for i=1:numel(L)

    % ---------------------------------------------------------------------
    % Get code / mask / missing modalities
    c   = L(i);
    io  = code2bin(c, P);
    im  = ~io;
    Pm  = sum(im);
    if Pm == 0, continue; end
    msk = (C == c);
    Nm  = sum(msk);
    if Nm == 0, continue; end

    % ---------------------------------------------------------------------
    % Initialise
    X(msk,im) = 0;

    % ---------------------------------------------------------------------
    % Compute posterior mean (expected value)
    % 1) t = sum_k {z * ( mu[m] + A[m]/A[m,o]*(mu[o]-g) ) }
    for k=1:K
        X1k = zeros(1, 'like', X);
        X1k = bsxfun(@plus,X1k,MU(im,k).');
        X1k = bsxfun(@plus,X1k,bsxfun(@minus, MU(io,k).', X(msk,io)) * (A(io,im,k) / A(im,im,k)));
        if sample
           %Smk = inv(A(im,im,k));
            X1k = X1k + chol(A(im,im,k))\randn(sum(im),Nm); % mvnrnd(zeros(1,Pm),Smk,Nm); %% mvnrnd is part of the statistics toolbox
        end
        X(msk,im) = X(msk,im) + bsxfun(@times, X1k, Z(msk,k));
    end
end

% =========================================================================
function logpX = marginal(X, cluster, const, L, E)
% logp = spm_gmm_lib('marginal', X, {MU,A},   const, [L], [E])
% logp = spm_gmm_lib('marginal', X, {MU,V,n}, const, [L], [E])
%
% X         - {NoxP} Observed values
% MU        -  PxK   (Expected) means
% A         -  PxPxK (Expected) precision matrices
% V         -  PxPxK Wishart scale matrices
% n         -  1xK   Wishart degrees of freedom
% const     -  MxK   Constant terms.
% L         -  MxP   Mask of missing patterns
% E         - {NoxP} Uncertainty (or 1xP)
%
% logpX     - {NoxK} (Expected) log-likelihood of belonging to each class
%
% Compute the expected log-likelihood of each observation belonging to each
% cluster: logpx(i,k) = E[ln p(g(i) | MU_k,A_k)]

%--------------------------------------------------------------------------
% Read input arguments
arraymode = ~iscell(X);
if ~iscell(X)
    X = {X};
    L = ones(1,size(X,2));
end
MU = double(cluster{1});
A  = double(cluster{2});
n  = [];
if numel(cluster) >= 3
    V = A;
    n = double(cluster{3});
end

% -------------------------------------------------------------------------
% Allocate output
P     = size(MU,1);
K     = size(MU,2);
logpX = cell(1,numel(X));
if nargin < 5, E = zeros(1,P); end
if isscalar(E), E = E*ones(1,P); end

% -------------------------------------------------------------------------
% For each combination of missing voxels
for i=1:size(L,1)

    % ---------------------------------------------------------------------
    % Get mask of missing values and modalities (with this particular code)
    Nm = size(X{i},1);          % Number of voxels with that code
    if Nm == 0, continue; end
    io = L(i,:);                % Observed channels
    Po = sum(io);               % Number of observed channels
    im = ~io;                   % Missing channels
    Pm = sum(im);               % Number of missing channels
    if Po == 0, continue; end
    if iscell(E), E1 = E{i}; else E1 = E(:,io); end % Uncertainty

    % ---------------------------------------------------------------------
    % Allocate logpX
    logpX{i} = zeros([Nm K], 'like', X{i});

    % ---------------------------------------------------------------------
    % Non constant terms
    for k=1:K

        % /!\ Sub-covariance is different from the inverse sub-precision
        % ML case:
        %   inv(S(o,o)) = A(o,o) - A(o,m)*A(m,m)\A(m,o)
        % Bayesian case:
        %   inv(S(o,o)) ~ W(V(o,o) - V(o,m)*V(m,m)\V(m,o), n - Pm)
        %   >> See Theorem 3.4.6 in:
        %      Mardia, K.V., Kent, J.T., Bibby, J.M., 1980.
        %      Multivariate Analysis, 1st edition. Academic Press.

        if sum(n) > 0
            Ao = V(io,io,k) - V(io,im,k)*(V(im,im,k)\V(im,io,k));
            Ao = (n(k)-Pm) * Ao;
        else
            Ao = A(io,io,k) - A(io,im,k)*(A(im,im,k)\A(im,io,k));
        end

        % Quadratic term in observed values: (obs-mean) x (obs-mean)
        l = bsxfun(@minus, X{i}, 2*MU(io,k)') * Ao;
        l = -0.5 * dot(l, X{i}, 2);

        % Binning uncertainty
        if ~isempty(E1) && any(any(E1))
            if size(E1,1)==1
                l = l - 0.5 * trace(diag(E1) * Ao);
            else
                l = l - 0.5 * sum(bsxfun(@times,diag(Ao)',E1),2);
            end
        end

        % Reshape as a column vector
        logpX{i}(:,k) = const(i,k) + l;
    end

end

if arraymode, logpX = logpX{1}; end

% =========================================================================
function [Z,lb] = responsibility(varargin)
% FORMAT [Z,lb] = spm_gmm_lib('responsibility', logpX, logPI, ...)
% logpX - {NoxK} Marginal log-likelihood
% logPI - {NoxK} Prior log-probabilities (or 1XK)
% ...   -        Other log-priors
%
% Compute responsibilities.
% Responsibilities are the posterior expected value of class-indexing
% vectors z_n.
% The posterior is computed as:
%   r_nk = exp(E[log Pi_k] + E[log p(x_n | Theta_k)]) / sum_k {r_nk}
%
% Extra terms can be added to the responsibilities prior to softmax by the
% varargin argument. These arguments need to be compatible (w.r.t. size) with
% the following function call: bsxfun(@plus, Z, varargin{i}).

arraymode = ~iscell(varargin{1});
if arraymode
    varargin{1} = varargin(1);
end

Z = cell(1,numel(varargin{1}));
[Z{:}] = deal(0);
for j=1:numel(varargin)
    Zj = varargin{j};
    if ~isempty(Zj)
        for i=1:numel(Z)
            if iscell(Zj), Zji = Zj{i};
            else           Zji = Zj; end
            Z{i} = bsxfun(@plus, Z{i}, Zji);
        end
    end
end

% Exponentiate and normalise
lb = 0;
for i=1:numel(Z)
    mx   = max(Z{i}, [], 2);
    Z{i} = bsxfun(@minus, Z{i}, mx);
    Z{i} = exp(Z{i});
    sz   = sum(Z{i}, 2);
    if nargout>=2
        lb = lb + sum(log(sz),1) + sum(mx,1);
    end
    Z{i} = bsxfun(@rdivide, Z{i}, sz);
end

if arraymode, Z = Z{1}; end

% =========================================================================
function varargout = suffstat(varargin)
% FORMAT [SS0,SS1,SS2] = spm_gmm_lib('suffstat', X, Z, W, [L])
% >> Compute sufficient statistics (per code)
% FORMAT [SS0,SS1,SS2] = spm_gmm_lib('suffstat', 'infer', SS0, SS1, SS2, {MU,A}, L)
% >> Compute expected sufficient statics (from per-code suff stat)
% FORMAT         [SS2] = spm_gmm_lib('suffstat', 'bin', E, Z, W, [L])
% >> Compute uncertainty-related statistics
%
% X    - {NoxP} Observed + Inferred values
% W    - {Nox1} Observation weights
% E    - {NoxP} Observation uncertainty (or 1xP)
% Z    - {NoxK} Responsibilities
% L    -   Mx1  List of missing codes
%
% SS0 - 1xK   0th order suff stat (sum of resp)
% SS1 - PxK   1st order suff stat (weighted sum of intensities)
% SS2 - PxPxK 2nd order suff stat (weighted sum of squared intensities)
%
% Compute sufficient statistics up to 2nd order, taking into account
% inferred values and their uncertainty.

if nargin == 0
    help spm_gmm_lib>suffstat
    error('Not enough argument. Type ''help spm_gmm_lib>suffstat'' for help.');
end
if ~ischar(varargin{1})
    id = 'base';
else
    id = varargin{1};
    varargin = varargin(2:end);
end
switch lower(id)
    case {'base'}
        if iscell(varargin{1})
            [varargout{1:nargout}] = suffstat_missing(varargin{:});
        else
            [varargout{1:nargout}] = suffstat_classic(varargin{:});
        end
    case {'infer'}
        [varargout{1:nargout}] = suffstat_infer(varargin{:});
    case {'uncertainty'}
        [varargout{1:nargout}] = suffstat_uncertainty(varargin{:});
    otherwise
        help spm_gmm_lib>suffstat
        error('Unknown function %s. Type ''help spm_gmm_lib>suffstat'' for help.', id)
end

% =========================================================================
function [SS0,SS1,SS2] = suffstat_classic(X, Z, W, varargin)
% FORMAT [SS0,SS1,SS2] = suffstat_classic(X, Z, W)
%
% X    - NxP Observed + Inferred values
% Z    - NxK Responsibilities
% W    - Nx1 Observation weights
%
% Compute sufficient statistics (up to 2nd order)

if nargin < 3, W = 1; end

%--------------------------------------------------------------------------
% Dimensions
N = size(X, 1);
P = size(X, 2);
K = size(Z, 2);

%--------------------------------------------------------------------------
% Weight responsibilities
Z = bsxfun(@times, Z, W);

%--------------------------------------------------------------------------
% Oth order
SS0 = sum(Z, 1, 'omitnan', 'double');

%--------------------------------------------------------------------------
% 1st order
SS1 = sum(bsxfun(@times, X, reshape(Z, [N 1 K])), 1, 'omitnan', 'double');
SS1 = reshape(SS1, [P K]);

%--------------------------------------------------------------------------
% 1nd order
SS2 = zeros(P,P,K, 'double');
for i=1:P
    SS2(i,i,:) = reshape(sum(bsxfun(@times, Z, X(:,i).^2),1,'omitnan', 'double'), [1 1 K]);
    for j=i+1:P
        SS2(i,j,:) = reshape(sum(bsxfun(@times, Z, X(:,i).*X(:,j)),1,'omitnan', 'double'), [1 1 K]);
        SS2(j,i,:) = SS2(i,j,:);
    end
end

% =========================================================================
function [SS0,SS1,SS2] = suffstat_missing(X, Z, W, L)
% FORMAT [{SS0},{SS1},{SS2}] = suffstat_missing(X, Z, W, L)
%
% X    - {NoxP} Observed + Inferred values
% Z    - {NoxK} Responsibilities
% W    - {Nox1} Observation weights
% L    -   MxP  Mask of missing channels
%
% Compute sufficient statistics for each missing pattern.

if nargin < 3, W = 1; end

%--------------------------------------------------------------------------
% Multiply responsibility with observation weight
for i=1:size(L,1)
    if iscell(W), Z{i} = bsxfun(@times, Z{i}, W{i});
    else         Z{i} = bsxfun(@times, Z{i}, W);    end
end

%--------------------------------------------------------------------------
% Dimensions
K = size(Z{1},2);

SS0 = cell(1,numel(X));
SS1 = cell(1,numel(X));
SS2 = cell(1,numel(X));

%--------------------------------------------------------------------------
% Sum missing data
for i=1:size(L,1)

    % ---------------------------------------------------------------------
    % Get mask of missing values and modalities (with this particular code)
    Nm = size(X{i},1);          % Number of voxels with that code
    if Nm == 0, continue; end
    io = L(i,:);                % Observed channels
    Po = sum(io);               % Number of observed channels
    if Po == 0, continue; end

    % ---------------------------------------------------------------------
    % Oth order moment
    SS0{i} = sum(Z{i}, 1, 'double') + eps;
    if nargout == 1, return; end

    % ---------------------------------------------------------------------
    % 1st and 2nd order moments
    SS1{i} = zeros(Po,K);
    SS2{i} = zeros(Po,Po,K);
    for k=1:K
        zk = double(Z{i}(:,k));
        for m=1:Po
            xm  = double(X{i}(:,m));
            zx  = zk.*xm;
            SS1{i}(m,k)   = sum(zx);
            SS2{i}(m,m,k) = sum(zx.*xm);
            for m1=(m+1):Po
                xm  = double(X{i}(:,m1));
                SS2{i}(m,m1,k) = sum(zx.*xm);
                SS2{i}(m1,m,k) = SS2{i}(m,m1,k);
            end
        end
    end
end

% =========================================================================
function [SS0,SS1,SS2] = suffstat_infer(lSS0, lSS1, lSS2, cluster, L)
% FORMAT [SS0,SS1,SS2] = suffstat_infer(SS0, SS1, SS2, {MU,A}, L)
%
% SS0  - {1xK}     Oth order sufficient statistics
% SS1  - {PoxK}    1st order sufficient statistics
% SS2  - {PoxPoxK} 2nd order sufficient statistics
% MU   - PxK       Clusters' mean
% A    - PxPxK     Clusters' precision matrix
% L    - MxP       Mask of observed channels
%
% Compute expected 1st/2nd order statistics.

MU = double(cluster{1});
A  = double(cluster{2});

%--------------------------------------------------------------------------
% Dimensions
P = size(MU,1);
K = size(MU,2);

SS0 = zeros(1,K, 'like', lSS0{1});
if nargout > 1
    SS1 = zeros(P,K, 'like', lSS1{1});
    if nargout > 2
        SS2 = zeros(P,P,K, 'like', lSS2{1});
    end
end

for k=1:K
    if nargout > 1
        ss1 = SS1(:,k);
        Ak  = A(:,:,k);
        if nargout > 2
            ss2 = SS2(:,:,k);
        end
    end

    for i=1:size(L,1)
        io = L(i,:);                % Observed channels
        Po = sum(io);               % Number of observed channels
        if Po == 0, continue; end
        im = ~io;                   % Missing channels

        % -----------------------------------------------------------------
        % 0th order moment
        SS0k   = lSS0{i}(k);
        SS0(k) = SS0(k) + SS0k;

        if nargout > 1
        % -----------------------------------------------------------------
        % 1st order moment
            SS1k = lSS1{i}(:,k);

            % 0) precompute stuff
            ss1o = ss1(io,:);
            ss1m = ss1(im,:);
            MUo  = MU(io,k);
            MUm  = MU(im,k);

            iAkmm = inv(Ak(im,im));
            SA    = iAkmm*Ak(im,io);

            % 1) observed
            ss1o = ss1o + SS1k;

            % 2) missing
            % > t = mu(m) + A(m,m) \ A(m,o) * (mu(o) - g)
            ss1m = ss1m + SS0k * MUm;
            ss1m = ss1m + SA * (SS0k * MUo - SS1k);

            ss1(io) = ss1o;
            ss1(im) = ss1m;
        end

        if nargout > 2
        % -----------------------------------------------------------------
        % 2nd order moment: quadratic terms
            SS2k = lSS2{i}(:,:,k);

            % 0) precompute stuff
            ss2oo = ss2(io,io);
            ss2mo = ss2(im,io);
            ss2mm = ss2(im,im);
            MUMUm = SS0k * (MUm * MUm.');
            MUMUo = SS0k * (MUo * MUo.');
            GMUo  = SS1k * MUo.';
            GMUm  = SS1k * MUm.';

            % 1) observed x observed
            ss2oo = ss2oo + SS2k;

            % 2) missing x observed
            ss2mo = ss2mo + GMUm.' + SA * (GMUo.' - SS2k);

            % 3) missing x missing
            tmp   = SA * (SS0k * MUo - SS1k) * MUm.';
            ss2mm = ss2mm + MUMUm ...
                          + (tmp+tmp') ...
                          + SA * (SS2k + MUMUo - GMUo.' - GMUo) * SA.';

            % 4) uncertainty ~ missing
            ss2mm = ss2mm + SS0k*iAkmm;

            ss2(io,io) = ss2oo;
            ss2(im,io) = ss2mo;
            ss2(io,im) = ss2mo.';
            ss2(im,im) = ss2mm;
        end
    end

    if nargout > 1
        SS1(:,k) = ss1;
        if nargout > 2
            SS2(:,:,k) = ss2;
        end
    end

end

% =========================================================================
function SS2 = suffstat_uncertainty(E, Z, W, L)

if iscell(Z)
    SS2 = suffstat_uncertainty_missing(E, Z, W, L);
else
    SS2 = suffstat_uncertainty_classic(E, Z, W);
end

% =========================================================================
function SS2 = suffstat_uncertainty_missing(E, Z, W, L)
% FORMAT SS2 = suffstat_uncertainty_missing(E, Z, W, L)
%
% E  {NoxP} - Variance in each modality due to binning
% Z  {NoxK} - Responisbilities
% W  {Nox1} - Observation weights
% L   MxP   - Mask of observed channel
%
% Compute "uncertainty" 2nd order statistics based on the posterior
% precision matrix about inferred values.

if nargin < 3 || isempty(W)
    W = 1;
end

%--------------------------------------------------------------------------
% Dimensions
K = size(Z{1},2);
P = size(L,2);
SS2 = zeros(P,P,K);
if sum(E) == 0, return; end

%--------------------------------------------------------------------------
% Multiply responsibility with observation weight
for i=1:size(L,1)
    if iscell(W), Z{i} = bsxfun(@times, Z{i}, W{i});
    else          Z{i} = bsxfun(@times, Z{i}, W);    end
end

% -------------------------------------------------------------------------
% 2nd order moment: uncertainty ~ binning
for i=1:size(L,1)
    io = L(i,:);                % Observed channels
    Po = sum(io);               % Number of observed channels
    if Po == 0, continue; end
    if iscell(E), E1 = E{i}; else E1 = E(:,io); end % Uncertainty

    list_p = 1:P;
    list_p = list_p(io);
    for p=list_p
        if size(E1,1)==1
            SS2(p,p,:) = SS2(p,p,:) ...
                + bsxfun(@times, reshape(sum(Z{i}, 1, 'double'), [1 1 K]), E1(p));
        else
            SS2(p,p,:) = SS2(p,p,:) ...
                + reshape(sum(bsxfun(@times, Z{i}, E1(:,p)), 1, 'double'), [1 1 K]);
        end
    end
end

% =========================================================================
function SS2 = suffstat_uncertainty_classic(E, Z, W, varargin)
% FORMAT SS2 = suffstat_uncertainty_classic(E, Z, W)
%
% E  NxP - Variance in each modality due to binning
% Z  NxK - Responisbilities
% W  Nx1 - Observation weights
%
% Compute "uncertainty" 2nd order statistics based on the posterior
% precision matrix about inferred values.

if nargin < 3 || isempty(W)
    W = 1;
end

%--------------------------------------------------------------------------
% Dimensions
K = size(Z,2);
P = size(E,2);
SS2 = zeros(P,P,K);
if sum(E) == 0, return; end

%--------------------------------------------------------------------------
% Multiply responsibility with observation weight
Z = bsxfun(@times, Z, W);

% -------------------------------------------------------------------------
% 2nd order moment: uncertainty ~ binning
for p=1:P
    if size(E,1)==1
        SS2(p,p,:) = SS2(p,p,:) ...
            + bsxfun(@times, reshape(sum(Z, 1, 'double'), [1 1 K]), E(p));
    else
        SS2(p,p,:) = SS2(p,p,:) ...
            + reshape(sum(bsxfun(@times, Z, E(:,p)), 1, 'double'), [1 1 K]);
    end
end

% =========================================================================
function c = normalisation(mean,prec,L)
% FORMAT c = spm_gmm_lib('normalisation', {MU,b}, {V,n})
% FORMAT c = spm_gmm_lib('normalisation', {MU,b}, {A})
% FORMAT c = spm_gmm_lib('normalisation', {MU},   {A})
% FORMAT c = spm_gmm_lib('normalisation', ..., L)
% MU - (Expected) mean
% b  - Mean df (if isempty or 0 -> no Bayesian prior)
% V  - Scale matrix     (if not n isempty or 0)
% A  - Precision matrix (if n isempty or 0)
% n  - Precision df (if isempty or 0 -> no Bayesian prior)
%
% L - If provided, compute one term for each combination of missing data
%
% Compute the constant term (w.r.t. voxels) of each Gaussian
% (expected) log-likelihood.

MU = [];
b  = [];
V  = []; % It can actually be A (when n == 0)
n  = [];

%--------------------------------------------------------------------------
% Read input arguments
if ~iscell(mean)
    MU = mean;
else
    if numel(mean) >= 1
        MU = mean{1};
        if numel(mean) >= 2
            b = mean{2};
        end
    end
end
if ~iscell(prec)
    V = prec;
else
    if numel(prec) >= 1
        V = prec{1};
        if numel(prec) >= 2
            n = prec{2};
        end
    end
end


%--------------------------------------------------------------------------
% Use double precision
MU = double(MU);
b  = double(b);
V  = double(V);
n  = double(n);

%--------------------------------------------------------------------------
% Dimensions
P = size(MU,1);
K = size(MU,2);

if nargin == 3 && ~isempty(L) && ~all(L(:)==1)
%--------------------------------------------------------------------------
% Marginal distribution (do not use missing dimensions)

% Assume that:
% . A is a KxK positive-definite matrix
% . A ~ W_K(V,n)
% . S = inv(A)
% . A is partitioned as [A11 A12: A12' A22] with A11 PxP
% . V is partitioned as [V11 V12: V12' V22] with V11 PxP
% . S is partitioned as [S11 S12: S12' S22] with S11 PxP
% Then
% . A11 ~ W_P(V11,n)
% . inv(S11) = A11 - A12*A22\A12'
% . inv(S11) ~ W_P(V11 - V12*V22\V12', n - K + P)
% This allows us to compute E[inv(S11)] and E[ln|S11]], which are needed to
% compute the expected marginal distribution within each cluster.

    c = zeros(size(L,1), K, 'like', MU);
    for i=1:size(L,1)
        io = L(i,:);                % Observed channels
        Po = sum(io);               % Number of observed channels
        im = ~io;                   % Missing channels
        Pm = sum(im);               % Number of missing channels
        for k=1:K
            Vo = V(io,io,k) - V(io,im,k)*(V(im,im,k)\V(im,io,k));
            c(i,k) = - 0.5 * Po * log(2*pi);
            if sum(n) > 0
                no = n(k) - Pm;
                c(i,k) = c(i,k) + 0.5 * wishart_elogdet(Vo,no) ...
                                - 0.5 * no * MU(io,k).' * Vo * MU(io,k);
            else
                c(i,k) = c(i,k) + 0.5 * logdet(Vo) ...
                                - 0.5 * MU(io,k).' * Vo * MU(io,k);
            end
            if sum(b) > 0
                c(i,k) = c(i,k) - 0.5 * Po / b(k);
            end
        end

    end

else
%--------------------------------------------------------------------------
% No missing dimensions
    c = zeros(1,K, 'like', MU);
    for k=1:K
        c(k) = - 0.5 * P * log(2*pi);
        if sum(n) > 0
            c(k) = c(k) + 0.5 * wishart_elogdet(V(:,:,k),n(k)) ...
                        - 0.5 * n(k) * MU(:,k)' * V(:,:,k) * MU(:,k);
        else
            c(k) = c(k) + 0.5 * logdet(V(:,:,k)) ...
                        - 0.5 * MU(:,k)' * V(:,:,k) * MU(:,k);
        end
        if sum(b) > 0
            c(k) = c(k) - 0.5 * P / b(k);
        end
    end
end

% =========================================================================
function [MU,A,b,V,n] = updateclusters(SS0,SS1,SS2,pr)
% FORMAT [MU,A,b,V,n] = updateclusters(SS0,SS1,SS2,{MU0,b0,V0,n0})
% SS0 - 0th order sufficient statistics (sum Z_i)
% SS1 - 1st order sufficient statistics (sum Z_i * X_i)
% SS2 - 2nd order sufficient statistics (sum Z_i * (X_i * X_i'))
% pr  - List of prior Gauss-Wishart parameters.
%
% Compute posterior GMM parameters from suff stats.

if nargin<4, pr =[]; end

K  = numel(SS0);
MU0 = [];
b0  = [];
V0  = [];
n0  = [];
if numel(pr) >= 1
    MU0 = pr{1};
    if numel(pr) >= 2
        b0 = pr{2};
        if numel(pr) >= 3
            V0 = pr{3};
            if numel(pr) >= 4
                n0 = pr{4};
            end
        end
    end
end

%--------------------------------------------------------------------------
% Use double precision
MU0 = double(MU0);
b0  = double(b0);
V0  = double(V0);
n0  = double(n0);

% -------------------------------------------------------------------------
% Mean
if sum(b0) == 0
    % ---------------------------------------------------------------------
    % Without prior
    b  = [];
    MU = bsxfun(@rdivide, SS1, SS0);
else
    % ---------------------------------------------------------------------
    % With prior
    b  = b0 + SS0;
    MU = bsxfun(@rdivide, SS1 + bsxfun(@times,b0,MU0), b);
end

% -------------------------------------------------------------------------
% Scale/Precision
V = zeros(size(SS2));
if sum(n0) == 0
    % ---------------------------------------------------------------------
    % Without prior
    n   = [];
    for k=1:K
        V(:,:,k) = inv(SS2(:,:,k) / SS0(k) - (MU(:,k) * MU(:,k).'));
    end
else
    % ---------------------------------------------------------------------
    % With prior
    n = n0 + SS0;
    for k=1:K
        V(:,:,k) = inv(SS2(:,:,k) + b0(k) * MU0(:,k) * MU0(:,k).' ...
                                  -  b(k) * MU(:,k)  * MU(:,k).' ...
                                  + inv(V0(:,:,k)));
    end
end
if sum(n) > 0
    A = bsxfun(@times, V, reshape(n, [1 1 K]));
else
    A = V;
end

% =========================================================================
function [PI,logPI,a] = updateproportions(SS0, a0)
% FORMAT [PI,logPI,a,ll] = updateproportions(SS0, a0)
%
% SS0 - 1xK 0th order sufficient statistics (sum of responsibilities)
% a0  - 1xK Dirichlet prior (can be 0)
%
% PI    - 1xK Cluster proportion posterior expected value
% logPI - 1xK ln(PI) or E[ln(PI)] (if Bayesian)
% a     - Dirichlet posterior (if Bayesian)
% ll    - 1x1 Lower bound: E[ln p(PI|a)] - E[ln q(PI)]
%
% Bayesian or ML update of cluster proportions.

a = a0 + SS0;
if sum(a0(:))
% Bayesian
    % expected values
    logPI = psi(a) - psi(sum(a));
    PI    = a ./ sum(a(:));
else
% Maximum Likelihood
    a     = max(a, eps);
    PI    = a ./ sum(a(:));
    logPI = log(PI);
end
% =========================================================================

% =========================================================================
function [GaussPrior,extras] = updatehyperpars(cluster,GaussPrior,varargin)
% FORMAT [GaussPrior,extras] = updatehyperpars(cluster,GaussPrior,varargin)
%
% REQUIRED
% cluster    - 1xS cell array where cluster{s} = {{MU,b},{V,n}}
% GaussPrior - {MU0,b0,V0,n0}
%
% OPTIONAL
% constrained - Optimise hierarchical prior on V [false]
% figname     - Postfix added to figure name ['']
% verbose     - Verbosity level: [false]=quiet, true=display
%
% OUTPUT
% GaussPrior - New {MU0,b0,V0,n0}
% extras     - Struct with lower bound information, etc.
%
% Update of VB-GMM hyper-parameters (m,b,W,n).

% Parse optional arguments
%--------------------------------------------------------------------------
p = inputParser;
p.FunctionName = 'updatehyperpars';
p.addParameter('constrained',0,@islogical);
p.addParameter('figname','',@ischar);
p.addParameter('verbose',0,@islogical);
p.addParameter('lkp',[],@isnumeric);
p.parse(varargin{:});
constrained = p.Results.constrained;
figname     = p.Results.figname;
verbose     = p.Results.verbose;
lkp         = p.Results.lkp;

% Parameters
S = numel(cluster); % Number of posteriors

m0 = GaussPrior{1};
b0 = GaussPrior{2};
W0 = GaussPrior{3};
n0 = GaussPrior{4};

N = size(m0,1);
K = size(m0,2);

% pre-allocate
LogDetW0  = zeros(size(n0));
V         = zeros(size(W0));
p         = zeros(size(n0));
p0        = 0;

% -------------------------------------------------------------------------
%   Gauss-Wishart "mean" parameters
% -------------------------------------------------------------------------

for k=1:K

    % ---------------------------------------------------------------------
    % Update m0 (mode, closed-form)

    Lambda   = 0;
    LambdaMu = 0;
    for s=1:S
        [m,~,W,n] = get_posteriors(cluster,s);
        Lambda    = Lambda   + n(k)*W(:,:,k);
        LambdaMu  = LambdaMu + n(k)*W(:,:,k)*m(:,k);
    end
    m0(:,k) = Lambda \ LambdaMu;

    % ---------------------------------------------------------------------


    % ---------------------------------------------------------------------
    % Update b0 (mode, closed-form)

    b0(k)= 0;
    for s=1:S
        [m,b,W,n] = get_posteriors(cluster,s);
        m1 = m(:,k) - m0(:,k);
        b0(k) = b0(k) + m1.' * (n(k)*W(:,:,k)) * m1 + N/b(k);
    end
    b0(k) = N*S/b0(k);

    % ---------------------------------------------------------------------

end


% =========================================================================
% NOT CONSTRAINED
if ~constrained

    % ---------------------------------------------------------------------
    %   Gauss-Wishart "precision" parameters
    % ---------------------------------------------------------------------

    for k=1:K

        % ---
        % Set up some constants
        sumLogDet = 0;
        sumPsi    = 0;
        Wn        = 0;
        for s=1:S
            [~,~,W,n] = get_posteriors(cluster,s);
            sumLogDet = sumLogDet + logdet( W(:,:,k));
            sumPsi    = sumPsi    + DiGamma(n(k)/2, N);
            Wn        = Wn        + n(k)*W(:,:,k);
        end
        sumLogDet = sumLogDet/S;
        sumPsi    = sumPsi/S;
        Wn        = Wn/S;

        % -----------------------------------------------------------------
        % Update n0 (mode, Gauss-Newton [convex])
        E     = inf;
        for gniter=1:1000

            % -------------------------------------------------------------
            % Update W0 (mode, closed-form)
            W0(:,:,k)   = Wn/n0(k);
            LogDetW0(k) = logdet( W0(:,:,k));
            % -------------------------------------------------------------

            % ---
            % Objective function
            Eprev = E;
            E = 0.5*S*n0(k)*( LogDetW0(k) - sumLogDet - sumPsi ) ...
                + S*LogGamma(n0(k)/2, N);

            if E == Eprev
                break;
            end

            % ---
            % Gradient & Hessian
            g = 0.5*S*( LogDetW0(k) - sumLogDet - sumPsi ...
                         + DiGamma(n0(k)/2, N) );
            H = S/4*DiGamma(n0(k)/2, N, 1);

            % ---
            % Update
            n0(k) = max(n0(k) - H\g, N-1+2*eps);

        end
        % -----------------------------------------------------------------

    end

    % ---------------------------------------------------------------------
    %   Save results
    % ---------------------------------------------------------------------
    extras.b   = b0;
    extras.m   = m0;
    extras.n   = n0;
    extras.W   = W0;
    extras.ldW = LogDetW0;
    extras.lb  = 0;

    GaussPrior{1} = m0;
    GaussPrior{2} = b0;
    GaussPrior{3} = W0;
    GaussPrior{4} = n0;

% =========================================================================
% CONSTRAINED
else

    lb = -inf;
    for em=1:50
        % ---
        % Starting estimate
        if p0 == 0
            p0 = 0;
            V0 = 0;
            for k=1:K
                p0 = p0 + S*n0(k);
                for s=1:S
                    [~,~,W,n] = get_posteriors(cluster,s);
                    V0 = V0 + inv(n(k)*W(:,:,k));
                end
            end
            p0 = p0/K;
            V0 = V0/K;
        end

        % -----------------------------------------------------------------
        %   Gauss-Wishart "precision" parameters
        % -----------------------------------------------------------------

        for k=1:K

            % ---
            % Set up some constants
            % > compute sum E[logdet W] and sum psi(nu/2)
            logDetW  = 0;
            psiN     = 0;
            Lambda   = 0;
            for s=1:S
                [~,~,W,n] = get_posteriors(cluster,s);
                logDetW = logDetW  + logdet( W(:,:,k));
                psiN    = psiN     + DiGamma(n(k)/2, N);
                Lambda  = Lambda   + n(k)*W(:,:,k);
            end
            logDetW  = logDetW/S;
            psiN = psiN/S;


            % -------------------------------------------------------------
            % Update n0 (mode, Gauss-Newton [convex])
            E = inf;

            % ---------------------------------------------------------
            % Update {p,V} for W0 (posterior, closed form)
            p(k)       = p0 + S*n0(k);
           %V(:,:,k)   = inv(inv(V0) + Lambda); % NO NEED TO DO THIS MULTIPLE TIMES
            V(:,:,k)   = (V0*Lambda + eye(size(Lambda)))\V0;

            for gniter=1:100

                % Useful values
                W0(:,:,k)   = inv(wishart_e(V(:,:,k), p(k)));
                LogDetW0(k) = -wishart_elogdet(V(:,:,k), p(k));
                % ---------------------------------------------------------

                % ---
                % Objective function
                E1 = S*n0(k)/2 * (LogDetW0(k) - logDetW - psiN) ...
                     + S*LogGamma(n0(k)/2, N);
                E = [E E1];

                subgain = get_gain(E);
                if subgain < 1e-6
                    % Finished
                    break
                end

                % ---
                % Gradient & Hessian
                g = S/2*(LogDetW0(k) - logDetW - psiN + DiGamma(n0(k)/2, N));
                H = S/4 * DiGamma(n0(k)/2, N, 1);

                % ---
                % Update
                n0(k) = max(n0(k) - H\g, N-1+2*eps);
            end
            % ------------------------------------------------------------

        end


        % -----------------------------------------------------------------
        %   Inverse-Wishart parameters
        % -----------------------------------------------------------------

        % ---
        % Set up some constants
        % > compute sum Logdet(psi) and sum psi(m/2)
        sumlogV = 0;
        sumPsi  = 0;
        pV      = 0;
        for k=1:K
            sumlogV = sumlogV + logdet( V(:,:,k));
            sumPsi  = sumPsi  + DiGamma(p(k)/2, N);
            pV      = pV      + p(k)*V(:,:,k);
        end
        sumlogV = sumlogV/K;
        sumPsi  = sumPsi/K;
        pV      = pV/K;


        % -----------------------------------------------------------------
        % Update p0 (mode, Gauss-Newton [convex])
        E = inf;
        for gniter=1:1000

            % -------------------------------------------------------------
            % Update V0 (closed-form)
            V0 = pV/p0;
            LogDetV0 = logdet( V0);
            % -------------------------------------------------------------

            % ---
            % Objective function
            Eprev = E;
            E = p0*K/2*( N*LogDetV0 - sumlogV - sumPsi ) + K*LogGamma(p0/2, N);
            if E == Eprev
                break;
            end

            % ---
            % Gradient & Hessian
            g = K/2*( LogDetV0 - sumlogV - sumPsi + DiGamma(p0/2, N) );
            H = K/4*DiGamma(p0/2, N, 1);

            % ---
            % Update
            p0 = max(p0 - H\g, N-1+2*eps);

        end
        % -----------------------------------------------------------------


        % ---
        % Objective function
        nlb  = 0;
        for k=1:K
            nlb  = nlb - wishart_kl(V(:,:,k), p(k), V0, p0);
        end

        lb   = [lb nlb];
        gain = get_gain(lb);

        if gain < 1e-3
            % Finished
            break
        end

    end % < "EM" loop

    % ---------------------------------------------------------------------
    %   Save results
    % ---------------------------------------------------------------------
    extras.b   = b0;
    extras.m   = m0;
    extras.n   = n0;
    extras.W   = W0;
    extras.ldW = LogDetW0;
    extras.V   = V;
    extras.p   = p;
    extras.V0  = V0;
    extras.p0  = p0;
    extras.lb  = 0;
    for k=1:K
        extras.lb  = extras.lb - wishart_kl(V(:,:,k), p(k), V0, p0);
    end

    GaussPrior{1} = m0;
    GaussPrior{2} = b0;
    GaussPrior{3} = W0;
    GaussPrior{4} = n0;
end

if verbose
    % Visualise results
    plot_GaussPrior(GaussPrior,lkp,figname);
end
% =========================================================================

%--------------------------------------------------------------------------
% Lower bound
% -------------------------------------------------------------------------

% =========================================================================
function [lb,const] = marginalsum(SS0, SS1, SS2, mean, prec, L, SS2b)
% [lb,const] = spm_gmm_lib('marginalsum', SS0, SS1, SS2, MU, A, L, SS2b)
% [lb,const] = spm_gmm_lib('marginalsum', SS0, SS1, SS2, {MU,b}, {V,n}, L, SS2b)
%
% SS0       - {1xK}   Zero-th order moment (per config)
% SS1       - {PxK}   First   order moment (per config)
% SS1       - {PxPxK} Second  order moment (per config)
% MU        - PxK     Means
% A/V       - PxPxK   Precision/Scale matrices
% b         - 1xK     Mean degrees of freedom
% n         - 1xK     Precision degrees of freedom
% L         - Mx1     List of existing codes
% SS2b      - PxPxK   Binning uncertainty
%
% lb        -         Sum of (expected) marginal likelihoods
% const     - MxK     Constant terms
%
% Compute the expected log-likelihood of each observation belonging to each
% cluster: lb = sum_{i,k} E[z_ik] E[ln p(g(i) | MU_k,A_k)]


MU = [];
A  = [];
V  = [];
n  = [];

%--------------------------------------------------------------------------
% Read input arguments
if ~iscell(mean)
    MU = mean;
else
    if numel(mean) >= 1
        MU = mean{1};
        if numel(mean) >= 2
            b = mean{2};
        end
    end
end
if ~iscell(prec)
    A = prec;
else
    if numel(prec) >= 1
        A = prec{1};
        if numel(prec) >= 2
            V = A;
            n = prec{2};
        end
    end
end
if nargin < 6
    L = [];
end
if nargin < 7
    SS2b = 0;
end

% -------------------------------------------------------------------------
% Dimensions
P  = size(MU,1);
K  = size(MU,2);
if isempty(L), L = ones(1,P); end % None missing

% -------------------------------------------------------------------------
% Constant term
const = normalisation(mean, prec, L);

lb  = 0;

% -------------------------------------------------------------------------
% For each combination of missing voxels
for i=1:size(L,1)

    % ---------------------------------------------------------------------
    % Get mask of missing values and modalities (with this particular code)
    io = L(i,:);                % Observed channels
    Po = sum(io);               % Number of observed channels
    if Po == 0, continue; end
    im = ~io;                   % Missing channels
    Pm = sum(im);               % Number of missing channels

    % ---------------------------------------------------------------------
    % Initialise with constant term
    lb = lb + sum(const(i,:) .* SS0{i}, 'double');

    % ---------------------------------------------------------------------
    % Non constant terms
    for k=1:K

        % /!\ Sub-covariance is different from the inverse sub-precision
        % ML case:
        %   inv(S(o,o)) = A(o,o) - A(o,m)*A(m,m)\A(m,o)
        % Bayesian case:
        %   inv(S(o,o)) ~ W(V(o,o) - V(o,m)*V(m,m)\V(m,o), n - Pm)
        if sum(n) > 0
            Ao = V(io,io,k) - V(io,im,k)*(V(im,im,k)\V(im,io,k));
            Ao = (n(k)-Pm) * Ao;
        else
            Ao = A(io,io,k) - A(io,im,k)*(A(im,im,k)\A(im,io,k));
        end

        % 1) obs x mean
        lb = lb + SS1{i}(:,k).' * Ao * MU(io,k);

        % 1) obs x obs
        lb = lb - 0.5 * trace(Ao * SS2{i}(:,:,k));

        % 3) Binning uncertainty
        lb = lb - 0.5 * trace(Ao * SS2b);

    end

end

% =========================================================================
function varargout = kl(varargin)
% Useful KL-divergences for Gaussian Mixture modelling
%
% [klMU,klA] = spm_gmm_lib('kl', 'GaussWishart', {MU,b}, {V,n}, {MU0,b0}, {V0,n0})
% > KL divergence between two Gauss-Wishart distributions
%
% klP = spm_gmm_lib('kl', 'Dirichlet', a, a0)
% > KL divergence between two Dirichlet distributions
%
% klZ = spm_gmm_lib('kl', 'Categorical', Z, W, logPI, labels, logmg_w)
% > KL divergence between two Categorical distributions

if nargin == 0
    help spm_gmm_lib>kl
    error('Not enough argument. Type ''help spm_gmm_lib>kl'' for help.');
end
id = varargin{1};
varargin = varargin(2:end);
switch lower(id)
    case {'gausswishart','gw'}
        [varargout{1:nargout}] = kl_gausswishart(varargin{:});
    case {'dirichlet','d'}
        [varargout{1:nargout}] = kl_dirichlet(varargin{:});
    case {'categorical','cat','c'}
        [varargout{1:nargout}] = kl_categorical(varargin{:});
    otherwise
        help spm_gmm_lib>kl
        error('Unknown function %s. Type ''help spm_gmm_lib>kl'' for help.', id)
end

% =========================================================================
function klZ = kl_categorical(Z, W, logPI, labels, logmg_w)
if nargin < 4, labels = 0; end
if nargin < 5, logmg_w = 0; end

if ~iscell(Z)
    Z = {Z};
end

% Initialise
klZ = zeros(1, 'like', Z{1});

for i=1:numel(Z)
    Z1 = Z{i};
    if iscell(W), W1 = W{i}; else W1 = W; end
    if iscell(logPI), logPI1 = logPI{i}; else logPI1 = logPI; end
    if iscell(labels), labels1 = labels{i}; else labels1 = labels; end

    if ~isempty(labels1)
        % E[ln p(Z|labels)]
        klZ = klZ + sum(sum(bsxfun(@times,Z1,labels1), 2) .* W1, 'double');
    end

    % E[ln p(Z|PI)] (prior ~ responsibilities)
    klZ = klZ + sum(sum(bsxfun(@times,Z1,bsxfun(@plus,logPI1,logmg_w)), 2) .* W1, 'double');

    % -E[ln q(Z)] (posterior ~ responsibilities))
    klZ = klZ - sum(sum(Z1 .* log(max(Z1,eps)), 2) .* W1, 'double');
end

% =========================================================================
function klP = kl_dirichlet(a, a0)

klP = zeros(1, 'like', a);
K   = numel(a);
if sum(a0) > 0
    % prior
    klP = gammaln(sum(a0)) - sum(gammaln(a0));
    klP = klP + sum((a0-1) .* (psi(a) - K*psi(sum(a))));
    % posterior
    klP = klP - gammaln(sum(a)) - sum(gammaln(a));
    klP = klP - sum((a-1) .* (psi(a) - K*psi(sum(a))));
end

% =========================================================================
function [klMU,klA] = kl_gausswishart(mean,prec,mean0,prec0)

MU  = [];
b   = [];
V   = []; % It can actually be A (when n == 0)
n   = [];
MU0 = [];
b0  = [];
V0  = [];
n0  = [];

%--------------------------------------------------------------------------
% Read input arguments
if ~iscell(mean)
    MU = mean;
else
    if numel(mean) >= 1
        MU = mean{1};
        if numel(mean) >= 2
            b = mean{2};
        end
    end
end
if ~iscell(prec)
    V = prec;
else
    if numel(prec) >= 1
        V = prec{1};
        if numel(prec) >= 2
            n = prec{2};
        end
    end
end
if nargin >= 3
    if ~iscell(mean0)
        MU0 = mean0;
    else
        if numel(mean0) >= 1
            MU0 = mean0{1};
            if numel(mean0) >= 2
                b0 = mean0{2};
            end
        end
    end
end
if nargin >=4
    if ~iscell(prec0)
        V0 = prec0;
    else
        if numel(prec0) >= 1
            V0 = prec0{1};
            if numel(prec0) >= 2
                n0 = prec0{2};
            end
        end
    end
end

%--------------------------------------------------------------------------
% Read input arguments
P = size(MU,1);
K = size(MU,2);
LogDetA = zeros(1,K, 'like', V);
if sum(n) > 0
    A = bsxfun(@times, V, reshape(n, [1 1 K]));
    for k=1:K
        LogDetA(k) = wishart_elogdet(V(:,:,k),n(k));
    end
else
    A = V;
    for k=1:K
        LogDetA(k) = logdet(A(:,:,k));
    end
end

% Lower bound
klMU = zeros(1, 'like', MU);
klA  = zeros(1, 'like', A);
for k=1:K
    % + prior
    if sum(b0) > 0
        % prior
        klMU = klMU - P*log(2*pi) ...
                    + P*log(b0(k)) ...
                    + LogDetA(k) ...
                    - b0(k)*(MU(:,k)-MU0(:,k)).'*A(:,:,k)*(MU(:,k)-MU0(:,k)) ...
                    - P*b0(k)/b(k);
        % posterior
        klMU = klMU + P*log(2*pi) ...
                    - P*log(b(k)) ...
                    - LogDetA(k) ...
                    + P;
    end
    if sum(n0) > 0
        klA = klA - wishart_kl(V(:,:,k), n(k), V0(:,:,k), n0(k));
    end
end
klMU = 0.5 * klMU;

%--------------------------------------------------------------------------
% "Missing code" image
% -------------------------------------------------------------------------

% =========================================================================
function [Xo,C,L] = obs2cell(X,C,prune)
% FORMAT [Xo,C,L] = obs2cell(Xi,[C],[prune])
% Xi - NxP          observation matrix
% prune -           only extract observed channels [true]
% Xo - Mx{NoxPo}    cell of observation matrices per missing code
% C  - Nx1          missing code image
% L  - MxP          mask of observed channels per code
%
% Prepare input observation matrix for processing through GMM loop with
% missing data handling.

if nargin < 3
    prune = true;
end
if nargin < 2
    C  = obs2code(X);
end
codes = unique(C);
codes = codes(codes ~= 0);
Xo = cell(1,numel(codes));
L  = false(numel(codes),size(X,2));
for i=1:numel(codes)
    msk    = (C == codes(i));
    io     = code2bin(codes(i),size(X,2));
    L(i,:) = io;
    if prune
        Xo{i}  = X(msk,io);
    else
        Xo{i}  = X(msk,:);
    end
end

% =========================================================================
function Xo = cell2obs(X,C,L)
% FORMAT Xo = cell2obs(Xi,C,L)
% Xi - Mx{NoxPo}    cell of observation matrices per missing code
% C  - Nx1          missing code image
% L  - MxP          mask of observed channels per code
% Xo - NxP          observation matrix

ispruned = true;
for i=1:numel(X)
    ispruned = ispruned && (size(X{i},2) == sum(L(i,:)));
    if ~ispruned, break; end
end

if ispruned
    Xo = NaN(numel(C),size(L,2), 'like', X{1});
else
    Xo = NaN(numel(C),size(X{1},2), 'like', X{1});
end
for i=1:size(L,1)
    io   = L(i,:);
    code = bin2code(io);
    msk  = (C == code);
    if ispruned
        Xo(msk,io) = X{i};
    else
        Xo(msk,:) = X{i};
    end
end

% =========================================================================
function code = obs2code(X)
% FORMAT code = obs2code(X)
%
% Compute a "missing code" image from the input observation matrix.

code = bin2code(~isnan(X));

% =========================================================================
function code = bin2code(X)
% FORMAT code = bin2code(X)
%
% Compute a "missing code" image from a mask of observed channels.

type = str2func(nbits2type(size(X,2)));
code = type(0);
for n=1:size(X,2)
    bin  = type(X(:,n));
    code = bitor(code, bitshift(bin, n-1));
end

% =========================================================================
function bin = code2bin(code, length)
% FORMAT bin = code2bin(code, length)
%
% Convert a "missing code" to a mask of observed channels

base = uint64(2).^uint64(0:(length-1));
bin  = bitand(uint64(code),base) > 0;

% =========================================================================
function type = nbits2type(nbits,signed)
% FORMAT type = nbits2type(nbits,signed)
%
% Find the best suited integer type to hold a given number of bits.
if nargin < 2
    signed = false;
end
type = 'int';
if ~signed
    type   = ['u' type];
end
if nbits <= 8
    type = [type '8'];
elseif nbits <= 16
    type = [type '16'];
elseif nbits <= 32
    type = [type '32'];
elseif nbits <= 64
    type = [type '64'];
else
    type = 'double';
    warning('No exact integer type for nbits > %d. I''ll use double instead', nbits);
end

% =========================================================================
function varargout = gmmplot(varargin)
% Custom visualisation tools for Gaussian Mixture modelling
%
% spm_gmm_lib('plot', 'lb', lb, (wintitle))
% > Plot lower bound
%
% spm_gmm_lib('plot', 'gmm', {X,W}, {MU,A}, PI, (wintitle))
% > Plot mixture fit
%
% spm_gmm_lib('plot', 'cat', dm, Z, Template, (wintitle))
% > Plot (categorical) responsibilities and template (if available)
%
% spm_gmm_lib('plot', 'gaussprior', GaussPrior, (wintitle))
% > Plot VB-GMM hyper-parameters
%
% c = spm_gmm_lib('plot', 'cat2rgb', f, pal)
% > Generate an RGB volume from a categorical (e.g. responsibilities) volume.
%
% c = spm_gmm_lib('plot', 'show_cat_img', img, title_nam)
% > Show categorical images
%

if nargin == 0
    help spm_gmm_lib>gmmplot
    error('Not enough argument. Type ''help spm_gmm_lib>gmmplot'' for help.');
end
id = varargin{1};
varargin = varargin(2:end);
switch lower(id)
    case {'lowerbound','lb'}
        [varargout{1:nargout}] = plot_lowerbound(varargin{:});
    case {'gmm'}
        [varargout{1:nargout}] = plot_gmm(varargin{:});
    case {'cat'}
        [varargout{1:nargout}] = plot_cat(varargin{:});
    case {'gaussprior'}
        [varargout{1:nargout}] = plot_GaussPrior(varargin{:});
    case {'cat2rgb'}
        [varargout{1:nargout}] = cat2rgb(varargin{:});
    case {'showcatimg'}
        [varargout{1:nargout}] = show_cat_img(varargin{:});
    otherwise
        help spm_gmm_lib>plot
        error('Unknown function %s. Type ''help spm_gmm_lib>gmmplot'' for help.', id)
end

% =========================================================================
function plot_lowerbound(lb, figname)

% -------------------------------------------------------------------------
% Get figure (create if it does not exist)
if nargin < 2
    figname = '(SPM) Plot GMM Lower Bound';
end
f = findobj('Type', 'Figure', 'Name', figname);
if isempty(f)
    f = figure('Name', figname, 'NumberTitle', 'off');
end
set(0, 'CurrentFigure', f);
clf(f);

% -------------------------------------------------------------------------
% Choose type
if isfield(lb, 'B')
    nrow = 2;
    ncol = 4;
else
    nrow = 2;
    ncol = 3;
end

% -------------------------------------------------------------------------
% Plots
subplot(nrow, ncol, sub2ind([ncol nrow], 1, 1));
plot(lb.sum)
title('Lower Bound')
subplot(nrow, ncol, sub2ind([ncol nrow], 2, 1));
if isfield(lb, 'B')
    plot(sum(lb.X,1) + sum(lb.XB,1));
else
    plot(sum(lb.X,1))
end
box on
title('Observations (E)')
subplot(nrow, ncol, sub2ind([ncol nrow], 3, 1));
plot(sum(lb.Z,1))
box on
title('Responsibilities (KL)')
subplot(nrow, ncol, sub2ind([ncol nrow], 1, 2));
plot(sum(lb.P,1))
box on
title('Proportions (KL)')
subplot(nrow, ncol, sub2ind([ncol nrow], 2, 2));
plot(sum(lb.MU,1))
box on
title('Means (KL)')
subplot(nrow, ncol, sub2ind([ncol nrow], 3, 2));
plot(sum(lb.A,1))
box on
title('Precisions (KL)')
if isfield(lb, 'B')
    subplot(nrow, ncol, sub2ind([ncol nrow], 4, 2));
    plot(sum(lb.B,1))
    box on
    title('Bias Prior')
end
drawnow

% =========================================================================
function plot_cat(Z,Template,ticklabels,figname)
if nargin < 3, ticklabels = {}; end

% -------------------------------------------------------------------------
% Get figure (create if it does not exist)
if nargin<4
    figname = '(SPM) Plot GMM Categorical';
end
f = findobj('Type', 'Figure', 'Name', figname);
if isempty(f)
    f = figure('Name', figname, 'NumberTitle', 'off');
end
set(0, 'CurrentFigure', f);
clf(f);

show_cat_img({Z,Template},{'Z','Template'},ticklabels);
% =========================================================================

% =========================================================================
function plot_gmm(X, W, L, cluster, PI, part, figname)
% spm_gmm_lib('plot', 'gmm', X, W, L, {MU,A}, PI, part)
% X        {NoxPo}  Observations (or NxP)
% W        {Nox1}   Weights      (or Nx1 or 1x1)
% L         MxP     Mask of observed channels per code [ones]
% MU        PxK     Means
% A         PxPxK   Precisions
% PI        1xC     Class proportions
% part.lkp  1xC     Mapping from C classes to K clusters
% part.mg   1xK     Cluster proportions

MU = cluster{1};
A  = cluster{2};

if ~iscell(X), X = {X}; end
if ~iscell(W), W = {W}; end

% -------------------------------------------------------------------------
% Get figure (create if it does not exist)
if nargin < 7
    figname = '(SPM) Plot GMM';
end
f = findobj('Type', 'Figure', 'Name', figname);
if isempty(f)
    f = figure('Name', figname, 'NumberTitle', 'off');
end
set(0, 'CurrentFigure', f);
clf(f);

% -------------------------------------------------------------------------
% Sizes / colors
P = size(MU, 1);
K = size(MU, 2);

if nargin < 6
    lkp = 1:K;
    mg  = ones(1,K);
else
    lkp = part.lkp;
    mg  = part.mg;
end

PI = mg.*PI(:,lkp);

% Set colors
colors = hsv(max(lkp));

% -------------------------------------------------------------------------
% Bin data (if needed)
centres = cell(1,P);
weights = cell(1,P);
for p=1:P
    X1 = [];
    W1 = [];
    for i=1:numel(X)
        pp = [];
        if isempty(L)
            pp = p;
        elseif L(i,p)
            pp = 1:P;
            pp = pp(L(i,:));
            pp = find(pp == p);
        end
        if ~isempty(pp)
            X1 = [X1; X{i}(:,pp)];
            if iscell(W)
                iw = min(i,numel(W));
                if ~isscalar(W{iw})
                    W1 = [W1; W{iw}(:)];
                else
                    W1 = [W1; W{iw}*ones(size(X{i},1),1)];
                end
            end
        end
    end
    [centres{p},weights{p},~,width] = spm_histN(X1, 64, 'Weights', W1);
    weights{p} = weights{p}./(sum(weights{p})*width);
    clear X1
end

% -------------------------------------------------------------------------
% For each input dimension
for p=1:P
    % ---------------------------------------------------------------------
    % Plot histogram and marginal density
    subplot(2, P, p)
    hold on
    % ---------
    % Histogram
    bar(centres{p}, weights{p}, 'EdgeColor', 'none', 'FaceColor', [0.7 0.7 0.7]);
    ymax = max(weights{p});
    xlims = [inf -inf];
    % -----------
    % GMM Density
    for k=1:K
        x = linspace(MU(p,k)-3*A(p,p,k)^(-0.5),MU(p,k)+3*A(p,p,k)^(-0.5),100);
        y = PI(k)*spm_Npdf(x, MU(p,k), A(p,p,k)^(-1));
        plot(x, y, 'Color', colors(lkp(k),:), 'LineWidth', 1)
        xlims = [min([xlims(1) x]) max([xlims(2) x])];
    end
    xlabel(sprintf('x%d',p))
    ylabel('density')
    xlim(xlims);
    if ymax == 0
        ylim([0 1.1]);
    else
        ylim([0 1.1*ymax]);
    end
    box on
    hold off

    % ---------------------------------------------------------------------
    % Plot joint density (X1 vs Xj)
    if p > 1
        subplot(2, P, P+p)
        hold on
        for k=1:K
            Mu1     = MU([1 p],k);
            Sigma2  = inv( A([1 p],[1 p],k));
            Sigma   = sqrt(Sigma2);
            [x1,x2] = meshgrid(linspace(Mu1(1)-3*Sigma(1,1),Mu1(1)+3*Sigma(1,1),100)', ...
                               linspace(Mu1(2)-3*Sigma(2,2),Mu1(2)+3*Sigma(2,2),100)');
            y = spm_mvNpdf([x1(:) x2(:)]', Mu1, Sigma2);
            contour(x2, x1, reshape(y, [100 100])', 1, 'color', colors(lkp(k),:), 'LineWidth', 1);
        end
        xlabel(sprintf('x%d',p))
        ylabel('x1')
        xlim(xlims); % use same scale as histogram plot for comparison
        box on
        hold off
    end
end

% -------------------------------------------------------------------------
% Plot proportions in the remaining square
subplot(2, P, P+1)
hold on
for k=1:K
    bar(k, PI(k), 'FaceColor', colors(lkp(k),:));
end
xlabel('class')
ylabel('proportion')
box on
hold off
drawnow
% =========================================================================

% =========================================================================
function plot_GaussPrior(GaussPrior,lkp,figname)

if nargin < 2, lkp = []; end

if nargin==3
    figname0 = figname;
else
    figname0 = '(SPM) GaussPrior';
end

% ---------------------------------------------------------------------
% Get figure (create if it does not exist)
f = findobj('Type', 'Figure', 'Name', figname0);
if isempty(f)
    f = figure('Name', figname0, 'NumberTitle', 'off');
end
set(0, 'CurrentFigure', f);
clf(f);

m0 = GaussPrior{1};
W0 = GaussPrior{3};
n0 = GaussPrior{4};

P = size(m0,1);
K = size(m0,2);

if isempty(lkp)
    lkp = 1:K;
end

colors = hsv(max(lkp));

MU = m0;
A  = bsxfun(@times,W0,reshape(n0,[1 1 K]));

% ---------------------------------------------------------------------
% For each input dimension
for p=1:P
    % -----------------------------------------------------------------
    % Plot histogram and marginal density
    if P>1
        subplot(2, P, p)
    end
    hold on

    xlims = [inf -inf];
    % -----------
    % GMM Density
    for k=1:K
        x = linspace(MU(p,k)-3*A(p,p,k)^(-0.5),MU(p,k)+3*A(p,p,k)^(-0.5),100);
        y = 1/K*spm_Npdf(x, MU(p,k), A(p,p,k)^(-1));
        plot(x, y, 'Color', colors(lkp(k),:), 'LineWidth', 1)
        xlims = [min([xlims(1) x]) max([xlims(2) x])];
    end
    xlabel(sprintf('x%d',p))
    ylabel('density')
    xlim(xlims);
    box on
    hold off

    % -----------------------------------------------------------------
    % Plot joint density (X1 vs Xj)
    if p > 1
        subplot(2, P, P+p)
        hold on
        for k=1:K
            Mu1     = MU([1 p],k);
            Sigma2  = inv( A([1 p],[1 p],k));
            Sigma   = sqrt(Sigma2);
            [x1,x2] = meshgrid(linspace(Mu1(1)-3*Sigma(1,1),Mu1(1)+3*Sigma(1,1),100)', ...
                               linspace(Mu1(2)-3*Sigma(2,2),Mu1(2)+3*Sigma(2,2),100)');
            Sigma2(1,2) = Sigma2(2,1); % Make sure symmetric, spm_mvNpdf complains otherwise - but it is just numerical accuracy
            y = spm_mvNpdf([x1(:) x2(:)]', Mu1, Sigma2);
            contour(x2, x1, reshape(y, [100 100])', 1, 'color', colors(lkp(k),:), 'LineWidth', 1);
        end
        xlabel(sprintf('x%d',p))
        ylabel('x1')
        xlim(xlims); % use same scale as histogram plot for comparison
        box on
        hold off
    end
end

drawnow
%==========================================================================

%==========================================================================
function c = cat2rgb(f, pal)
% FORMAT c = cat2rgb(f, pal)
% f   - categorical (4D) image.
% pal - palette (Mx3 array or handle to palette function) [hsv]
%
% Generate an RGB volume from a categorical (e.g. responsibilities) volume.

if nargin < 2
    pal = @hsv;
end

if size(f,3)>1
    z = floor(size(f,3)/2) + 1;
    f = f(:,:,z,:);
end

tri = false;
if numel(size(f)) == 4 && size(f, 3) == 1
    tri = true;
    dm  = [size(f) 1 1];
    f   = reshape(f, [dm(1:2) dm(4)]);
end
if isa(pal, 'function_handle')
    pal = pal(size(f,3));
end

dm = [size(f) 1 1];
c  = zeros([dm(1:2) 3]); % output RGB image
s  = zeros(dm(1:2));     % normalising term

for k=1:dm(3)
    s = s + f(:,:,k);
    color = reshape(pal(k,:), [1 1 3]);
    c = c + bsxfun(@times, f(:,:,k), color);
end
if dm(3) == 1
    c = c / max(1, max(s(:)));
else
    c = bsxfun(@rdivide, c, s);
end

if tri
    c = reshape(c, [size(c, 1) size(c, 2) 1 size(c, 3)]);
end
%==========================================================================

%==========================================================================
function show_cat_img(img,title_nam,ticklabels)
if isnumeric(img)
    img = {img};
end
N       = numel(img);
if nargin < 2, title_nam = cell(1,N); end

dm0    = size(img{1});
K      = dm0(4);
colors = hsv(K);

if nargin < 3 || isempty(ticklabels), ticklabels = 1:K; end

if numel(ticklabels) ~= K
    ticklabels = 1:K;
end

if dm0(3)==1
    % 2d
    %----------------------------------------------------------------------

    for n=1:N

        subplot(1,N,n)

        slice = img{n}(:,:,1,:);
        slice = spm_gmm_lib('plot', 'cat2rgb', slice, colors);
        slice = squeeze(slice(:,:,:,:));
        slice = permute(slice,[2 1 3]);
        imagesc(slice); axis off image xy;

        if ~isempty(title_nam{n})
            title(title_nam{n})
        end
    end

    colormap(colors);
    cb = colorbar;
    set(gca, 'clim', [0.5 K+0.5]);
    set(cb, 'ticks', 1:K, 'ticklabels', ticklabels);
else
    % 3d
    %----------------------------------------------------------------------

    for n=1:N
        subplot(N,3,3*(n - 1) + 1)

        slice = img{n}(:,:,floor(dm0(3)/2) + 1,:);
        slice = spm_gmm_lib('plot', 'cat2rgb', slice, colors);
        slice = squeeze(slice(:,:,:,:));
        slice = permute(slice,[2 1 3]);
        imagesc(slice); axis off image xy;

        subplot(N,3,3*(n - 1) + 2)

        slice = permute(img{n}(:,floor(dm0(2)/2) + 1,:,:),[3 1 2 4]);
        slice = spm_gmm_lib('plot', 'cat2rgb', slice, colors);
        slice = squeeze(slice(:,:,:,:));
        imagesc(slice); axis off image xy;

        if ~isempty(title_nam{n})
            title(title_nam{n})
        end

        subplot(N,3,3*(n - 1) + 3)

        slice = permute(img{n}(floor(dm0(1)/2) + 1,:,:,:),[2 3 1 4]);
        slice = spm_gmm_lib('plot', 'cat2rgb', slice, colors);
        slice = squeeze(slice(:,:,:,:));
        slice = permute(slice,[2 1 3]);
        imagesc(slice); axis off image xy;

        colormap(colors);
        cb = colorbar;
        set(gca, 'clim', [0.5 K+0.5]);
        set(cb, 'ticks', 1:K, 'ticklabels', ticklabels);
    end
end

drawnow
%==========================================================================

% =========================================================================
function varargout = gmm_extras(varargin)
% Additional GMM functions
%
% gmm = spm_gmm_lib('extras', 'more_gmms', gmm, part)
% > A crude heuristic to replace a single Gaussian by a bunch of Gaussians.
% gmm = spm_gmm_lib('extras', 'collapse_gmms', gmm, part)
% > A crude heuristic to replace a bunch of Gaussian with a single Gaussian.
%

if nargin == 0
    help spm_gmm_lib>extras
    error('Not enough argument. Type ''help spm_gmm_lib>plot'' for help.');
end
id = varargin{1};
varargin = varargin(2:end);
switch lower(id)
    case {'more_gmms'}
        [varargout{1:nargout}] = more_gmms(varargin{:});
    case {'collapse_gmms'}
        [varargout{1:nargout}] = collapse_gmms(varargin{:});
    otherwise
        help spm_gmm_lib>extras
        error('Unknown function %s. Type ''help spm_gmm_lib>extras'' for help.', id)
end
% =========================================================================

% =========================================================================
function [gmm,mg_w] = more_gmms(gmm,lkp)
% FORMAT gmm = spm_gmm_lib('extras', 'more_gmms', gmm, lkp)
%
% gmm  - Cell with the following format {m,b,W,n}, where there are K
%        Gaussians.
% lkp - [1,K_p] vector partitioning a K GMM into a K_p GMM (K_P>=K). E.g.
%        [1 1 1 2 3 4 5 6 6] means that the first Gaussian will be divided
%        into 3 and the last into 2. The rest will remain the same.
%
% gmm  - Cell with the following format {m,b,W,n}, where there are K_p
%        Gaussians.
%
% A crude heuristic to replace a single VB Gaussian by a bunch of VB Gaussians.
% If there is only one Gaussian, then it should be the same as the
% original distribution.
%__________________________________________________________________________
% Copyright (C) 2018 Wellcome Centre for Human Neuroimaging

MU0 = gmm{1};
b0  = gmm{2};
W0  = gmm{3};
n0  = gmm{4};

Kb = numel(n0);
K  = numel(lkp);
C  = size(MU0,1);

A0 = bsxfun(@times,W0,reshape(n0,[1 1 Kb])); % E[Lambda]

m  = zeros(C,K);
b  = zeros(1,K);
W  = zeros(C,C,K);
n  = zeros(1,K);
mg_w = ones(1,K);

for k=1:Kb
    kk = sum(lkp==k);
    w  = 1./(1 + exp(-(kk - 1)*0.25)) - 0.5;
    mn = MU0(:,k);
    vr = inv(A0(:,:,k));

    mn = sqrtm(vr)*sort(randn(C,kk),2)*w + repmat(mn,[1 kk]);
    vr = vr*(1 - w);
    pr = inv(vr);
    W1 = (1/n0(k))*pr;

    m(:,lkp==k)   = mn;
    b(lkp==k)     = b0(k);
    W(:,:,lkp==k) = repmat(W1,[1 1 kk]);
    n(lkp==k)     = n0(k);

    mg_w(lkp==k) = 1/kk;
end

gmm{1} = m;
gmm{2} = b;
gmm{3} = W;
gmm{4} = n;
%==========================================================================

% =========================================================================
function [gmm,mg] = collapse_gmms(gmm,lkp)
% FORMAT gmm = spm_gmm_lib('extras', 'collapse', gmm, lkp)
%
% A crude heuristic to replace a bunch of Gaussian with a single Gaussian.
%__________________________________________________________________________
% Copyright (C) 2018 Wellcome Centre for Human Neuroimaging

MU0 = gmm{1};
b0  = gmm{2};
W0  = gmm{3};
n0  = gmm{4};

Kb = max(lkp);
K  = numel(lkp);
C  = size(MU0,1);

A0 = bsxfun(@times,W0,reshape(n0,[1 1 K]));
vr = zeros(size(A0));
for k=1:K
   vr(:,:,k) = inv(A0(:,:,k));
end

m  = zeros(C,Kb);
b  = zeros(1,Kb);
W  = zeros(C,C,Kb);
n  = zeros(1,Kb);
mg = ones(1,Kb);

for k=1:Kb
    kk = sum(lkp==k);

    mn  = mean(MU0(:,lkp == k),2);
    vr1 = mean(vr(:,:,lkp == k),3);

    pr = inv(vr1);
    W1 = (1/n0(k))*pr;

    m(:,k)   = mn;
    W(:,:,k) = W1;

    b(k) = mean(b0(lkp == k))/kk;
    n(k) = mean(n0(lkp == k))/kk;
end

gmm{1} = m;
gmm{2} = b;
gmm{3} = W;
gmm{4} = n;
%==========================================================================

%==========================================================================
% HELPER FUNCTIONS
%==========================================================================

%==========================================================================
function [m,b,W,n] = get_posteriors(cluster,s)
m = cluster{s}{1}{1};
b = cluster{s}{1}{2};
W = cluster{s}{2}{1};
n = cluster{s}{2}{2};
%==========================================================================

% === logdet =============================================================
function ld = logdet(A)
% A  - A postive-definite square matrix
% ld - Logarithm of determinant of A
%
% Log-determinant of a positive-definite matrix.
% Cholesky factorisation is used to compute a more stable log-determinant.
%__________________________________________________________________________
% Copyright (C) 2017 Wellcome Trust Centre for Neuroimaging

% John Ashburner
% $Id$

% Cholseki decomposition of A (A = C' * C, with C upper-triangular)
[C, p] = chol(A);

if p > 0
   % A should usually be positive definite, but check anyway.
   warning(['Attempting to compute log determinant of matrix ' ...
            'that is not positive definite (p=%d).'], p);
end

% Because C is triangular, |C| = prod(diag(C))
% Hence: log|C| = sum(log(diag(C)))
% And:   log|A| = log|C'*C| = log(|C|^2) = 2 * sum(log(diag(C)))
ld = 2 * sum(log(diag(C)));

% === loaddiag ===========================================================
function A = loaddiag(A)
% A  - A square matrix
%
% Load A's diagonal until it is well conditioned for inversion.
%__________________________________________________________________________
% Copyright (C) 2017 Wellcome Trust Centre for Neuroimaging

factor = 1e-7;
while rcond(A) < 1e-5
    A = A + factor * max([diag(A); eps]) * eye(size(A));
    factor = 10 * factor;
end

% === inv ================================================================
function A = inv_stable(A)
% A  - A positive-definite square matrix
% iA - Its inverse
%
% Stable inverse of a positive-definite matrix.
% Eigendecomposition is used to compute a more stable inverse.
%__________________________________________________________________________
% Copyright (C) 2017 Wellcome Trust Centre for Neuroimaging
% I           = eye(size(A,1));
% opts        = struct;
% opts.POSDEF = true;
% opts.SYM    = true;
% A           = linsolve(A, I, opts);
% A           = (A + A.')/2;
% [V,D] = eig(A);
% if any(diag(D) <= 0)
%     warning('[spm_gmm_lib::inv_stable] Matrix has negative eigenvalues')
%     D(D <= 0) = eps; % Threshold negative eigenvalues
% end
% D     = loaddiag(D);
% A     = real(V * (D \ V'));
% A     = (A+A.')/2;

A     = inv(A);
return

[V,D] = eig(A);
D     = max(diag(D), eps);
A     = real(V * bsxfun(@ldivide, D, V'));
A     = 0.5*(A + A.'); % Ensure symetric inverse

% === wishart_elogdet =====================================================
function out = wishart_elogdet(V, n, mode)
if nargin < 3 || mode(1) ~= 'n'
    K   = size(V, 1);
    out = DiGamma(0.5*n, K) + K*log(2) + logdet(V);
else
    out = wishart_elogdet(V/n, n);
end

% === wishart_kl ==========================================================
function kl = wishart_kl(varargin)
% FORMAT kl = wishart_kl(V1,      n1, V0,      n0)
% FORMAT kl = wishart_kl(lambda1, n1, lambda0, n0, 'normal')

% Check if we are in the reparameterised case
if nargin == 5
    kl = wishart_kl(varargin{1}/varargin{2}, varargin{2}, ...
                  varargin{3}/varargin{4}, varargin{4});
    return
end

% Usual KL
V1 = varargin{1};
n1 = varargin{2};
V0 = varargin{3};
n0 = varargin{4};
K  = size(V1, 1);
kl =   0.5*n0*(logdet(V0) - logdet(V1)) ...
     + 0.5*n1*(trace(V0\V1) - K) ...
     + 0.5*(n1 - n0)*DiGamma(0.5*n1, K) ...
     + LogGamma(0.5*n0, K) - LogGamma(0.5*n1, K);

 % === wishart_e ==========================================================
function out = wishart_e(V, n, mode)
if nargin < 3 || mode(1) ~= 'n'
    out = n*V;
else
    out = V;
end

% === LogGamma ============================================================
function lg = LogGamma(a, p)
if nargin < 2
    p = 1;
end
lg = (p*(p-1)/4)*log(pi);
for i=1:p
    lg = lg + gammaln(a + (1-p)/2);
end

% === DiGamma =============================================================
function dg = DiGamma(a, p, k)
if nargin < 3
    k = 0;
    if nargin < 2
        p = 1;
    end
end
dg = 0;
for i=1:p
    dg = dg + psi(k, a + (1-i)/2);
end

% === Gain ================================================================
function gain = get_gain(vals)
% FORMAT gain = get_gain(vals)
%
% vals - A vector of values
%
% gain - Computed gain
%
% Compute gain --- usually used to determine a stopping criteria when
% optimising
%__________________________________________________________________________
% Copyright (C) 2018 Wellcome Trust Centre for Neuroimaging
vals = vals(:);
gain = abs((vals(end - 1) - vals(end))/(max(vals(isfinite(vals))) - min(vals(isfinite(vals)))));
%==========================================================================

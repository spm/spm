function [L F df p] = spm_wilks(X, Y, C, persist)
% wilks -- Wilks' Lambda for multivariate General Linear Model
%
% Example:
%   [L F df p] = spm_wilks(X, Y, C)
%   [L F df p] = spm_wilks(X, Y, C, persist)
% For multivariate data in the rows of Y, compute the Wilks' Lambda, L,
% comparing the full model (Y = X*beta + error) to a reduced one defined by
% the null hyp. C'beta = 0, for contrast matrix C and multivariate beta,
% error is assumed to be multivariate Gaussian, i.i.d. over the rows of Y.
%
% If persist is true, multiple calls to wilks with the same design matrix
% and contrast, but different data Y, will be more efficiently computed.
%
% If requested, F (on df degrees of freedom) is the transformation of L to 
% Rao's F-statistic. This is exact in special cases (e.g. Hotelling's
% T-square test), and should be a good approximation in other cases.
% p is the p-value from this (possibly approximate) F statistic.
%
% NB: L is the reciprocal of Wilks' Lambda as more commonly defined;
% larger L indicate greater evidence to reject the null hypothesis. L is a
% monotonically increasing function of the likelihood ratio of the full to
% the reduced model, calculated by the ratio of the determinants of the
% matrices of residual sums of squares and products for the two models.
% F is a monotonically increasing function of L; for univariate data it
% reproduces the standard F-statistic.
% 
% References:
%   [1] Johnson, R. A. & Wichern, D. W.
%       Applied Multivariate Statistical Analysis.
%       5th Ed, 2002, Prentice Hall.
%   [2] Christensen, R
%       Advanced Linear Modeling
%       2nd Ed, 2001, Springer
%   [3] Anderson, T. W.
%       An introduction to multivariate statistical analysis
%       3rd Ed, 2003, Wiley
%% Ged Ridgway 2012
%
% $Id: spm_wilks.m 5548 2013-06-12 12:14:55Z gareth $



persistent N rnk Rx Rz
if nargin < 4 || ~persist
    % Clear any previous persistant definition, to re-evaluate...
    % (if persist is true and this is the first run, Rx will also be empty)
    Rx = [];
end

% Variables empty if first use of persistant, or if emptied due to ~persist
if isempty(Rx)
    [N(1) N(2)] = size(Y); N(3) = size(X, 2);
    if size(X, 1) ~= N(1);
        error('wilks:datadesrows', ...
            '#rows in data (%d) and des-mat. (%d) must match', ...
            N(1), size(X, 1));
    end
    if ~exist('C', 'var') || isempty(C);
        % Generate contrast for H0 of all columns being equal:
        C = diff(eye(N(3)))'; % (SPM C'*beta=0 convention)
    end
    if size(C, 1) ~= N(3)
        error('wilks:concols', ...
            '#cols in contrast (%d) and des. (%d) must match', ...
            size(C, 1), N(3));
    end

    In = eye(N(1));
    Px = X * pinv(X);
    rnk(1) = round(trace(Px)); % rank(X)
    Rx = In - Px;

    C0 = eye(N(3)) - C * pinv(C);
    Z = X * C0; % Reduced design matrix
    Pz = Z * pinv(Z);
    rnk(2) = round(trace(Pz)); % rank(Z)
    Rz = In - Pz;
end

XSS = Y' * Rx * Y; % full-model Sum of Squares (and products)
ZSS = Y' * Rz * Y; % reduced-model SS (and products)

L = det(ZSS) / det(XSS); % (reciprocal of) Wilks' Lambda statistic
% note that ZSS = XSS + HSS, where HSS is the Hypothesis SS (and products)

if nargout == 1,
    return
end
[F df] = transform(L, rnk(1), rnk(2), N(1), N(2));

if nargout > 3
    if exist('fcdf', 'file')
        p = 1 - fcdf(F, df(1), df(2));
    elseif exist('spm_Fcdf', 'file')
        p = 1 - spm_Fcdf(F, df);
    else
        warning('wilks:fcdf', 'No stats toolbox or SPM, returning p = []');
        p = [];
    end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [F df] = transform(L, Xr, Zr, Nn, Nm)
% Determine appropriate transformation of L.
% Rao's F approximation from [2,3], exact if min(Nm, Xr - Zr) is 1 or 2.
d = Xr - Zr;
% s = Nm * d / 2 + 1; % mistake in [2],
s = Nm * d / 2 - 1; % as given in [3]
f = (Nn - Xr) + d - (d + Nm + 1)/2;
if min(Nm, d) == 1
    t = 1;
else
    t = sqrt((Nm^2*d^2 - 4) / (Nm^2 + d^2 - 5));
end
df = [Nm*d f*t-s];
F = (L^(1/t) - 1) * df(2) / df(1);

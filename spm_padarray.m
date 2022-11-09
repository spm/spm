function Y = spm_padarray(X, padsize, method, direction)
% FORMAT Y = spm_padarray(X, padsize, [method], [direction])
% X         - numeric array
% padsize   - padding size along each dimension of the array (>= 0)
% method    - 'circular', 'replicate', 'symmetric' or a value [0]
% direction - 'pre'/'post'/['both']
%
% Note that:
% * 'circular'  corresponds to the boundary condition of an FFT
% * 'symmetric' corresponds to the boundary condition of a DCT-II
%
% If padsize < 0, it is set to 0 instead.
%__________________________________________________________________________

% Yael Balbastre
% Copyright (C) 2019-2022 Wellcome Centre for Human Neuroimaging


% Possible extensions:
% * method and direction could have a different value per dimension
% * other methods: antisymmetric (DST/Dirichlet)
%                  symmetry wrt voxel center rather than voxel side (DCT-I)


if nargin < 3, method = 0;         end
if nargin < 4, direction = 'both'; end

% Compute output dimensions
%--------------------------------------------------------------------------
dim     = size(X);
dim     = [dim ones(1,max(numel(padsize)-numel(dim),0))];
padsize = reshape(padsize, 1, []);
padsize = [padsize zeros(1,numel(dim)-numel(padsize))];
padsize(padsize < 0 | ~isfinite(padsize)) = 0;
newdim  = dim + padsize * (strcmpi(direction, 'both') + 1);
if strcmpi(direction,'post')
    padsize = zeros(size(padsize));
end

% Pad array
%--------------------------------------------------------------------------
S = substruct('()', repmat({':'}, [1 numel(dim)]));

if ~ischar(method)
    % Allocate output volume
    Y = repmat(cast(method, class(X)), newdim);
    
    % Copy input volume
    for d=1:numel(dim)
        S.subs{d} = padsize(d) + (1:dim(d));
    end
    Y = subsasgn(Y,S,X);
else
    % Padding method
    switch lower(method)
        case 'circular'
            fillidx = @circidx;
        case 'replicate'
            fillidx = @repidx;
        case 'symmetric'
            fillidx = @symidx;
        otherwise
            error('Unknown method %s.', method);
    end
    
    % Sample input volume
    for d=1:numel(dim)
        S.subs{d} = fillidx((1:newdim(d)) - padsize(d), dim(d));
    end
    Y = subsref(X,S);
end


%==========================================================================
% Padding method functions
%==========================================================================
function i = circidx(i,n)

i          = i - 1;
nonneg     = (i >= 0);
i(nonneg)  = mod(i(nonneg), n);
i(~nonneg) = mod(n + mod(i(~nonneg), n), n);
i          = i + 1;

function i = repidx(i,n)

pre        = (i <= 0);
post       = (i > n);
i(pre)     = 1;
i(post)    = n;

function i = symidx(i,n)

i          = i - 1;
n2         = n*2;
pre        = (i < 0);
i(pre)     = n2 - 1 - mod(-i(pre)-1, n2);
i(~pre)    = mod(i(~pre), n2);
post       = (i >= n);
i(post)    = n2 - i(post) - 1;
i          = i + 1;

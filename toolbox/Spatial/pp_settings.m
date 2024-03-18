function sett = pp_settings(deg,bnd,ext)
% Settings for push/pull
% FORMAT sett = pp_settings(deg,bnd,ext)
% deg - interpolation degree in each dimension (3 elements)
%       0 - nearest neighbour
%       1 - trilinear
%       2 - cubic B-spline
%       3 - 3rd degree B-spline
%       4 - 4th degree B-spline
% bnd - boundary conditions in each dimension (3 elements)
%       0 - circulant
%       1 - reflected
%       2 - reflected negative
% ext - extrapolation flag 0/1
%__________________________________________________________________________

% John Ashburner
% Copyright (C) 2023 Wellcome Centre for Human Neuroimaging


if true
    usize = @uint64;
else
    usize = @uint32;
end

% Interpolation degree
if nargin<1
    deg = usize([1 1 1]);
else
    deg = usize(deg);
    if isscalar(deg)
        deg = repmat(deg,[1 3]);
    else
        if numel(deg)~=3, error('"deg" must be scalar or have three elements.'); end
    end
    if any(deg>4) || any(deg<0), error('"deg" out of range.'); end
end


% Boundary conditions
if nargin<2
    bnd = int32([1 1 1]);
else
    bnd = int32(bnd);
    if isscalar(bnd)
        bnd = repmat(bnd,[1 3]);
    else
        if numel(bnd)~=3, error('"bnd" must be scalar or have three elements.'); end
    end
    if any(bnd>2) || any(bnd<0), error('"bnd" out of range.'); end
end


% Extrapolate flag
if nargin<3
    ext = int32(1);
else
    if ~isscalar(ext), error('"ext" must be scalar.'); end
    if ext<0 || ext>1, error('"ext" out of range.');   end
    ext = int32(ext);
end

sett = struct('deg',deg, 'bnd',bnd, 'ext',ext);

function [Y, xY] = spm_summarise(V,xY,fhandle,keepNaNs)
% Summarise data within a Region of Interest
% FORMAT [Y, xY] = spm_summarise(V,xY,fhandle)
% V       - [1 x n] vector of mapped image volumes to read (from spm_vol)
%           Or a char array of filenames
% xY      - VOI structure (from spm_ROI)
%           Or a VOI_*.mat (from spm_regions) or a mask image filename
%           Or the keyword 'all' to summarise all voxels in the images
%           Or a [3 x m] matrix of voxel coordinates {mm}
% fhandle - function handle to be applied on image data within VOI
%           Must transform a [1 x m] array into a [1 x p] array
%           Default is Identity (returns raw data, vectorised into rows).
%           Can also use keyword 'litres' to compute the total volume,
%           within the region of interest, for a tissue segment image.
%
% Y       - [n x p] data summary
% xY      - (updated) VOI structure
%__________________________________________________________________________
%
% Example:
% spm_summarise('beta_0001.nii',...
%               struct('def','sphere', 'spec',8, 'xyz',[10 20 30]'),...
%               @mean)
%__________________________________________________________________________
% Copyright (C) 2010-2015 Wellcome Trust Centre for Neuroimaging

% Guillaume Flandin, Ged Ridgway
% $Id: spm_summarise.m 6490 2015-06-26 11:51:46Z guillaume $

%-Argument checks
%--------------------------------------------------------------------------
if nargin < 1 || isempty(V)
    [V, sts] = spm_select([1 Inf], 'image', 'Specify Images');
    if ~sts, error('Must select 1 or more images'), end
end
if iscellstr(V), V = char(V); end
if ischar(V), V = spm_data_hdr_read(V); end
spm_check_orientations(V);

if nargin < 2 || isempty(xY), xY = struct; end
if ischar(xY)
    if strcmpi(xY, 'all')
        xY = struct('def', 'all');
    elseif any(regexpi(xY, '\.mat$'))
        try
            load(xY,'xY'); % VOI_*.mat file. Warns if .mat has no xY ...
            xY = rmfield(xY,'XYZmm'); % ...  error if .mat has no xY
        catch
            xY = struct; % GUI specification in spm_ROI 
        end
    else % assume mask image filename
        xY = struct('def','mask', 'spec',xY);
    end
elseif isnumeric(xY) && any(size(xY, 1) == [3 4])
    xY = struct('XYZmm', xY(1:3, :));
elseif isstruct(xY) && isfield(xY,'fname')
    xY = struct('def','mask', 'spec',xY);
elseif ~isstruct(xY)
    error('Incorrect xY specified')
end
if ~isfield(xY,'XYZmm'), [xY, xY.XYZmm] = spm_ROI(xY,V(1)); end

if nargin < 3 || isempty(fhandle), fhandle = @(x) x; end
if ischar(fhandle) && strcmp(fhandle, 'litres')
    vsz     = abs(det(V(1).mat));  % voxel size in mm^3
    fhandle = @(x) sum(x) * vsz / 1e6;
end
if ischar(fhandle), fhandle = str2func(fhandle); end

% Undocumented option in case anyone wants to keep (e.g. to check for) NaNs
if nargin < 4, keepNaNs = false; end
if keepNaNs
    dropNaNs = @(x) x;
else
    dropNaNs = @(x) x(~isnan(x));
end

%-Summarise data
%--------------------------------------------------------------------------
XYZ = round(V(1).mat \ [xY.XYZmm; ones(1, size(xY.XYZmm, 2))]);

% Run on first volume to determine p, and transpose if column vector
Y = fhandle(dropNaNs(spm_data_read(V(1), 'xyz', XYZ)));
if ndims(Y) > 2
    error('Function must return a [1 x p] array')
elseif size(Y, 1) ~= 1
    if size(Y, 2) == 1
        Y = Y';
    else
        error('Function returned a [%d x %d] array instead of [1 x p]', ...
            size(Y, 1), size(Y, 2))
    end
end

% Preallocate space and then run on remaining volumes
Y(2:numel(V), :) = 0;
for i = 2:numel(V)
    Y(i, :) = fhandle(dropNaNs(spm_data_read(V(i), 'xyz', XYZ)));
end

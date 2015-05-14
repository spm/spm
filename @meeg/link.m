function this = link(this, fnamedat, dtype, slope, offset)
% Links the object to data file (only if exists)
% FORMAT this = link(this) 
%   Will try to find the datafile based on fname and path     
% FORMAT this = link(this, fnamedat) 
%   Will find the datafile using the provided name and path
% FORMAT this = link(this, fnamedat, dtype, slope, offset)
%   Additional parameters for non-float data files
% _________________________________________________________________________
% Copyright (C) 2011 Wellcome Trust Centre for Neuroimaging

% Vladimir Litvak
% $Id: link.m 6437 2015-05-14 12:27:21Z vladimir $

if isempty(this)
   error('All header dimensions should be >0');
end

if nargin == 1 || isempty(fnamedat)
    [p, f] = fileparts(fullfile(this));
    fnamedat = fullfile(p, [f '.dat']);
else
   [p, f, x] = fileparts(fnamedat);
   if isempty(p)
       p = path(this);
   end
   if isempty(x)
       x = '.dat';
   end
   fnamedat = fullfile(p, [f x]); 
end

% This is to re-use existing settings for non-floats. Use unlink to clear
% these settings.
if (nargin < 3) && isa(this.data, 'file_array') && ~isequal(lower(this.data.dtype), 'float32-le')
    dtype  = this.data.dtype;
    offset = this.data.offset;
    slope  = this.data.scl_slope;
else
    if nargin < 3, dtype  = 'float32-le';  end
    if nargin < 4, slope  = 1;             end
    if nargin < 5, offset = 0;             end
end

if ~exist(fnamedat, 'file')
    error('Data file not found');
end

% Size determination here should worked for unlinked objects and be insensitive
% to online montages.
if ~strncmpi(transformtype(this), 'TF', 2)
    siz = [length(this.channels), nsamples(this), ntrials(this)];
else
    siz = [length(this.channels), nfrequencies(this), nsamples(this), ntrials(this)];
end

this.data = file_array(fnamedat, siz, dtype, offset, slope);

siz = num2cell(size(this));

try    
    this.data(siz{:});
catch
    error('Dimensions mismatch. Could not link to the data file');
end

this = check(this);

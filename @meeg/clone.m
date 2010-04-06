function new = clone(this, fnamedat, dim, reset)
% Creates a copy of the object with a new, empty data file,
% possibly changing dimensions
% FORMAT new = clone(this, fnamedat, dim, reset)
% reset - 0 (default) do not reset channel or trial info unless dimensions
%          change, 1 - reset channels only, 2 - trials only, 3 both
% Note that when fnamedat comes with a path, the cloned meeg object uses
% it. Otherwise, its path is by definition that of the meeg object to be
% cloned.
% _________________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% Stefan Kiebel, Vladimir Litvak
% $Id: clone.m 3807 2010-04-06 19:29:08Z vladimir $

if nargin < 4
    reset = 0;
end

if nargin < 3
    if ~strcmp(transformtype(this), 'TF')
        dim = [nchannels(this), nsamples(this), ntrials(this)];
    else
        dim = [nchannels(this), nfrequencies(this), nsamples(this), ntrials(this)];
    end
end

new = this;

% check file path first
[pth,fname,ext] = fileparts(fnamedat);
if isempty(pth)
    pth = this.path;
end
newFileName = [fullfile(pth,fname),ext];
% copy the file_array
d = this.data.y; % 
d.fname = newFileName;
d.dim = dim; 

% physically initialise file
if length(dim) == 3
    d(end, end, end) = 0;
    nsampl = dim(2);
    ntrial = dim(3);
elseif length(dim) == 4
    d(end, end, end, end) = 0;
    nsampl = dim(3);
    ntrial = dim(4);
    
    if ~strncmpi(transformtype(new), 'TF',2)
        new = transformtype(new, 'TF');
    end        
else
   error('Dimensions different from 3 or 4 are not supported.');
end

% link into new meeg object
new.data.y = d;

% change filenames
new.data.fnamedat = [fname,'.dat'];
new.fname = [fname,'.mat'];
new.path = pth;

% ensure consistency 
if (dim(1) ~= nchannels(this)) || ismember(reset, [1 3])
    new.channels = [];
    for i = 1:dim(1)
        new.channels(i).label = ['Ch' num2str(i)];
    end
end

if ntrial ~= ntrials(this) || ismember(reset, [2 3])
    new.trials = repmat(struct('label', 'Undefined'), 1, ntrial);
end
    
if (nsampl ~= nsamples(this))
    new.Nsamples = nsampl;
end

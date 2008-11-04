function new = clone(this, fnamedat, dim)
% Creates a copy of the object with a new, empty data file,
% possibly changing dimensions
% FORMAT new = clone(this, fnamedat, dim)
% Note that when fnamedat comes with a path, the cloned meeg object uses
% it. Otherwise, its path is by definition that of the meeg object to be
% cloned.
% _________________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% Stefan Kiebel, Vladimir Litvak
% $Id: clone.m 2436 2008-11-04 10:46:27Z stefan $

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
if isempty(path)
    pth = this.path;
end
newFileName = [fullfile(pth,fname),ext];
% initialise new file_array
d = file_array(newFileName, dim, dtype(this));

% physically initialise file
if length(dim) == 3
    d(end, end, end) = 0;
    nsampl = dim(2);
    ntrial = dim(3);
elseif length(dim) == 4
    d(end, end, end, end) = 0;
    nsampl = dim(3);
    ntrial = dim(4);
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
if dim(1) ~= nchannels(this)
    new.channels = [];
    for i = 1:dim(1)
        new.channels(i).label = ['Ch' num2str(i)];
    end
end

if ntrial ~= ntrials(this)
    new.trials = repmat(struct('label', 'Undefined'), 1, ntrial);
end
    
if (nsampl ~= nsamples(this))
    new.Nsamples = nsampl;
end

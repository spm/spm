function new = clone(this, fnamedat, dim)
% Creates a copy of the object with a new, empty data file,
% possibly changing dimensions
% FORMAT new = clone(this, fnamedat, dim)
% _________________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% Stefan Kiebel, Vladimir Litvak
% $Id: clone.m 1263 2008-03-28 10:30:10Z stefan $

if nargin < 3
    dim = [nchannels(this), nsamples(this), ntrials(this)];
end

new = this;

% initialise new file_array thisect
d = file_array(fullfile(this.path,fnamedat), dim, dtype(this));

% physically initialise file
if length(dim) == 3
    d(end, end, end) = 0;
    nsampl = dim(2);
    ntrial = dim(3);
elseif length(dim) == 2
    d(end, end, end, end) = 0;
    nsampl = dim(3);
    ntrial = dim(4);
else
   error('Dimensions different from 3 or 4 are not supported.');
end

% link into new meeg object
new.data.y = d;

% change filenames
new.fname = [spm_str_manip(fnamedat, 'rt') '.mat'];
new.data.fnamedat = fnamedat;

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

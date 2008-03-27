function new = clone(this, fnamedat, dim)
% Creates a copy of the object with a new, empty data file,
% possibly changing dimensions
% FORMAT new = clone(this, fnamedat, dim)
% _________________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% Stefan Kiebel, Vladimir Litvak
% $Id: clone.m 1254 2008-03-27 18:41:42Z vladimir $

if nargin < 3
    dim = [nchannels(this), nsamples(this), ntrials(this)];
end

nchan = dim(1);
nsampl = dim(2);
ntrial = dim(3);

new = this;

% initialise new file_array thisect
d = file_array(fullfile(this.path,fnamedat),...
    [nchan nsampl ntrial], dtype(this));

% physically initialise file
d(end,end,end) = 0;

% link into meeg thisect
new.data.y = d;

% change filenames
new.fname = [spm_str_manip(fnamedat, 'r') '.mat'];
new.data.fnamedat = fnamedat;

% ensure consistency 
if nchan ~= nchannels(this)
    new.channels =[];
    for i = 1:nchan
        new.channels(i).label = ['Ch' num2str(i)];
    end
end

if ntrial ~= ntrials(this)
    new.trials = repmat(struct('label', 'Undefined'), 1, ntrial);
end
    
if (nsampl ~= nsamples(this))
    new.Nsamples = nsampl;
end
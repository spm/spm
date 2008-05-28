function raw = ftraw(this, memmap)
% Method for converting meeg object to Fieldtrip raw struct
% FORMAT raw = ftraw(this, memmap)
%        memmap - 1 (default) memory map the data with file_array
%                 0 load the data into memory
% _______________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% Vladimir Litvak
% $Id: ftraw.m 1742 2008-05-28 11:58:04Z vladimir $

if nargin < 2
    memmap = 1;
end

nchans  = nchannels(this);
ntime   = nsamples(this);
ntrl =    ntrials(this);

raw = [];
raw.fsample = fsample(this);
raw.label   = chanlabels(this)';

wsize=wordsize(this.data.datatype);

for i=1:ntrl
    if memmap
        offset = (i-1)*nchans*ntime*wsize;
        trialdim = [nchans ntime];
        raw.trial{i} = file_array(fullfile(this.path, this.data.fnamedat),trialdim,this.data.datatype,offset,1,0,'ro');
    else
        raw.trial{i} = this.data.y(:, :, i);
    end
end

raw.time = repmat({time(this)}, 1, ntrl);

function s = wordsize(datatype)

datatype = strtok(datatype, '-');

switch datatype
    case 'int8'
        s = 1;
    case 'int16'
        s = 2;
    case 'int32'
        s = 4;
    case 'uint8'
        s = 1;
    case 'uint16'
        s = 2;
    case 'uint32'
        s = 4;
    case 'float'
        s = 4;
    case 'double'
        s = 8;
    case 'float32'
        s = 4;
    case 'float64'
        s = 8;
    otherwise
        error('unrecognized datatype');
        % FIXME, add support for le and be versions
end
function raw = ftraw(this)
% Method for converting meeg object to Fieldtrip raw struct
% FORMAT raw = ftraw(this)
% _______________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% Vladimir Litvak
% $Id: ftraw.m 1390 2008-04-14 16:08:09Z vladimir $


nchans  = nchannels(this);
ntime   = nsamples(this);
ntrl =    ntrials(this);

raw = [];
raw.fsample = fsample(this);
raw.label   = chanlabels(this)';

wsize=wordsize(this.data.datatype);

for i=1:ntrl
    offset = (i-1)*nchans*ntime*wsize;
    trialdim = [nchans ntime];
    raw.trial{i} = file_array(fullfile(this.path, this.data.fnamedat),trialdim,this.data.datatype,offset,1,0,'ro');
end

raw.time = repmat({time(this)}, 1, ntrl);

function s = wordsize(datatype)
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
function spm_copy(source, dest, varargin)
% Copy file(s)
% FORMAT spm_copy(source, dest [,opts])
%__________________________________________________________________________
% Copyright (C) 2017 Wellcome Trust Centre for Neuroimaging

% Guillaume Flandin
% $Id: spm_copy.m 7121 2017-06-21 16:35:40Z guillaume $


%-Source and destination
%--------------------------------------------------------------------------
source = cellstr(source);
if nargin < 2
    dest = pwd;
end
dest = cellstr(dest);
if numel(source) == 1
    source = repmat(source,numel(dest),1);
elseif numel(dest) == 1
    dest = repmat(dest,numel(source),1);
elseif numel(source) ~= numel(dest)
    error('Number of elements in source and dest must be one or equal.');
end

%-Options (struct array or key/value pairs)
%--------------------------------------------------------------------------
opts = struct('gzip', false, 'gunzip', false, 'nifti', false, 'mode', {{}});
if nargin > 2
    if isstruct(varargin{1})
        opt = varargin{1};
    else
        opt = struct;
        for i=1:2:numel(varargin)
            opt.(varargin{i}) = varargin{i+1};
        end
    end
else
    opt = struct([]);
end
fn = fieldnames(opt);
for i=1:numel(fn)
    if ~isfield(opts,lower(fn{i}))
        warning('Unknown option "%s".',fn{i});
    end
    opts.(lower(fn{i})) = opt.(fn{i});
end

%-Actual copy
%--------------------------------------------------------------------------
for i=1:numel(source)
    protocol = source{i}(1:find(source{i}==':',1)-1);
    if ismember(protocol,{'file','http','https','ftp'})
        urlwrite(source{i}, dest{i}); % dest{i} has to be a filename...
    else
        sts = copyfile(source{i}, dest{i}, opts.mode{:});
    end
    if opts.nifti
        if strcmp(spm_file(source{i},'ext'),'img')
            s = copyfile(spm_file(source{i},'ext','hdr'), dest{i}, opts.mode{:});
            s = copyfile(spm_file(source{i},'ext','mat'), dest{i}, opts.mode{:});
        elseif strcmp(spm_file(source{i},'ext'),'hdr')
            s = copyfile(spm_file(source{i},'ext','img'), dest{i}, opts.mode{:});
            s = copyfile(spm_file(source{i},'ext','mat'), dest{i}, opts.mode{:});
        elseif strcmp(spm_file(source{i},'ext'),'nii')
            s = copyfile(spm_file(source{i},'ext','mat'), dest{i}, opts.mode{:});
        end
    end
    if opts.gzip && ~strcmp(spm_file(source{i},'ext'),'gz')
        gzip(spm_file(source{i},'path',dest{i}));
        spm_unlink(spm_file(source{i},'path',dest{i}));
    end
    if opts.gunzip && strcmp(spm_file(source{i},'ext'),'gz')
        gunzip(spm_file(source{i},'path',dest{i}));
        spm_unlink(spm_file(source{i},'path',dest{i}));
    end
end

function spm_copy(source, dest, varargin)
% Copy file(s)
% FORMAT spm_copy(source, dest [,opts])
% source  - pathnames of files or directories to be copied
%           character vector or cellstr
% dest    - pathnames of destination files or directories [default: pwd]
%           character vector or cellstr
% opts    - structure or list of name/value pairs of optional parameters:
%             gzip: compress uncompressed copied files at destination
%             gunzip: uncompress compressed copied files at destination
%             nifti: also copy sidecar .hdr/.img/.mat/.json if present
%             gifti: also copy sidecar .dat file if present
%             mode: copy mode (see copyfile's help)
%__________________________________________________________________________

% Guillaume Flandin
% Copyright (C) 2017-2023 Wellcome Centre for Human Neuroimaging


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
opts = struct('gzip', false, 'gunzip', false, ...
    'nifti', false, 'gifti', false, ...
    'mode', {{}});
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
    ext = file_ext(source{i});
    if opts.nifti
        if strcmp(ext,'img') || strcmp(ext,'img.gz')
            sts = copyfile(file_ext(source{i},'hdr'), dest{i}, opts.mode{:});
            sts = copyfile(file_ext(source{i},'mat'), dest{i}, opts.mode{:});
            sts = copyfile(file_ext(source{i},'json'), dest{i}, opts.mode{:});
        elseif strcmp(ext,'hdr') || strcmp(ext,'hdr.gz')
            sts = copyfile(file_ext(source{i},'img'), dest{i}, opts.mode{:});
            sts = copyfile(file_ext(source{i},'mat'), dest{i}, opts.mode{:});
            sts = copyfile(file_ext(source{i},'json'), dest{i}, opts.mode{:});
        elseif strcmp(ext,'nii') || strcmp(ext,'nii.gz')
            sts = copyfile(file_ext(source{i},'mat'), dest{i}, opts.mode{:});
            sts = copyfile(file_ext(source{i},'json'), dest{i}, opts.mode{:});
        end
    end
    if opts.gifti
        if strcmp(ext,'gii') || strcmp(ext,'gii.gz')
            sts = copyfile(spm_file(source{i},'ext','dat'), dest{i}, opts.mode{:});
            % edit ExternalFileName in .gii file if full path or renaming
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


%==========================================================================
function ext = file_ext(file,ext)
% Equivalent to spm_file(file,'ext') with special case for .gz extension
if nargin == 1
    if strcmp(spm_file(file,'ext'),'gz')
        ext = [spm_file(file(1:end-3),'ext') '.gz'];
    else
        ext = spm_file(file,'ext');
    end
else
    if strcmp(spm_file(file,'ext'),'gz')
        ext = spm_file(file(1:end-3),'ext',ext);
    else
        ext = spm_file(file,'ext',ext);
    end
end

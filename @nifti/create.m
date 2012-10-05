function create(obj,varargin)
% Create a NIFTI-1 file
% FORMAT create(obj)
% This writes out the header information for the nifti object
%
% create(obj,wrt)
% This also writes out an empty image volume if wrt==1
%__________________________________________________________________________
% Copyright (C) 2005-2012 Wellcome Trust Centre for Neuroimaging

%
% $Id: create.m 4986 2012-10-05 17:35:09Z guillaume $

for i=1:numel(obj)
    create_each(obj(i),varargin{:});
end

%==========================================================================
%-function create_each(obj, wrt)
%==========================================================================
function create_each(obj, wrt)
if ~isa(obj.dat,'file_array')
    error('Data must be a file_array');
end

fname = obj.dat.fname;
if isempty(fname)
    error('No filename to write to.');
end

dt = obj.dat.dtype;
sts = write_hdr_raw(fname,obj.hdr,dt(end-1)=='B');
if ~sts
    error('Unable to write header for "%s".',fname);
end

write_extras(fname,obj.extras);

if nargin>1 && any(wrt==1)
    % Create an empty image file if necessary
    d   = findindict(obj.hdr.datatype, 'dtype');
    dim = double(obj.hdr.dim(2:end));
    dim((double(obj.hdr.dim(1))+1):end) = 1;
    nbytes = ceil(d.size*d.nelem*prod(dim(1:2)))*prod(dim(3:end))+double(obj.hdr.vox_offset);
    
    [pth,nam] = fileparts(obj.dat.fname);
    if any(strcmp(obj.hdr.magic(1:3),{'n+1','n+2'}))
        ext = '.nii';
    else
        ext = '.img';
    end
    iname = fullfile(pth,[nam ext]);
    fp    = fopen(iname,'a+');
    if fp==-1
        error(['Unable to create image for "' fname '".']);
    end

    fseek(fp,0,'eof');
    pos = ftell(fp);
    if pos<nbytes
        bs      = 2048; % Buffer-size
        nbytes  = nbytes - pos;
        buf     = uint8(0);
        buf(bs) = 0;
        while(nbytes>0)
            if nbytes<bs, buf = buf(1:nbytes); end
            nw = fwrite(fp,buf,'uint8');
            if nw<min(bs,nbytes)
                fclose(fp);
                error(['Problem while creating image for "' fname '".']);
            end
            nbytes = nbytes - nw;
        end
    end
    fclose(fp);
end

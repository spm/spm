function create(obj)
% Create a NIFTI-1 file
% This writes out the header information, and an empty image
% volume if necesary.
% _______________________________________________________________________
% %W% John Ashburner %E%

if ~isa(obj.dat,'file_array'),
    error('Data must be a file-array');
end;
fname = obj.dat.fname;
if isempty(fname),
    error('No filename to write to.');
end;
dt = obj.dat.dtype;
ok = write_hdr_raw(fname,obj.hdr,dt(end-1)=='B');
if ~ok,
    error(['Unable to write header for "' fname '".']);
end;

write_extras(fname,obj.extras);

% Create an empty image file if necessary
d   = findindict(obj.hdr.datatype, 'dtype');
dim = double(obj.hdr.dim(2:end));
dim((double(obj.hdr.dim(1))+1):end) = 1;
nbytes = ceil(d.size*d.nelem*prod(dim(1:2)))*prod(dim(3:end))+double(obj.hdr.vox_offset);
[pth,nam,ext] = fileparts(obj.dat.fname);
if any(strcmp(deblank(obj.hdr.magic),{'n+1','nx1'})),
    ext = '.nii';
else
    ext = '.img';
end;
iname = fullfile(pth,[nam ext]);
fp    = fopen(iname,'a+');
if fp==-1,
    error(['Unable to create image for "' fname '".']);
end;

if 0,
    fseek(fp,0,'eof');
    pos = ftell(fp);
    if pos<nbytes,
        bs      = 2048; % Buffer-size
        nbytes  = nbytes - pos;
        buf     = uint8(0);
        buf(bs) = 0;
        while(nbytes>0)
            if nbytes<bs, buf = buf(1:nbytes); end;
            nw = fwrite(fp,buf,'uint8');
            if nw<min(bs,nbytes),
                fclose(fp);
                error(['Problem while creating image for "' fname '".']);
            end;
            nbytes = nbytes - nw;
        end;
    end;
end;

fclose(fp);
return;


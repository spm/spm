function sts = write_hdr_raw(fname,hdr,be)
% Write a NIFTI-1 header
% FORMAT sts = write_hdr_raw(fname,hdr,be)
% fname      - filename of image
% hdr        - a structure containing hdr info
% be         - whether big-endian or not [Default: native]
%
% sts        - status (1=good, 0=bad)
%__________________________________________________________________________
% Copyright (C) 2005-2017 Wellcome Trust Centre for Neuroimaging

%
% $Id: write_hdr_raw.m 7370 2018-07-09 10:44:51Z guillaume $


[pth,nam] = fileparts(fname);
if isempty(pth), pth = pwd; end

nifti1_bytes = 348;
nifti2_bytes = 540;

if isfield(hdr,'magic')
    switch hdr.magic(1:3)
        case {'ni1'}
            org   = niftistruc('nifti1');
            hname = fullfile(pth,[nam '.hdr']);
            bytes = zeros(nifti1_bytes,1,'uint8');
        case {'ni2'}
            org   = niftistruc('nifti2');
            hname = fullfile(pth,[nam '.hdr']);
            bytes = zeros(nifti2_bytes,1,'uint8');
        case {'n+1'}
            org   = niftistruc('nifti1');
            hname = fullfile(pth,[nam '.nii']);
            bytes = zeros(nifti1_bytes,1,'uint8');
        case {'n+2'}
            org   = niftistruc('nifti2');
            hname = fullfile(pth,[nam '.nii']);
            bytes = zeros(nifti2_bytes,1,'uint8');
        otherwise
            error('Bad header.');
    end
else
    org   = mayostruc;
    hname = fullfile(pth,[nam '.hdr']);
    bytes = zeros(nifti1_bytes,1,'uint8');
end

if nargin >= 3
    if be, mach = 'ieee-be';
    else   mach = 'ieee-le';
    end
else       mach = 'native';
end

sts = true;
try
    is_file = spm_existfile(hname);
catch
    is_file = exist(hname,'file') > 0;
end
if is_file
    [fp,msg] = fopen(hname,'r+',mach);
else
    [fp,msg] = fopen(hname,'w+',mach);
end
if fp == -1
    sts = false;
    fprintf('Error: %s\n',msg);
end

pos = 0;

if sts
    for i=1:length(org)
        if isfield(hdr,org(i).label)
            dat = hdr.(org(i).label);
            if length(dat) ~= org(i).len
                if length(dat)< org(i).len
                    if ischar(dat), z = char(0); else z = 0; end
                    dat = [dat(:) ; repmat(z,org(i).len-length(dat),1)];
                else
                    dat = dat(1:org(i).len);
                end
            end
        else
            dat = org(i).def;
        end

        if be
            d = typecast(swapbytes(cast(dat,org(i).dtype.prec)),'uint8');
        else
            d = typecast(cast(dat,org(i).dtype.prec),'uint8');
        end
        len = numel(d);
        bytes((1:len) + pos) = d;
        pos = pos + len;
    end
    sts = sts && (fwrite(fp,bytes,'uint8') == numel(bytes));
    fclose(fp);
end

if ~sts
     fprintf('There was a problem writing to the header of\n');
     fprintf('  "%s"\n', fname);
end


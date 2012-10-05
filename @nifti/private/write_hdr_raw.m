function sts = write_hdr_raw(fname,hdr,be)
% Write a NIFTI-1 header
% FORMAT ok = write_hdr_raw(fname,hdr,be)
% fname     - filename of image
% hdr       - a structure containing hdr info
% be        - whether big-endian or not [Default: native]
%
% sts       - status (1=good, 0=bad)
%__________________________________________________________________________
% Copyright (C) 2005-2012 Wellcome Trust Centre for Neuroimaging

%
% $Id: write_hdr_raw.m 4986 2012-10-05 17:35:09Z guillaume $


[pth,nam,ext] = fileparts(fname);
if isempty(pth), pth = pwd; end

if isfield(hdr,'magic')
    switch hdr.magic(1:3)
        case {'ni1'}
            org = niftistruc('nifti1');
            hname = fullfile(pth,[nam '.hdr']);
        case {'ni2'}
            org = niftistruc('nifti2');
            hname = fullfile(pth,[nam '.hdr']);
        case {'n+1'}
            org = niftistruc('nifti1');
            hname = fullfile(pth,[nam '.nii']);
        case {'n+2'}
            org = niftistruc('nifti2');
            hname = fullfile(pth,[nam '.nii']);
        otherwise
            error('Bad header.');
    end
else
    org   = mayostruc;
    hname = fullfile(pth,[nam '.hdr']);
end

if nargin >=3
    if be, mach = 'ieee-be';
    else   mach = 'ieee-le';
    end
else       mach = 'native';
end

sts = true;
if spm_existfile(hname)
    fp = fopen(hname,'r+',mach);
else
    fp = fopen(hname,'w+',mach);
end
if fp == -1
    sts = false;
    return;
end

for i=1:length(org)
    if isfield(hdr,org(i).label)
        dat = hdr.(org(i).label);
        if length(dat) ~= org(i).len
            if length(dat)< org(i).len
                dat = [dat(:) ; zeros(org(i).len-length(dat),1)];
            else
                dat = dat(1:org(i).len);
            end
        end
    else
        dat = org(i).def;
    end
    % fprintf('%s=\n',org(i).label)
    % disp(dat)
    len = fwrite(fp,dat,org(i).dtype.prec);
    if len ~= org(i).len
        sts = false;
    end
end

fclose(fp);
if ~sts
     fprintf('There was a problem writing to the header of\n');
     fprintf('"%s"\n', fname);
end

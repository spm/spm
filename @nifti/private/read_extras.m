function extras = read_extras(fname)
% Read extra bits of information
%__________________________________________________________________________
% Copyright (C) 2005-2017 Wellcome Trust Centre for Neuroimaging

%
% $Id: read_extras.m 7147 2017-08-03 14:07:01Z spm $


extras = struct;
[pth,nam,ext] = fileparts(fname);
switch ext
    case {'.hdr','.img','.nii'}
        mname = fullfile(pth,[nam '.mat']);
    case {'.HDR','.IMG','.NII'}
        mname = fullfile(pth,[nam '.MAT']);
    otherwise
        mname = fullfile(pth,[nam '.mat']);
end

if spm_existfile(mname)
    try
        extras = load(mname);
    catch
        warning('Can not load "%s" as a binary MAT file.', mname);
    end
end

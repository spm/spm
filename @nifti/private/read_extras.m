function extras = read_extras(fname)
% Read extra bits of information
%_______________________________________________________________________
% Copyright (C) 2005 Wellcome Department of Imaging Neuroscience

%
% $Id: read_extras.m 443 2006-02-17 16:44:54Z john $


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

if exist(mname,'file'),
    try,
        extras = load(mname);
    catch,
        warning('Can not load "%s" as a binary MAT file.\n', mname);
    end;
end;

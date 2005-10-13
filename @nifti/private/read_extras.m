function extras = read_extras(fname)
% Read extra bits of information
%_______________________________________________________________________
% Copyright (C) 2005 Wellcome Department of Imaging Neuroscience

%
% $Id: read_extras.m 253 2005-10-13 15:31:34Z guillaume $


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
    extras = load(mname);
end;

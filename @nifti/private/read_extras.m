function extras = read_extras(fname)
% Read extra bits of information
% _______________________________________________________________________
% @(#)read_extras.m	1.2 John Ashburner 04/11/26

extras = [];
[pth,nam,ext] = fileparts(fname);
switch ext
case {'.hdr','.img','.nii'}
    mname = fullfile(pth,[nam '.mat']);
case {'.HDR','.IMG','.NII'}
    mname = fullfile(pth,[nam '.MAT']);
otherwise
    mname = fullfile(pth,[nam '.mat']);
end

if exist(mname) == 2,
    extras = load(mname);
end;

function extras = write_extras(fname,extras)
% Write extra bits of information
%__________________________________________________________________________
% Copyright (C) 2005-2011 Wellcome Trust Centre for Neuroimaging

%
% $Id: write_extras.m 4492 2011-09-16 12:11:09Z guillaume $


[pth,nam,ext] = fileparts(fname);
switch ext
    case {'.hdr','.img','.nii'}
        mname = fullfile(pth,[nam '.mat']);
    case {'.HDR','.IMG','.NII'}
        mname = fullfile(pth,[nam '.MAT']);
    otherwise
        mname = fullfile(pth,[nam '.mat']);
end
if isstruct(extras) && ~isempty(fieldnames(extras))
    save(mname,'-struct',extras, spm_get_defaults('mat.format'));
end

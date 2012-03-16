function extras = write_extras(fname,extras)
% Write extra bits of information
%__________________________________________________________________________
% Copyright (C) 2005-2012 Wellcome Trust Centre for Neuroimaging

%
% $Id: write_extras.m 4692 2012-03-16 16:27:51Z guillaume $


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
    try
        opt = spm_get_defaults('mat.format');
    catch
        opt = '-v6';
    end
    save(mname,'-struct','extras', opt);
end

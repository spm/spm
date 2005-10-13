function extras = write_extras(fname,extras)
% Write extra bits of information
%_______________________________________________________________________
% Copyright (C) 2005 Wellcome Department of Imaging Neuroscience

%
% $Id: write_extras.m 253 2005-10-13 15:31:34Z guillaume $


[pth,nam,ext] = fileparts(fname);
switch ext
case {'.hdr','.img','.nii'}
    mname = fullfile(pth,[nam '.mat']);
case {'.HDR','.IMG','.NII'}
    mname = fullfile(pth,[nam '.MAT']);
otherwise
    mname = fullfile(pth,[nam '.mat']);
end
if isstruct(extras) && ~isempty(fieldnames(extras)),
    savefields(mname,extras);
end;

function savefields(fnam,p)
if length(p)>1, error('Can''t save fields.'); end;
fn = fieldnames(p);
for i_=1:length(fn),
    eval([fn{i_} '= p.' fn{i_} ';']);
end;
if numel(fn)>0,
    save(fnam,fn{:});
end;
return;


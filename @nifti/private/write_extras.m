function extras = write_extras(fname,extras)
% Write extra bits of information
% _______________________________________________________________________
% @(#)write_extras.m	1.2 John Ashburner 04/11/26

[pth,nam,ext] = fileparts(fname);
switch ext
case {'.hdr','.img','.nii'}
    mname = fullfile(pth,[nam '.mat']);
case {'.HDR','.IMG','.NII'}
    mname = fullfile(pth,[nam '.MAT']);
otherwise
    mname = fullfile(pth,[nam '.mat']);
end

if isstruct(extras),
    savefields(mname,extras);
end;

function savefields(fnam,p)
if length(p)>1, error('Can''t save fields.'); end;
fn = fieldnames(p);
for i_=1:length(fn),
    eval([fn{i} '= p.' fn{i} ';']);
end;
save(fnam,fn{:});
return;


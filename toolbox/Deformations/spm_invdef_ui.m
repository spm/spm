function spm_invdef_ui
% Writes the inverse of a deformation field.
% It requires a deformation field in the form of three images, and
% a further image from which to derive various dimensions etc.
% The inverse deformation field is written to the same directory
% as the original deformation field, but with "i" prefixed to the
% filenames.
%_______________________________________________________________________
% %W% John Ashburner %E%

n       = spm_input('Number of subjects','+0', 'n', '1', 1)';
for i=1:n,
	P{i}    = spm_get(3,'*y?_*.img',['Select deformation field ' num2str(i)]);
	PT{i}   = spm_get(1,'*.img',['Image to base inverse (' num2str(i) ') on']);
end;

spm_progress_bar('Init',n,'Inverting deformations','volumes completed');
for i=1:n,
	doit(P{i},PT{i});
	spm_progress_bar('Set',i);
end;
spm_progress_bar('Clear')
return;
%_______________________________________________________________________

%_______________________________________________________________________
function doit(P,PT)
V  = spm_vol(P);
y1 = spm_load_float(V(1));
y2 = spm_load_float(V(2));
y3 = spm_load_float(V(3));

VT            = spm_vol(PT);
[iy1,iy2,iy3] = spm_invdef(y1,y2,y3,VT.dim(1:3),inv(VT.mat),V(1).mat);

VO         = VT;
VO.fname   = prepend(V(1).fname, 'i');
VO.dim(4)  = spm_type('float');
VO.pinfo   = [1 0 0]';
VO.descrip = 'Inverse deformation field - X';
spm_write_vol(VO,iy1);

VO         = VT;
VO.fname   = prepend(V(2).fname, 'i');
VO.dim(4)  = spm_type('float');
VO.pinfo   = [1 0 0]';
VO.descrip = 'Inverse deformation field - Y';
spm_write_vol(VO,iy2);

VO         = VT;
VO.fname   = prepend(V(3).fname, 'i');
VO.dim(4)  = spm_type('float');
VO.pinfo   = [1 0 0]';
VO.descrip = 'Inverse deformation field - Y';
spm_write_vol(VO,iy3);
return;
%_______________________________________________________________________

%_______________________________________________________________________
function out = prepend(in, pre)
[pth,nme,ext,ver] = fileparts(in);
out = fullfile(pth,[pre nme ext ver]);
return;
%_______________________________________________________________________

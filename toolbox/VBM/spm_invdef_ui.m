function spm_invdef_ui
% Writes the inverse of a deformation field.
% It requires a deformation field in the form of three images, and
% a further image from which to derive various dimensions etc.
% The inverse deformation field is written to the same directory
% as the original deformation field, but with "i" prefixed to the
% filenames.
%_______________________________________________________________________
% @(#)spm_invdef_ui.m	1.2 John Ashburner 02/08/16

P    = spm_get(Inf,{'*y_*.img','noexpand'},'Select deformation fields');
n    = size(P,1);
PT   = spm_get(n,'*.img',['Images to base inverses on']);

spm_progress_bar('Init',n,'Inverting deformations','volumes completed');
for i=1:n,
	Pi = [repmat([deblank(P(i,:)) ','],3,1) num2str([1 2 3]')];
	doit(Pi,PT(i,:));
	spm_progress_bar('Set',i);
end;
spm_progress_bar('Clear')
return;
%_______________________________________________________________________

%_______________________________________________________________________
function doit(V,VT)
if ischar(V),  V  = spm_vol(V);  end;
if ischar(VT), VT = spm_vol(VT); end;

y1 = spm_load_float(V(1));
y2 = spm_load_float(V(2));
y3 = spm_load_float(V(3));

[iy1,iy2,iy3] = spm_invdef(y1,y2,y3,VT.dim(1:3),inv(VT.mat),V(1).mat);

VO              = V;
VO(1).fname     = prepend(VO(1).fname, 'i');
VO(1).dim(1:3)  = VT.dim(1:3);
VO(1).pinfo     = [1 0 0]';
VO(1).mat       = VT.mat;
VO(1).descrip   = 'Inverse deformation field';
spm_write_vol(VO(1),iy1);

VO(2).fname     = prepend(VO(2).fname, 'i');
VO(2).dim(1:3)  = VT.dim(1:3);
VO(2).pinfo     = [1 0 0]';
VO(2).mat       = VT.mat;
VO(2).descrip   = 'Inverse deformation field';
spm_write_vol(VO(2),iy2);                                                                                   

VO(3).fname     = prepend(VO(3).fname, 'i');
VO(3).dim(1:3)  = VT.dim(1:3);
VO(3).pinfo     = [1 0 0]';
VO(3).mat       = VT.mat;
VO(3).descrip   = 'Inverse deformation field';
spm_write_vol(VO(3),iy3);
return;
%_______________________________________________________________________

%_______________________________________________________________________
function out = prepend(in, pre)
[pth,nme,ext,ver] = fileparts(in);
out = fullfile(pth,[pre nme ext ver]);
return;
%_______________________________________________________________________

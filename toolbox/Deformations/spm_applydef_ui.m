function spm_applydef_ui
% Applies a deformation field to an image
%_______________________________________________________________________
% %W% John Ashburner %E%

n       = spm_input('Number of subjects','+0', 'n', '1', 1)';
for i=1:n,
        P{i}    = spm_get(3,'*y?_*.img',['Select deformation field ' num2str(i)]);
;
        PT{i}   = spm_get(Inf,'*.img',['Image(s) to warp (' num2str(i) ')'
]);
end;

spm_progress_bar('Init',n,'Applying deformations','subjects completed');
for i=1:n,
        doit(P{i},PT{i});
        spm_progress_bar('Set',i);
end;
spm_progress_bar('Clear')
return;
%_______________________________________________________________________

%_______________________________________________________________________
function doit(PD,PI)
VD = spm_vol(PD);
VI = spm_vol(PI);
VO = VI;
for i=1:length(VO),
	VO(i).fname    = prepend(VO(i).fname,'n');
	VO(i).dim(1:3) = VD(1).dim(1:3);
	VO(i).mat      = VD(1).mat;
	if ~isfield(VO,'descrip'), VO(i).descrip = ''; end;
	VO(i).descrip  = ['warped' VO(i).descrip];
end;
VO = spm_create_vol(VO(i));

for p=1:VD(1).dim(3),
	M  = spm_matrix([0 0 p]);
	x1 = spm_slice_vol(VD(1), M, VD(1).dim(1:2),1);
	x2 = spm_slice_vol(VD(2), M, VD(1).dim(1:2),1);
	x3 = spm_slice_vol(VD(3), M, VD(1).dim(1:2),1);
	for i=1:length(VI),
		M     = inv(VI(i).mat);
		y1    = M(1,1)*x1+M(1,2)*x2+M(1,3)*x3+M(1,4);
		y2    = M(2,1)*x1+M(2,2)*x2+M(2,3)*x3+M(2,4);
		y3    = M(3,1)*x1+M(3,2)*x2+M(3,3)*x3+M(3,4);
		img   = spm_sample_vol(VI(i),y1,y2,y3,1);
		VO(i) = spm_write_plane(VO(i),img,p);
	end;
end;
VO = spm_close_vol(VO);
return;
%_______________________________________________________________________

%_______________________________________________________________________
function out = prepend(in, pre)
[pth,nme,ext,ver] = fileparts(in);
out = fullfile(pth,[pre nme ext ver]);
return;
%_______________________________________________________________________

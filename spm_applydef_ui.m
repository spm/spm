function spm_applydef_ui(P,PT)
% Applies a deformation field to an image
%_______________________________________________________________________
% Copyright (C) 2005 Wellcome Department of Imaging Neuroscience

% John Ashburner
% $Id: spm_applydef_ui.m 112 2005-05-04 18:20:52Z john $


if nargin<2
	n = spm_input('Number of subjects','+0', 'n', '1', 1)';
	for i=1:n
		P{i}    = spm_select(1,'.*y_.*\.img',['Select deformation field ' num2str(i)]);
		PT{i}   = spm_select(Inf,'image',['Image(s) to warp (' num2str(i) ')']);
	end;
else
	n = length(P);
	if n ~=length(PT)
		error('Must have matching deformation / files list cell arrays as inputs to spm_applydef_ui')
	end
end

spm_progress_bar('Init',n,'Applying deformations','subjects completed');
for i=1:length(P),
	Pi = [repmat([P{i} ','],3,1) num2str([1 2 3]')];
        spm_applydef(Pi,PT{i});
        spm_progress_bar('Set',i);
end;
spm_progress_bar('Clear')
return;
%_______________________________________________________________________

%_______________________________________________________________________
function spm_applydef(VD,VI)
if ischar(VD), VD = spm_vol(VD); end;
if ischar(VI), VI = spm_vol(VI); end;
VO = VI;
for i=1:length(VO),
	VO(i).fname    = prepend(VO(i).fname,'w');
	VO(i).dim(1:3) = VD(1).dim(1:3);
	VO(i).mat      = VD(1).mat;
	if ~isfield(VO,'descrip'), VO(i).descrip = ''; end;
	VO(i).descrip  = ['warped ' VO(i).descrip];
end;
VO = spm_create_vol(VO);

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
return;
%_______________________________________________________________________

%_______________________________________________________________________
function out = prepend(in, pre)
[pth,nme,ext,ver] = fileparts(in);
out = fullfile(pth,[pre nme ext ver]);
return;
%_______________________________________________________________________

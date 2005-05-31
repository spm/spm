function spm_combdef_ui
% Combines deformation fields
%_______________________________________________________________________
% Copyright (C) 2005 Wellcome Department of Imaging Neuroscience

% John Ashburner
% $Id: spm_combdef_ui.m 184 2005-05-31 13:23:32Z john $


P = spm_select(Inf,'.*y_.*\.img$','Select deformation fields');
V = cell(size(P,1),1);
for i=1:size(P,1),
	V{i}   = spm_vol([repmat([deblank(P(i,:)) ','],3,1) num2str([1 2 3]')]);
end;
hld = 1;

VO = V{1};
for i=1:length(VO),
	VO(i).fname = prepend(VO(i).fname,'c');
	VO(i).desc  = 'Combined deformation field';
end;

VO = spm_create_vol(VO);

spm_progress_bar('Init',VO(1).dim(3),'Combining deformations','planes completed');
for p=1:VO(1).dim(3),
	M  = spm_matrix([0 0 p]);
	y1 = spm_slice_vol(V{1}(1), M, V{1}(1).dim(1:2),1);
	y2 = spm_slice_vol(V{1}(2), M, V{1}(1).dim(1:2),1);
	y3 = spm_slice_vol(V{1}(3), M, V{1}(1).dim(1:2),1);
	for i=2:length(V)
		M   = inv(V{i}(1).mat);
		ty1 = M(1,1)*y1+M(1,2)*y2+M(1,3)*y3+M(1,4);
		ty2 = M(2,1)*y1+M(2,2)*y2+M(2,3)*y3+M(2,4);
		ty3 = M(3,1)*y1+M(3,2)*y2+M(3,3)*y3+M(3,4);
		y1 = spm_sample_vol(V{i}(1),ty1,ty2,ty3,[hld NaN]);
		y2 = spm_sample_vol(V{i}(2),ty1,ty2,ty3,[hld NaN]);
		y3 = spm_sample_vol(V{i}(3),ty1,ty2,ty3,[hld NaN]);
	end;
	VO(1) = spm_write_plane(VO(1),y1,p);
	VO(2) = spm_write_plane(VO(2),y2,p);
	VO(3) = spm_write_plane(VO(3),y3,p);
        spm_progress_bar('Set',p);
end;
spm_progress_bar('Clear')
return;
%_______________________________________________________________________

%_______________________________________________________________________
function out = prepend(in, pre)
[pth,nme,ext,ver] = fileparts(in);
out = fullfile(pth,[pre nme ext ver]);
return;
%_______________________________________________________________________

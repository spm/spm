function spm_procrustes_ui
% Removes pose (and global size effects) from deformation fields.
% Corrected fields are written prefixed by 'p'.
%_______________________________________________________________________
% Example Reference:
% F. L. Bookstein (1997).  "Landmark Methods for Forms Without
% Landmarks: Morphometrics of Group Differences in Outline Shape"
% Medical Image Analysis 1(3):225-243
%_______________________________________________________________________
% %W% John Ashburner %E%

P   = spm_get(Inf,{'*y_*.img','noexpand'},'Select deformation fields');
PW  = spm_get(1,'*.img','Weighting image');
flg = spm_input('Remove size?','+1','y/n',[1,0],1);
n   = size(P,1);
spm_progress_bar('Init',n,'Tweaking deformations','volumes completed');
for i=1:n,
	Pi = spm_vol([repmat([deblank(P(i,:)) ','],3,1) num2str([1 2 3]')]);
	doit(Pi,PW,flg);
	spm_progress_bar('Set',i);
end;
spm_progress_bar('Clear')
return;
%_______________________________________________________________________

%_______________________________________________________________________
function doit(V,VW,flg)
if ischar(V),  V  = spm_vol(V);  end;
if ischar(VW), VW = spm_vol(VW); end;
y1 = spm_load_float(V(1));
y2 = spm_load_float(V(2));
y3 = spm_load_float(V(3));

T = procrustes(y1,y2,y3,V(1).mat,VW,flg);
spm_affdef(y1,y2,y3,T);
disp(T)

VO         = V(1);
VO.fname   = prepend(V(1).fname, 'p');
VO.dim(4)  = spm_type('float');
VO.pinfo   = [1 0 0]';
VO.descrip = 'Shape field';
spm_write_vol(VO,y1);

VO         = V(2);
VO.fname   = prepend(V(2).fname, 'p');
VO.dim(4)  = spm_type('float');
VO.pinfo   = [1 0 0]';
VO.descrip = 'Shape field';
spm_write_vol(VO,y2);

VO         = V(3);
VO.fname   = prepend(V(3).fname, 'p');
VO.dim(4)  = spm_type('float');
VO.pinfo   = [1 0 0]';
VO.descrip = 'Shape field';
spm_write_vol(VO,y3);
return;
%_______________________________________________________________________

%_______________________________________________________________________
function out = prepend(in, pre)
[pth,nme,ext,ver] = fileparts(in);
out = fullfile(pth,[pre nme ext ver]);
return;
%_______________________________________________________________________

%_______________________________________________________________________
function T = procrustes(Y1,Y2,Y3,M1,V,flg)
% Compute the affine transformation that minimises the
% Procrustes distance between transformed and un-transformed
% points.
%
% FORMAT T = procrustes(Y1,Y2,Y3,M1,V,flg)
% Y1, Y2 and Y3 - spatial mapping
% M1            - voxel to mm mapping of Y1, Y2 and Y3.
% V             - a handle for a weighting image.
% flg           - a flag to say if zooming should be included.
%                 1 means no
% T             - the determined affine transformation matrix.
%_______________________________________________________________________
% Example Reference:
% F. L. Bookstein (1997).  "Landmark Methods for Forms Without
% Landmarks: Morphometrics of Group Differences in Outline Shape"
% Medical Image Analysis 1(3):225-243
%_______________________________________________________________________
% %W% John Ashburner %E%

[z1,z2] = ndgrid(1:size(Y1,1),1:size(Y1,2));
z1      = z1(:);
z2      = z2(:);

%_______________________________________________________________________

% Compute centres of mass.
sm = 0;
c1 = [0 0 0];
c2 = [0 0 0];
d  = size(Y1);
for j=1:size(Y1,3),
	M  = V.mat\M1*spm_matrix([0 0 j]);
	wt = spm_slice_vol(V,M,d(1:2),1);
	wt = wt(:);

	y1 = double(Y1(:,:,j)); y1 = y1(:);
	y2 = double(Y2(:,:,j)); y2 = y2(:);
	y3 = double(Y3(:,:,j)); y3 = y3(:);

	x1 = M1(1,1)*z1 + M1(1,2)*z2 + (M1(1,3)*j + M1(1,4)); 
	x2 = M1(2,1)*z1 + M1(2,2)*z2 + (M1(2,3)*j + M1(2,4));
	x3 = M1(3,1)*z1 + M1(3,2)*z2 + (M1(3,3)*j + M1(3,4));

	sm = sm + sum(wt);
	c1 = c1 + sum([x1.*wt x2.*wt x3.*wt]);
	c2 = c2 + sum([y1.*wt y2.*wt y3.*wt]);
end;
c1 = c1/sm;
c2 = c2/sm;

%_______________________________________________________________________

if nargin<6 | flg,
	% Compute centroid size.
	sm   = 0;
	mom1 = 0;
	mom2 = 0;
	for j=1:size(Y1,3),
		M  = V.mat\M1*spm_matrix([0 0 j]);
		wt = spm_slice_vol(V,M,d(1:2),1);
		wt = wt(:);

		y1 = double(Y1(:,:,j)); y1 = y1(:);
		y2 = double(Y2(:,:,j)); y2 = y2(:);
		y3 = double(Y3(:,:,j)); y3 = y3(:);

		x1 = M1(1,1)*z1 + M1(1,2)*z2 + (M1(1,3)*j + M1(1,4)); 
		x2 = M1(2,1)*z1 + M1(2,2)*z2 + (M1(2,3)*j + M1(2,4));
		x3 = M1(3,1)*z1 + M1(3,2)*z2 + (M1(3,3)*j + M1(3,4));

		sm   = sm + sum(wt);
		mom1 = mom1 + sum([(x1-c1(1)).^2.*wt  (x2-c1(2)).^2.*wt (x3-c1(3)).^2.*wt]);
		mom2 = mom2 + sum([(y1-c2(1)).^2.*wt  (y2-c2(2)).^2.*wt (y3-c2(3)).^2.*wt]);
	end;
	zm = sqrt(mom1/mom2);
else,
	zm = 1;
end;
ZM = spm_matrix([0 0 0  0 0 0  zm zm zm  0 0 0]);

%_______________________________________________________________________

% Compute rotations
C = zeros(3);
for j=1:size(Y1,3),
	M  = V.mat\M1*spm_matrix([0 0 j]);
	wt = spm_slice_vol(V,M,d(1:2),1);
	wt = wt(:);

	y1 = double(Y1(:,:,j)); y1 = y1(:);
	y2 = double(Y2(:,:,j)); y2 = y2(:);
	y3 = double(Y3(:,:,j)); y3 = y3(:);

	x1 = M1(1,1)*z1 + M1(1,2)*z2 + (M1(1,3)*j + M1(1,4)); 
	x2 = M1(2,1)*z1 + M1(2,2)*z2 + (M1(2,3)*j + M1(2,4));
	x3 = M1(3,1)*z1 + M1(3,2)*z2 + (M1(3,3)*j + M1(3,4));

	wt1 = wt*zm;
	C = C + [(x1-c1(1)).*wt  (x2-c1(2)).*wt  (x3-c1(3)).*wt ]' * ...
	        [(y1-c2(1)).*wt1 (y2-c2(2)).*wt1 (y3-c2(3)).*wt1];
end;
[u,s,v]    = svd(C);
R          = eye(4);
R(1:3,1:3) = u*v';

%_______________________________________________________________________

% Combine rotations, translations (and zooms) into a single
% affine transformation matrix.
T = spm_matrix(c1)*R*ZM*spm_matrix(-c2);
return;
%_______________________________________________________________________

function spm_sn2def(varargin)
% Deformation field utilities (works on sn.mat files).
%
% This toolbox includes a number of utilities for extracting information
% from _sn.mat files created during spatial normalisation.  Much of
% the information can be used for either deformation-based morphometry
% (DBM), or tensor-based morphometry (TBM).
%
% Deformations:
% Writes out the deformation fields as y_*.img.
% These fields are the voxel to voxel mappings between an image
% normalised with the specified bounding-box and voxel sizes, and the
% original images.  The deformations can be used for multivariate
% morphometric methods (DBM), but they will need some initial
% corrections for voxel-sizes and for some measure of global position
% and possibly size (Procrustes shape).
%
% Jacobian Matrices:
% Writes out the Jacobian matrix field as a single file.  The order
% of the elements are 11, 12, 13, 21, 22, 23, 31, 32, 33.
%
% Jacobian Determinants:
% Writes out the volume change at each location in j_*.img.  Fields
% involving a flip should have begative Jacobian determinants, but
% these are implicitly made positive to make things easier to
% understand.  The determinants can be used for morphometric methods
% (TBM) that characterise volumetric differences.
%
% Strain Tensors:
% These are intended for voxelwise multivariate morphometry (TBM).
% Much of the notation and ideas are from:
%    "Non-linear Elastic Deformations", by R. W. Ogden.  Dover
%    Publications, 1984.
%
% The transformation maps from elements in the template (x1,x2,x3)
% to elements in the original images (y1,y2,y3).
% An affine mapping from x to y is given by:
% y = A*x + c, where A is the deformation gradient (second order tensor
% which is the Jacobian matrix) and c is a constant representing
% translations.
% We wish to represent shape changes within the co-ordinate framework of
% the template images (Lagrangean framework, as opposed to the Eulerian
% framework where deformations are relative to the individual images).
% To do this, matrix A is decomposed (using polar decomposition), such
% that A = R*U.  Matrix R is a rigid body transformation matrix, and U
% is a matrix that is purely shears and zooms (no rotations).
% This gives a mapping from template voxels to image voxels of the form
% y = R*U*x + c,  which means zoom and shear the positions in the template,
% and then do a rigid body rotation (and also add the translation).
% We are not interested in the translations and rotations - only in the
% local zooms and shears.  Therefore all the information we need is in
% the matrix U, which is simply obtained by U=(A'*A)^(1/2).
% Note that the zooms and shears are done while in the orientation of
% the template, before rotations to the orientation of the images are
% introduced.
%
% Strain tensors are defined that model the amount of distortion.  If
% there is no strain, then the tensors are all zero. Generically,
% the family of Lagrangean strain tensors are given by:
% (U^m-eye(3))/m when m~=0, and logm(U) when m==0.
%
%
%_______________________________________________________________________
% @(#)spm_sn2def.m	1.4 John Ashburner 03/11/26

global defaults
defs = defaults.normalise.write;

if nargin==0,
	allowable = {'def','jacmat','jacdet','tensor'};
	a1 = spm_input('Write what?',1,'m',...
		['Deformations|Jacobian Matrices|Jacobian Determinants|' ...
		 'Strain Tensors'],...
		[1 2 3 4],1);
	arg = allowable{a1};

	if strcmp(arg,'tensor')
		tensorder = spm_input('Write what?',2,'m',...
			['Strain Tensors, m=-2 (Almansi)|'...
			 'Strain Tensors, m=-1          |'...
			 'Strain Tensors, m= 0 (Hencky) |'...
			 'Strain Tensors, m= 1 (Biot)   |'...
			 'Strain Tensors, m= 2 (Green)   '],...
			[-2 -1 0 1 2],3);
	end;
else,
	arg = lower(varargin{1});
	switch arg,
		case 'tensor',
			if nargin >1,
				tensorder = varargin{2};
			else,
				tensorder = 0;
			end;
		case {'def','jacmat','jacdet'},
		otherwise,
			error('Unknown argument');
	end;
end;

files = spm_get(Inf,'*_sn.mat','Select *_sn.mat files');

for i=1:size(files,1),
	sn(i) = load(deblank(files(i,:)));
end;
for i=1:size(files,1),
	sn(i).fname = deblank(files(i,:));
end;

switch arg,
	case 'tensor'
		spm_progress_bar('Init',length(sn),...
			['Computing Tensors (m= ' num2str(tensorder) ')'],...
			'volumes completed');
	case 'def'
		spm_progress_bar('Init',length(sn),...
			['Computing Deformations'],'volumes completed');
	case 'jacmat'
		spm_progress_bar('Init',length(sn),...
			['Computing Jacobian Matrices'],...
			'volumes completed');
	case 'jacdet'
		spm_progress_bar('Init',length(sn),...
			['Computing Jacobian Determinants'],...
			'volumes completed');
end;
spm('Pointer','Watch')
for i=1:length(sn),
	%spm('FigName',['Deformations: working on subj ' num2str(i)],Finter,CmdLine);drawnow;
	switch arg,
		case 'tensor',
			write_tensor(sn(i),defs.vox,defs.bb,tensorder);
		case 'def',
			spm_write_defs(sn(i),defs.vox,defs.bb);
		case 'jacmat',
			write_jacobianm(sn(i),defs.vox,defs.bb);
		case 'jacdet',
			write_det(sn(i),defs.vox,defs.bb);
	end;
	spm_progress_bar('Set',i);
end;
spm_progress_bar('Clear');
%spm('FigName','Deformations: done',Finter,CmdLine);
spm('Pointer');

return;
%_______________________________________________________________________

%_______________________________________________________________________
function write_jacobianm(sn,vox,bb)
vo = init_vo(sn,vox,bb);
[pth,nm,xt,vr]  = fileparts(deblank(sn.VF.fname));
vo.fname   = fullfile(pth,['a_' nm '.img']);
vo.dim(4)  = spm_type('float');
vo.pinfo   = [1 0 0]';
vo.descrip = ['Jacobian_Matrix'];
vo.n       = 1;
for i=1:3,
	for j=1:3,
		VO(i,j)          = vo;
		VO(i,j).n        = (i-1)*3+j;
		%VO(i,j).fname   = fullfile(pth,['a' num2str(i) num2str(j) '_' nm '.img']);
		%VO(i,j).descrip = ['Jacobian_Matrix(' num2str(i) ',' num2str(j) ')'];
		%VO(i,j).n       = 1;
	end;
end;
VO = spm_create_vol(VO);
for p=1:vo.dim(3),
	A = get_jacobian(sn, vo, p);
	for i=1:3,
		for j=1:3,
			VO(i,j) = spm_write_plane(VO(i,j),A(:,:,i,j),p);
		end;
	end;
end;
VO = spm_close_vol(VO);
return;
%_______________________________________________________________________

%_______________________________________________________________________
function write_det(sn,vox,bb)
VO = init_vo(sn,vox,bb);
[pth,nm,xt,vr]  = fileparts(deblank(sn.VF.fname));
VO.fname   = fullfile(pth,['j_' nm '.img']);
VO.dim(4)  = spm_type('float');
VO.pinfo   = [1 0 0]';
VO.descrip = ['Jacobian_Determinant'];
VO         = spm_create_vol(VO);

for p=1:VO.dim(3),
	A = get_jacobian(sn, VO(1,1), p);
	dtA = A(:,:,1,1).*(A(:,:,2,2).*A(:,:,3,3) - A(:,:,2,3).*A(:,:,3,2)) ...
	    - A(:,:,2,1).*(A(:,:,1,2).*A(:,:,3,3) - A(:,:,1,3).*A(:,:,3,2)) ...
	    + A(:,:,3,1).*(A(:,:,1,2).*A(:,:,2,3) - A(:,:,1,3).*A(:,:,2,2));
	VO = spm_write_plane(VO,dtA,p);
end;
VO = spm_close_vol(VO);
return;
%_______________________________________________________________________

%_______________________________________________________________________
function write_tensor(sn,vox,bb,m)
vo = init_vo(sn,vox,bb);
%ij = [1 1; 2 1; 3 1; 2 2; 3 2; 3 3];
[pth,nm,xt,vr] = fileparts(deblank(sn.VF.fname));
vo.fname    = fullfile(pth,['e' num2str(m) '_' nm '.img']);
vo.descrip  =  ['Strain_Tensor - m=' num2str(m)];
vo.dim(4)   = spm_type('float');
vo.pinfo    = [1 0 0]';
vo.n        = 1;
for i=1:6,
	VO(i)          = vo;
	VO(i).n        = i;
	%VO(i).fname   = fullfile(pth,['e' num2str(ij(i,1)) num2str(ij(i,2)) 'm' num2str(m) '_' nm '.img']);
	%VO(i).descrip = ['Strain_Tensor(' num2str(ij(i,1)) ',' num2str(ij(i,2)) ') - m=' num2str(m)];
	%VO(i).n       = 1;
end;
VO = spm_create_vol(VO);

for p=1:vo.dim(3),
	A = get_jacobian(sn, vo, p);
	E = tensors_from_jacobian(A,m);
	for i=1:6,
		VO(i) = spm_write_plane(VO(i),E(:,:,i),p);
	end;
end;
VO = spm_close_vol(VO);
return;
%_______________________________________________________________________

%_______________________________________________________________________
function VO = init_vo(sn,vox,bb)
if nargin>=2,
	[bb0,vox0] = bbvox_from_V(sn.VG);

	if any(~finite(bb)),  bb  = bb0;  end;
	if any(~finite(vox)), vox = vox0; end;

	bb      = sort(bb);
	vox     = abs(vox);

	M       = sn.VG(1).mat;
	ogn     = M\[0 0 0 1]';
	ogn     = ogn(1:3)';
	vxg     = sqrt(sum(M(1:3,1:3).^2));
	bb(:,1) = round(bb(:,1)/vox(1))*vox(1);
	bb(:,2) = round(bb(:,2)/vox(2))*vox(2);
	bb(:,3) = round(bb(:,3)/vox(3))*vox(3);
        x       = (bb(1,1):vox(1):bb(2,1))/vxg(1) + ogn(1);
        y       = (bb(1,2):vox(2):bb(2,2))/vxg(2) + ogn(2);
        z       = (bb(1,3):vox(3):bb(2,3))/vxg(3) + ogn(3);
	og      = -vxg.*ogn;
	of      = -vox.*(round(-bb(1,:)./vox)+1);
	M1      = [vxg(1) 0 0 og(1) ; 0 vxg(2) 0 og(2) ; 0 0 vxg(3) og(3) ; 0 0 0 1];
	M2      = [vox(1) 0 0 of(1) ; 0 vox(2) 0 of(2) ; 0 0 vox(3) of(3) ; 0 0 0 1];
        dim     = [length(x) length(y) length(z)];
	mat     = sn.VG.mat*inv(M1)*M2;
else,
        dim     = sn.VG.dim(1:3);
        x       = 1:dim(1);
        y       = 1:dim(2);
        z       = 1:dim(3);

        mat     = sn.VG.mat;
end;
VO     = struct('fname','',...
        'dim',[dim NaN],  'mat',mat,...
        'pinfo',[NaN NaN NaN]',...
        'descrip','');

return;
%_______________________________________________________________________

%_______________________________________________________________________
function A = get_jacobian(sn, VO, j)
% Each element of the Jacobian matrix field consists of something analagous to the following:
% maple diff( m11*(x1+p11*b1(x1,x2)+p12*b2(x1,x2)) + m12*(x2+p21*b1(x1,x2)+p22*b2(x1,x2)) + m13, x1)
% maple diff( m11*(x1+p11*b1(x1,x2)+p12*b2(x1,x2)) + m12*(x2+p21*b1(x1,x2)+p22*b2(x1,x2)) + m13, x2)
% maple diff( m21*(x1+p11*b1(x1,x2)+p12*b2(x1,x2)) + m22*(x2+p21*b1(x1,x2)+p22*b2(x1,x2)) + m23, x1)
% maple diff( m21*(x1+p11*b1(x1,x2)+p12*b2(x1,x2)) + m22*(x2+p21*b1(x1,x2)+p22*b2(x1,x2)) + m23, x2)


dim       = VO.dim(1:3);
mat       = VO.mat;
st        = size(sn.Tr);

% Assume transverse images, and obtain position of pixels in millimeters,
% and convert to voxel space of template.
%----------------------------------------------------------------------------
M   = sn.VG.mat\VO.mat;
x   = M(1,1)*(1:dim(1))+M(1,4);
y   = M(2,2)*(1:dim(2))+M(2,4);
z   = M(3,3)*j         +M(3,4);
X   = x'*ones(1,dim(2));
Y   = ones(dim(1),1)*y;

bbX = spm_dctmtx(sn.VG(1).dim(1),st(1),x-1);
bbY = spm_dctmtx(sn.VG(1).dim(2),st(2),y-1);
bbZ = spm_dctmtx(sn.VG(1).dim(3),st(3),z-1);

dbX = spm_dctmtx(sn.VG(1).dim(1),st(1),x-1,'diff');
dbY = spm_dctmtx(sn.VG(1).dim(2),st(2),y-1,'diff');
dbZ = spm_dctmtx(sn.VG(1).dim(3),st(3),z-1,'diff');

M   = sn.VF.mat*sn.Affine*inv(sn.VG.mat);

if prod(st)>0,
	% Nonlinear deformations
	%----------------------------------------------------------------------------
	% 2D transforms for each plane
	tbx = reshape( reshape(sn.Tr(:,:,:,1),st(1)*st(2),st(3)) *bbZ', st(1), st(2) );
	tby = reshape( reshape(sn.Tr(:,:,:,2),st(1)*st(2),st(3)) *bbZ', st(1), st(2) );
	tbz = reshape( reshape(sn.Tr(:,:,:,3),st(1)*st(2),st(3)) *bbZ', st(1), st(2) );

	tdx = reshape( reshape(sn.Tr(:,:,:,1),st(1)*st(2),st(3)) *dbZ', st(1), st(2) );
	tdy = reshape( reshape(sn.Tr(:,:,:,2),st(1)*st(2),st(3)) *dbZ', st(1), st(2) );
	tdz = reshape( reshape(sn.Tr(:,:,:,3),st(1)*st(2),st(3)) *dbZ', st(1), st(2) );

	% Jacobian of transformation from template
	% to affine registered image.
	%---------------------------------------------
	j11 = dbX*tbx*bbY' + 1;
	j12 = bbX*tbx*dbY';
	j13 = bbX*tdx*bbY';

	j21 = dbX*tby*bbY';
	j22 = bbX*tby*dbY' + 1;
	j23 = bbX*tdy*bbY';

	j31 = dbX*tbz*bbY';
	j32 = bbX*tbz*dbY';
	j33 = bbX*tdz*bbY' + 1;
else,
	j11 = ones(size(bbX,1),size(bbY,1));
	j12 = zeros(size(bbX,1),size(bbY,1));
	j13 = zeros(size(bbX,1),size(bbY,1));

	j21 = zeros(size(bbX,1),size(bbY,1));
	j22 = ones(size(bbX,1),size(bbY,1));
	j23 = zeros(size(bbX,1),size(bbY,1));

	j31 = zeros(size(bbX,1),size(bbY,1));
	j32 = zeros(size(bbX,1),size(bbY,1));
	j33 = ones(size(bbX,1),size(bbY,1));
end;

% Combine Jacobian of transformation from
% template to affine registered image, with
% Jacobian of transformation from affine
% registered image to original image.
%---------------------------------------------
A = zeros([size(j11) 3 3]);

A(:,:,1,1) = M(1,1)*j11 + M(1,2)*j21 + M(1,3)*j31;
A(:,:,1,2) = M(1,1)*j12 + M(1,2)*j22 + M(1,3)*j32;
A(:,:,1,3) = M(1,1)*j13 + M(1,2)*j23 + M(1,3)*j33;

A(:,:,2,1) = M(2,1)*j11 + M(2,2)*j21 + M(2,3)*j31;
A(:,:,2,2) = M(2,1)*j12 + M(2,2)*j22 + M(2,3)*j32;
A(:,:,2,3) = M(2,1)*j13 + M(2,2)*j23 + M(2,3)*j33;

A(:,:,3,1) = M(3,1)*j11 + M(3,2)*j21 + M(3,3)*j31;
A(:,:,3,2) = M(3,1)*j12 + M(3,2)*j22 + M(3,3)*j32;
A(:,:,3,3) = M(3,1)*j13 + M(3,2)*j23 + M(3,3)*j33;

return;
%_______________________________________________________________________

%_______________________________________________________________________
function tensor = tensors_from_jacobian(A, m)

d      = size(A);
tensor = zeros([d(1:2) 6]);
I      = eye(3);
if m==0,
	% Hencky
	for j=1:d(2),
		for i=1:d(1),
			J = squeeze(A(i,j,:,:));
			T = logm(J'*J)*0.5;
			tensor(i,j,:) = T([1 2 3 5 6 9]);
		end;
	end;
elseif m==2,
	% 
	for j=1:d(2),
		for i=1:d(1),
			J = squeeze(A(i,j,:,:));
			T = 0.5*(J'*J - I);
			tensor(i,j,:) = T([1 2 3 5 6 9]);
		end;
	end;
else,
	for j=1:d(2),
		for i=1:d(1),
			J = squeeze(A(i,j,:,:));
			T = ((J'*J)^(m/2)-I)/m;
			tensor(i,j,:) = T([1 2 3 5 6 9]);
		end;
	end;
end;
return;
%_______________________________________________________________________

%_______________________________________________________________________
function spm_write_defs(sn, vox,bb)
% Write deformation field.
% FORMAT spm_write_defs(sn, vox,bb)
% sn  - information from the `_sn.mat' file containing the spatial
%         normalization parameters.
% The deformations are stored in y1.img, y2.img and y3.img
%_______________________________________________________________________

[bb0,vox0] = bbvox_from_V(sn.VG);
if any(~finite(vox)), vox = vox0; end;
if any(~finite(bb)),  bb  = bb0;  end;
bb  = sort(bb);
vox = abs(vox);

if nargin>=3,

	if any(~finite(vox)), vox = vox0; end;
	if any(~finite(bb)),  bb  = bb0;  end;
	bb  = sort(bb);
	vox = abs(vox);

	% Adjust bounding box slightly - so it rounds to closest voxel.
	bb(:,1) = round(bb(:,1)/vox(1))*vox(1);
	bb(:,2) = round(bb(:,2)/vox(2))*vox(2);
	bb(:,3) = round(bb(:,3)/vox(3))*vox(3);
 
	M   = sn.VG(1).mat;
	vxg = sqrt(sum(M(1:3,1:3).^2));
	ogn = M\[0 0 0 1]';
	ogn = ogn(1:3)';
 
	% Convert range into range of voxels within template image
	x   = (bb(1,1):vox(1):bb(2,1))/vxg(1) + ogn(1);
	y   = (bb(1,2):vox(2):bb(2,2))/vxg(2) + ogn(2);
	z   = (bb(1,3):vox(3):bb(2,3))/vxg(3) + ogn(3);
 
	og  = -vxg.*ogn;
	of  = -vox.*(round(-bb(1,:)./vox)+1);
	M1  = [vxg(1) 0 0 og(1) ; 0 vxg(2) 0 og(2) ; 0 0 vxg(3) og(3) ; 0 0 0 1];
	M2  = [vox(1) 0 0 of(1) ; 0 vox(2) 0 of(2) ; 0 0 vox(3) of(3) ; 0 0 0 1];
	mat = sn.VG.mat*inv(M1)*M2; 
	dim = [length(x) length(y) length(z)];
else,
	dim    = sn.VG.dim(1:3);
	x      = 1:dim(1);
	y      = 1:dim(2);
	z      = 1:dim(3);
	mat    = sn.VG.mat;
end;

[pth,nm,xt,vr]  = fileparts(deblank(sn.VF.fname));
%VX = struct('fname',fullfile(pth,['y1_' nm '.img']),  'dim',[dim 16], ...
%	'mat',mat,  'pinfo',[1 0 0]',  'descrip','Deformation field - X');
%VY = struct('fname',fullfile(pth,['y2_' nm '.img']),  'dim',[dim 16], ...
%	'mat',mat,  'pinfo',[1 0 0]',  'descrip','Deformation field - Y');
%VZ = struct('fname',fullfile(pth,['y3_' nm '.img']),  'dim',[dim 16], ...
%	'mat',mat,  'pinfo',[1 0 0]',  'descrip','Deformation field - Z');

VX = struct('fname',fullfile(pth,['y_' nm '.img']),  'dim',[dim 16], ...
        'mat',mat,  'pinfo',[1 0 0]',  'descrip','Deformation field', 'n',1);
VY = VX; VY.n = 2;
VZ = VX; VZ.n = 3;

X = x'*ones(1,VX.dim(2));
Y = ones(VX.dim(1),1)*y;

st = size(sn.Tr);

if (prod(st) == 0),
	affine_only = 1;
	basX = 0; tx = 0;
	basY = 0; ty = 0;
	basZ = 0; tz = 0;
else,
	affine_only = 0;
	basX = spm_dctmtx(sn.VG(1).dim(1),st(1),x-1);
	basY = spm_dctmtx(sn.VG(1).dim(2),st(2),y-1);
	basZ = spm_dctmtx(sn.VG(1).dim(3),st(3),z-1); 
end,

if strcmp(computer,'PCWIN'),
	% Windows may not like more than one file handle
	% to the same file.

	VX = spm_create_vol(VX,'noopen');
	VY = spm_create_vol(VY,'noopen');
	VZ = spm_create_vol(VZ,'noopen');
else,
	VX = spm_create_vol(VX);
	VY = spm_create_vol(VY);
	VZ = spm_create_vol(VZ);
end;

% Cycle over planes
%----------------------------------------------------------------------------
for j=1:length(z)

	% Nonlinear deformations
	%----------------------------------------------------------------------------
	if (~affine_only)
		% 2D transforms for each plane
		tx = reshape( reshape(sn.Tr(:,:,:,1),st(1)*st(2),st(3)) *basZ(j,:)', st(1), st(2) );
		ty = reshape( reshape(sn.Tr(:,:,:,2),st(1)*st(2),st(3)) *basZ(j,:)', st(1), st(2) );
		tz = reshape( reshape(sn.Tr(:,:,:,3),st(1)*st(2),st(3)) *basZ(j,:)', st(1), st(2) );

		X1 = X    + basX*tx*basY';
		Y1 = Y    + basX*ty*basY';
		Z1 = z(j) + basX*tz*basY';
	end

	% Sample each volume
	%----------------------------------------------------------------------------
	Mult = sn.VF.mat*sn.Affine;
	if (~affine_only)
		X2= Mult(1,1)*X1 + Mult(1,2)*Y1 + Mult(1,3)*Z1 + Mult(1,4);
		Y2= Mult(2,1)*X1 + Mult(2,2)*Y1 + Mult(2,3)*Z1 + Mult(2,4);
		Z2= Mult(3,1)*X1 + Mult(3,2)*Y1 + Mult(3,3)*Z1 + Mult(3,4);
	else
		X2= Mult(1,1)*X + Mult(1,2)*Y + (Mult(1,3)*z(j) + Mult(1,4));
		Y2= Mult(2,1)*X + Mult(2,2)*Y + (Mult(2,3)*z(j) + Mult(2,4));
		Z2= Mult(3,1)*X + Mult(3,2)*Y + (Mult(3,3)*z(j) + Mult(3,4));
	end

	VX = spm_write_plane(VX,X2,j);
	VY = spm_write_plane(VY,Y2,j);
	VZ = spm_write_plane(VZ,Z2,j);

end;

if ~strcmp(computer,'PCWIN'),
	VX = spm_close_vol(VX);
	VY = spm_close_vol(VY);
	VZ = spm_close_vol(VZ);
end;

return;
%_______________________________________________________________________

%_______________________________________________________________________
function [bb,vx] = bbvox_from_V(V)
vx = sqrt(sum(V.mat(1:3,1:3).^2));
o  = V.mat\[0 0 0 1]';
o  = o(1:3)';
bb = [-vx.*(o-1) ; vx.*(V.dim(1:3)-o)];
return;


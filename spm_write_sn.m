function spm_write_sn(V,prm,flags)
% Write Out Normalized Images.
% FORMAT spm_write_sn(V,prm,Hold,bb,vxf)
% V         - Images to transform (filenames or volume structure).
% matname   - Transformation information (filename or structure).
% flags     - flags structure, with fields...
%           interp   - interpolation method (0-7)
%           wrap     - wrap edges (e.g., [1 1 0] for 2D MRI sequences)
%           vox      - voxel sizes (3 element vector - in mm)
%                      Non-finite values mean use template vox.
%           bb       - bounding box (2x3 matrix - in mm)
%                      Non-finite values mean use template bb.
%           preserve - either 0 or 1.  A value of 1 will "modulate"
%                      the spatially normalised images so that total
%                      units are preserved, rather than just
%                      concentrations.
%
% Non-finite vox or bounding box suggests that values should be derived
% from the template image.
%
% Don't use interpolation methods greater than one for data containing
% NaNs.
%_______________________________________________________________________
% %W% John Ashburner %E%

if ischar(prm), prm = load(prm);  end;
if ischar(V),   V   = spm_vol(V); end;                                                                                                                

def_flags = struct('interp',1,'vox',NaN,'bb',NaN,'wrap',[0 0 0],'preserve',0);
[def_flags.bb, def_flags.vox] = bbvox_from_V(prm.VG);

if nargin < 3,
	flags = def_flags;
else,
	fnms = fieldnames(def_flags);
	for i=1:length(fnms),
		if ~isfield(flags,fnms{i}),
			flags = setfield(flags,fnms{i},getfield(def_flags,fnms{i}));
		end;
	end;
end; 

if ~all(finite(flags.vox(:))), flags.vox = def_flags.vox; end;
if ~all(finite(flags.bb(:))),  flags.bb  = def_flags.bb;  end;

bb  = sort(flags.bb);
vox = abs(flags.vox);


% The old voxel size and origin notation is used here.
% This requires that the position and orientation
% of the template is transverse.  It would not be
% straitforward to account for templates that are
% in different orientations because the basis functions
% would no longer be seperable.  The seperable basis
% functions mean that computing the deformation field
% from the parameters is much faster.

% Adjust bounding box slightly - so it rounds to closest voxel.
bb(:,1) = round(bb(:,1)/vox(1))*vox(1);
bb(:,2) = round(bb(:,2)/vox(2))*vox(2);
bb(:,3) = round(bb(:,3)/vox(3))*vox(3);

M   = prm.VG(1).mat;
vxg = sqrt(sum(M(1:3,1:3).^2));
%if det(M(1:3,1:3))<0, vxg(1) = -vxg(1); end;
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
Dim = [length(x) length(y) length(z)];

X = x'*ones(1,Dim(2));
Y = ones(Dim(1),1)*y;

Tr = prm.Tr;
if (prod(size(Tr)) == 0),
	affine_only = 1;
	basX = 0; tx = 0;
	basY = 0; ty = 0;
	basZ = 0; tz = 0;
else,
	affine_only = 0;
	basX = spm_dctmtx(prm.VG(1).dim(1),size(Tr,1),x-1);
	basY = spm_dctmtx(prm.VG(1).dim(2),size(Tr,2),y-1);
	basZ = spm_dctmtx(prm.VG(1).dim(3),size(Tr,3),z-1);
	if flags.preserve,
		DX = spm_dctmtx(prm.VG(1).dim(1),size(Tr,1),x-1,'diff');
		DY = spm_dctmtx(prm.VG(1).dim(2),size(Tr,2),y-1,'diff');
		DZ = spm_dctmtx(prm.VG(1).dim(3),size(Tr,3),z-1,'diff');
	end;
end;


msk = cell(length(z),1);
t1 = cat(3,V.mat);
t2 = cat(1,V.dim);
t  = [reshape(t1,[16 length(V)])' t2(:,1:3)];

if prod(size(V))>1 & any(any(diff(t,1,1))),
	spm_progress_bar('Init',length(z),'Computing available voxels','planes completed');
	for j=1:length(z),   % Cycle over planes
		Count = zeros(Dim(1),Dim(2));
		if affine_only,
			% Generate a mask for where there is data for all images
			%----------------------------------------------------------------------------
			for i=1:prod(size(V)),
				[X2,Y2,Z2] = mmult(X,Y,z(j),V(i).mat\prm.VF.mat*prm.Affine);
				Count      = Count + getmask(X2,Y2,Z2,V(i).dim(1:3),flags.wrap);
			end;
		else,
			% Nonlinear deformations
			%----------------------------------------------------------------------------
			tx = reshape(reshape(Tr(:,:,:,1),size(Tr,1)*size(Tr,2),size(Tr,3)) *basZ(j,:)', size(Tr,1), size(Tr,2) );
			ty = reshape(reshape(Tr(:,:,:,2),size(Tr,1)*size(Tr,2),size(Tr,3)) *basZ(j,:)', size(Tr,1), size(Tr,2) );
			tz = reshape(reshape(Tr(:,:,:,3),size(Tr,1)*size(Tr,2),size(Tr,3)) *basZ(j,:)', size(Tr,1), size(Tr,2) );
			X1 = X    + basX*tx*basY';
			Y1 = Y    + basX*ty*basY';
			Z1 = z(j) + basX*tz*basY';

			% Generate a mask for where there is data for all images
			%----------------------------------------------------------------------------
			for i=1:prod(size(V)),
				[X2,Y2,Z2] = mmult(X1,Y1,Z1,V(i).mat\prm.VF.mat*prm.Affine);
				Count      = Count + getmask(X2,Y2,Z2,V(i).dim(1:3),flags.wrap);
			end;
		end;
		msk{j} = uint32(find(Count ~= prod(size(V))));
		spm_progress_bar('Set',j);
	end;
else,
	for j=1:length(z), msk{j} = uint32([]); end;
end;

% Loop over volumes.  This is an attempt to improve the efficiency of resampling the
% data, since the memory management should be able to keep the whole volume paged
% into RAM.  The only disadvantage is that the nonlinear deformations need to be
% computed over and over again - but this is pretty fast.  In theory, the deformations
% could be stored for a whole volume, but this would require quite a lot of memory.

d     = [flags.interp*[1 1 1]' flags.wrap(:)];

spm_progress_bar('Init',prod(size(V)),'Resampling','volumes completed');
for i=1:prod(size(V)),
	VO          = V(i);
	VO.fname    = prepend(V(i).fname,'n');
	VO.mat      = prm.VG.mat*inv(M1)*M2;
	if spm_flip_analyze_images,
		Flp    = [-1 0 0 (Dim(1)+1); 0 1 0 0; 0 0 1 0; 0 0 0 1];
		VO.mat = VO.mat*Flp;
	end;
	VO.dim(1:3) = Dim;
	VO.descrip  = ['spm - 3D normalized'];

	detAff      = det(prm.VF.mat*prm.Affine/prm.VG.mat);
	if flags.preserve & affine_only,
		VO.pinfo(1:2,:) = VO.pinfo(1:2,:)/detAff;
	end;
	if flags.preserve & ~affine_only,
		Dat = zeros(VO.dim(1:3));
	else,
		VO  = spm_create_vol(VO);
	end;

	C = spm_bsplinc(V(i),d);

	if affine_only,
		for j=1:length(z),   % Cycle over planes
			[X2,Y2,Z2]  = mmult(X,Y,z(j),V(i).mat\prm.VF.mat*prm.Affine);
			dat         = spm_bsplins(C,X2,Y2,Z2,d);
			if flags.preserve, dat = dat*detAff; end;
			dat(msk{j}) = NaN;
			if spm_flip_analyze_images, dat = flipud(dat); end;
			spm_write_plane(VO,dat,j);
			if prod(size(V))<5, spm_progress_bar('Set',i-1+j/length(z)); end;
		end;
	else, % if ~affine_only,
		for j=1:length(z),   % Cycle over planes
			% Nonlinear deformations
			%----------------------------------------------------------------------------
			tx = reshape(reshape(Tr(:,:,:,1),size(Tr,1)*size(Tr,2),size(Tr,3)) *basZ(j,:)', size(Tr,1), size(Tr,2) );
			ty = reshape(reshape(Tr(:,:,:,2),size(Tr,1)*size(Tr,2),size(Tr,3)) *basZ(j,:)', size(Tr,1), size(Tr,2) );
			tz = reshape(reshape(Tr(:,:,:,3),size(Tr,1)*size(Tr,2),size(Tr,3)) *basZ(j,:)', size(Tr,1), size(Tr,2) );
			X1 = X    + basX*tx*basY';
			Y1 = Y    + basX*ty*basY';
			Z1 = z(j) + basX*tz*basY';

			[X2,Y2,Z2]  = mmult(X1,Y1,Z1,V(i).mat\prm.VF.mat*prm.Affine);
			dat         = spm_bsplins(C,X2,Y2,Z2,d);
			dat(msk{j}) = NaN;
			if ~flags.preserve,
				if spm_flip_analyze_images, dat = flipud(dat); end;
				spm_write_plane(VO,dat,j);
			else,
				tdx = reshape(reshape(Tr(:,:,:,1),size(Tr,1)*size(Tr,2),size(Tr,3)) *DZ(j,:)', size(Tr,1), size(Tr,2) );
				tdy = reshape(reshape(Tr(:,:,:,2),size(Tr,1)*size(Tr,2),size(Tr,3)) *DZ(j,:)', size(Tr,1), size(Tr,2) );
				tdz = reshape(reshape(Tr(:,:,:,3),size(Tr,1)*size(Tr,2),size(Tr,3)) *DZ(j,:)', size(Tr,1), size(Tr,2) );

				j11 = DX*tx*basY' + 1; j12 = basX*tx*DY';     j13 = basX*tdx*basY';
				j21 = DX*ty*basY';     j22 = basX*ty*DY' + 1; j23 = basX*tdy*basY';
				j31 = DX*tz*basY';     j32 = basX*tz*DY';     j33 = basX*tdz*basY' + 1;

				% The determinant of the Jacobian reflects relative volume changes.
				%------------------------------------------------------------------
				detJ       = detAff * (j11.*(j22.*j33 - j23.*j32) - j21.*(j12.*j33 - j13.*j32) + j31.*(j12.*j23 - j13.*j22));
				dat        = dat.*detJ;
				if spm_flip_analyze_images, dat = flipud(dat); end;
				Dat(:,:,j) = dat;
			end;
			if prod(size(V))<5, spm_progress_bar('Set',i-1+j/length(z)); end;
		end;
	end;
	if flags.preserve & ~affine_only,
		VO = spm_write_vol(VO,Dat);
	else,
		VO = spm_close_vol(VO);
	end;
	spm_progress_bar('Set',i);
end;
spm_progress_bar('Clear');
return;
%_______________________________________________________________________

%_______________________________________________________________________
function PO = prepend(PI,pre)
[pth,nm,xt,vr] = fileparts(deblank(PI));
PO             = fullfile(pth,[pre nm xt vr]);
return;
%_______________________________________________________________________

%_______________________________________________________________________
function Mask = getmask(X,Y,Z,dim,wrp)
% Find range of slice
tiny = 5e-2;
Mask = logical(ones(size(X)));
if ~wrp(1), Mask = Mask & (X >= (1-tiny) & X <= (dim(1)+tiny)); end;
if ~wrp(2), Mask = Mask & (Y >= (1-tiny) & Y <= (dim(2)+tiny)); end;
if ~wrp(3), Mask = Mask & (Z >= (1-tiny) & Z <= (dim(3)+tiny)); end;
return;
%_______________________________________________________________________

%_______________________________________________________________________
function [X2,Y2,Z2] = mmult(X1,Y1,Z1,Mult);
if length(Z1) == 1,
	X2= Mult(1,1)*X1 + Mult(1,2)*Y1 + (Mult(1,3)*Z1 + Mult(1,4));
	Y2= Mult(2,1)*X1 + Mult(2,2)*Y1 + (Mult(2,3)*Z1 + Mult(2,4));
	Z2= Mult(3,1)*X1 + Mult(3,2)*Y1 + (Mult(3,3)*Z1 + Mult(3,4));
else,
	X2= Mult(1,1)*X1 + Mult(1,2)*Y1 + Mult(1,3)*Z1 + Mult(1,4);
	Y2= Mult(2,1)*X1 + Mult(2,2)*Y1 + Mult(2,3)*Z1 + Mult(2,4);
	Z2= Mult(3,1)*X1 + Mult(3,2)*Y1 + Mult(3,3)*Z1 + Mult(3,4);
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


function TabDat = spm_VOI(SPM,hReg)
% List of local maxima and adjusted p-values for a small Volume of Interest
% FORMAT TabDat = spm_VOI(SPM,hReg)
%
% SPM    - structure containing SPM, distribution & filtering details
%        - required fields are:
% .swd   - SPM working directory - directory containing current SPM.mat
% .Z     - minimum of n Statistics {filtered on u and k}
% .n     - number of conjoint tests
% .STAT  - distribution {Z, T, X or F}
% .df    - degrees of freedom [df{interest}, df{residual}]
% .u     - height threshold
% .k     - extent threshold {resels}
% .XYZ   - location of voxels {voxel coords}
% .XYZmm - location of voxels {mm}
% .S     - search Volume {voxels}
% .R     - search Volume {resels}
% .FWHM  - smoothness {voxels}
% .M     - voxels - > mm matrix
% .VOX   - voxel dimensions {mm}
% .DIM   - image dimensions {voxels} - column vector
% .Vspm  - Mapped statistic image(s)
% .Msk   - mask: a list of scalar indicies into image voxel space
%
% hReg   - Handle of results section XYZ registry (see spm_results_ui.m)
%
% TabDat - Structure containing table data
%        - see spm_list for definition
%
%_______________________________________________________________________
%
% spm_VOI is  called by the SPM results section and takes variables in
% SPM to compute p-values corrected for a specified volume of interest.
%
% The volume of interest may be defined as a box or sphere centred on
% the current voxel, by the cluster in which it is embedded, or by an
% external mask image.
%
% If the VOI is defined by the cluster in which it is embedded, this
% cluster must have been defined independently of the SPM using a mask
% based on an orthogonal contrast and u = -Inf (i.e. p = 1)
%
% External mask images should be in the same orientation as the SPM
% (i.e. as the input used in stats estimation). The VOI is defined by
% voxels with values greater than 0.
%
% See also: spm_list
%_______________________________________________________________________
% %W% Karl Friston %E%

%-Parse arguments
%-----------------------------------------------------------------------
if nargin < 1,   error('insufficient arguments'), end
if nargin < 2,	 hReg = []; end

Num     = 16;			% maxima per cluster
Dis     = 04;			% distance among maxima (mm)

%-Title
%-----------------------------------------------------------------------
spm('FigName',['SPM{',SPM.STAT,'}: Small Volume Correction']);

%-Get current location
%-----------------------------------------------------------------------
xyzmm   = spm_results_ui('GetCoords');
v2mm    = diag([VOL.M(1,1) VOL.M(2,2) VOL.M(3,3)]);
mm2v    = inv(v2mm);
xyz     = VOL.M \ [xyzmm; 1]; xyz = xyz(1:3);

%-Specify search volume
%-----------------------------------------------------------------------
str     = sprintf(' at [%.0f,%.0f,%.0f]',xyzmm(1),xyzmm(2),xyzmm(3));
SPACE   = spm_input('Search volume...',-1,'m',...
		{['Sphere',str],['Box',str],'Nearest cluster',...
		'Image'},['S','B','V','I']);
Q       = ones(1,size(SPM.XYZmm,2));
vsc     = [1 1 1];				%-Voxel size scaling for FWHM
DIM     = VOL.DIM;

switch SPACE, case 'S'                                          % Sphere
	%---------------------------------------------------------------
	D     = spm_input('radius of spherical VOI {mm}',-2);
	str   = sprintf('%0.1fmm sphere',D);
	j     = find(sum((SPM.XYZmm - xyzmm*Q).^2) <= D^2);

	d     = ceil(mm2v*[D D D]')+1;
	rxyz  = round(xyz);
	[jx jy jz]  = ndgrid(rxyz(1)+(-d(1):d(1)),...
			     rxyz(2)+(-d(2):d(2)),...
			     rxyz(3)+(-d(3):d(3)));
	XYZ   = [jx(:) jy(:) jz(:)]';
	XYZ(:,any(XYZ<=0)) = []; 
	XYZ(:,(XYZ(1,:)>DIM(1)|XYZ(2,:)>DIM(2)|XYZ(3,:)>DIM(3))) = []; 
	Qs    = ones(1,size(XYZ,2));
	js    = find(sum((v2mm*(XYZ - xyz*Qs)).^2) <= D^2);
	XYZ   = XYZ(:,js);

	% Convert to indicies; map T image, check image for zeros; delete
        % those entries from XYZ

	D     = D./SPM.VOX;
	S     = (4/3)*pi*prod(D);

case 'B'                                                           % Box
	%---------------------------------------------------------------
	D     = spm_input('box dimensions [k l m] {mm}',-2);
	str   = sprintf('%0.1f x %0.1f x %0.1f mm box',D(1),D(2),D(3));
	j     = find(all(abs(SPM.XYZmm - xyzmm*Q) <= D(:)*Q/2));

	d     = ceil(mm2v*D(:))+1;
	rxyz  = round(xyz);
	[jx jy jz]  = ndgrid(rxyz(1)+(-d(1):d(1)),...
			     rxyz(2)+(-d(2):d(2)),...
			     rxyz(3)+(-d(3):d(3)));
	XYZ   = [jx(:) jy(:) jz(:)]';
	XYZ(:,any(XYZ<=0)) = [];
	XYZ(:,(XYZ(1,:)>DIM(1)|XYZ(2,:)>DIM(2)|XYZ(3,:)>DIM(3))) = []; 
	Qs    = ones(1,size(XYZ,2));
	js    = find(all(abs(v2mm*(XYZ - xyz*Qs)) <= D(:)*Qs/2));
	XYZ   = XYZ(:,js);

	% Convert to indicies; map T image, check image for zeros; delete
        % those entries from XYZ

	D     = D(:)./SPM.VOX(:);
	S     = prod(D);

case 'V'                                                        %-Voxels
	%---------------------------------------------------------------
	if ~length(SPM.XYZ)
		spm('alert!','No suprathreshold clusters!',mfilename,0);
		spm('FigName',['SPM{',SPM.STAT,'}: Results']);
		return
	end

	[xyzmm,i] = spm_XYZreg('NearestXYZ',xyzmm,SPM.XYZmm);
	spm_results_ui('SetCoords',xyzmm);
	A     = spm_clusters(SPM.XYZ);
	j     = find(A == A(i));
	str   = sprintf('%0.0f voxel cluster',length(j));
	D     = SPM.XYZ(:,j);
	S     = length(j);

case 'I'                                                         % Image
	%---------------------------------------------------------------
	%-The VOI image must be in the same (real) space as the SPM:
	%-Although the masking code below does take account of the
	% .mat file for the VOI image, the orientation of the VOI image
	% .mat file is *not* used in the resel calculation,
	% (see spm_resels_vol.c) so any mat file should contain
	% no rotations / shears relative to the SPM, for accuracy.
	im    = spm_get(1,'img','Image defining volume subset');
	D     = spm_vol(im);
	tM    = D.mat \ SPM.M;
	if any(tM([2:5,7:10,12]))
		spm('alert!',{	'Mask image rotated/sheared!',...
				'(relative to SPM image)',...
				'Can''t use for SVC.'},mfilename,0);
		spm('FigName',['SPM{',SPM.STAT,'}: Results']);
		return
	end
	mXYZ  = D.mat \ [SPM.XYZmm; ones(1, size(SPM.XYZmm, 2))];
	j     = find(spm_sample_vol(D, mXYZ(1,:),mXYZ(2,:),mXYZ(3,:),0) > 0);
	str   = sprintf('image mask: %s',spm_str_manip(im,'a30'));
	%-Compute in-mask volume S:
	% Correct for differences in mask and SPM voxel sizes
	Y     = spm_read_vols(D);
	vsc   = sqrt(sum(D.mat(1:3,1:3).^2)) ./ SPM.VOX;
	S     = sum(Y(:)>0) * prod(vsc);

end

spm('Pointer','Watch')

%-Select voxels within subspace
%-----------------------------------------------------------------------
SPM.Z     = SPM.Z(j);
SPM.XYZ   = SPM.XYZ(:,j);
SPM.XYZmm = SPM.XYZmm(:,j);
SPM.R     = spm_resels(SPM.FWHM./vsc,D,SPACE);
SPM.S     = S;
if (SPACE=='S') | (SPACE=='B')
  SPM.Msk   = XYZ(1,:) + ...
             (XYZ(2,:)-1)*DIM(1) + ...
             (XYZ(3,:)-1)*DIM(1)*DIM(2);
  % Eliminate voxels with no data; assume zero masking of stat images
  Y         = spm_read_vols(SPM.Vspm);
  SPM.Msk(Y(SPM.Msk)==0) = [];
  % Note, we get a more accurate (smaller) S, but the RFT doesn't
  % make use of it.
  SPM.S = length(SPM.Msk);
elseif (SPACE=='V')
  SPM.Msk   = SPM.XYZ(1,:) + ...
             (SPM.XYZ(2,:)-1)*DIM(1) + ...
             (SPM.XYZ(3,:)-1)*DIM(1)*DIM(2);
elseif (SPACE=='I')
  SPM.Msk   = D;
end


%-Tabulate p values
%-----------------------------------------------------------------------
str    = sprintf('search volume: %s',str);
if any(strcmp(SPACE,{'S','B','V'}))
	str = sprintf('%s at [%.0f,%.0f,%.0f]',str,xyzmm(1),xyzmm(2),xyzmm(3));
end
TabDat = spm_list('List',SPM,hReg);

%-Reset title
%-----------------------------------------------------------------------
spm('FigName',['SPM{',SPM.STAT,'}: Results']);

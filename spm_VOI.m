function TabDat = spm_VOI(SPM,VOL,Dis,Num,hReg)
% List of local maxima and adjusted p-values for a small Volume of Interest
% FORMAT TabDat = spm_VOI(SPM,VOL,Dis,Num,hReg)
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
%
% VOL    - structure containing details of volume analysed
%        - required fields are:
% .S     - search Volume {voxels}
% .R     - search Volume {resels}
% .FWHM  - smoothness {voxels}
% .iM    - mm -> voxels matrix
% .VOX   - voxel dimensions {mm}
%
% Dis    - Minimum distance between maxima
%          (Passed to spm_list.m: Defaults on missing or empty)
% Num    - Maxiumum number of local maxima tabulated per cluster
%          (Passed to spm_list.m: Defaults on missing or empty)
% hReg   - Handle of results section XYZ registry (see spm_results_ui.m)
%
% TabDat - Structure containing table data
%        - see spm_list for definition
%
%_______________________________________________________________________
%
% spm_VOI is  called by the SPM results section and takes variables in
% SPM and VOL to compute p-values corrected for a specified volume of
% interest.
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
if nargin<2,     error('insufficient arguments'), end
if nargin<5,	 hReg = []; end
if nargin<4,     Num  = []; end
if isempty(Num), Num  = 16; end
if nargin<3,     Dis  = []; end
if isempty(Dis), Dis  = 04; end

%-Title
%-----------------------------------------------------------------------
spm('FigName',['SPM{',SPM.STAT,'}: Small Volume Correction']);

%-Get current location
%-----------------------------------------------------------------------
xyzmm   = spm_results_ui('GetCoords');


%-Specify search volume
%-----------------------------------------------------------------------
str     = sprintf(' at [%.0f,%.0f,%.0f]',xyzmm(1),xyzmm(2),xyzmm(3));
SPACE   = spm_input('Search volume...',-1,'m',...
		{['Sphere',str],['Box',str],'Nearest cluster',...
		'Image'},['S','B','V','I']);
Q       = ones(1,size(SPM.XYZmm,2));

switch SPACE, case 'S'                                          % Sphere
	%---------------------------------------------------------------
	D     = spm_input('radius of spherical VOI {mm}',-2);
	str   = sprintf('%0.1fmm sphere',D);
	j     = find(sum((SPM.XYZmm - xyzmm*Q).^2) <= D^2);
	D     = D./VOL.VOX;
	S     = (4/3)*pi*prod(D);

case 'B'                                                           % Box
	%---------------------------------------------------------------
	D     = spm_input('box dimensions [k l m] {mm}',-2);
	str   = sprintf('%0.1f x %0.1f x %0.1f mm box',D(1),D(2),D(3));
	j     = find(all(abs(SPM.XYZmm - xyzmm*Q) <= D(:)*Q/2));
	D     = D(:)./VOL.VOX;
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
	%-The image must be in the same (real) space as the SPM:
	%-Although the voxel thresholding does take account of the
	% .mat file for the VOI image, the orientation of the VOI image
	% .mat file is *not* used in the resel calculation,
	% (see spm_resels_vol.c) so any mat file should contain
	% no rotations / shears, for accuracy.
	im    = spm_get(1,'img','Image defining volume subset');
	D     = spm_vol(im);
	if any(D.mat([2:5,7:10,12]))
		spm('alert!',{	'Mask image rotated/sheared!',...
				'Can''t use for SVC.'},mfilename,0);
		spm('FigName',['SPM{',SPM.STAT,'}: Results']);
		return
	end
	mXYZ  = D.mat \ [SPM.XYZmm; ones(1, size(SPM.XYZmm, 2))];
	j     = spm_sample_vol(D, mXYZ(1,:),mXYZ(2,:),mXYZ(3,:),0) > 0;	
	str   = sprintf('image mask: %s',spm_str_manip(im,'a30'));
	%-Compute in-mask volume S:
	% Correct for differences in mask and SPM voxel sizes
	Y     = spm_read_vols(D);
	vsc   = prod(sqrt(sum(D.mat(1:3,1:3).^2))) / prod(VOL.VOX);
	S     = sum(Y(:)>0) * vsc;

end

spm('Pointer','Watch')

%-Select voxels within subspace
%-----------------------------------------------------------------------
SPM.Z     = SPM.Z(j);
SPM.XYZ   = SPM.XYZ(:,j);
SPM.XYZmm = SPM.XYZmm(:,j);
VOL.R     = spm_resels(VOL.FWHM,D,SPACE);
VOL.S     = S;


%-Tabulate p values
%-----------------------------------------------------------------------
str = sprintf('search volume: %s',str);
if any(strcmp(SPACE,{'S','B','V'}))
	str = sprintf('%s at [%.0f,%.0f,%.0f]',str,xyzmm(1),xyzmm(2),xyzmm(3));
end
TabDat = spm_list('List',SPM,VOL,Dis,Num,str,hReg);

%-Reset title
%-----------------------------------------------------------------------
spm('FigName',['SPM{',SPM.STAT,'}: Results']);

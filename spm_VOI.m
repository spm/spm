function TabDat = spm_VOI(SPM,VOL)
% List of local maxima and adjusted p-values for a small Volume of Interest
% FORMAT TabDat = spm_VOI(SPM,VOL)
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
% TabDat - Structure containing table data
%        - see spm_list for definition
%
%_______________________________________________________________________
%
% spm_VOI is  called by the SPM results section and takes variables in
% SPM and VOL to compute p-values corrected for a specified volume of
% interest centred on the current voxel.
%
% This volume may be defined by the cluster in which it is embedded. Clearly
% this cluster must have been defined independently of the SPM using a
% mask based on an orthogonal contrast and u = -Inf (i.e. p = 1)
%
% See also: spm_list
%_______________________________________________________________________
% %W% Karl Friston %E%

%-Title
%-----------------------------------------------------------------------
spm('FigName',['SPM{',SPM.STAT,'}: Small Volume Correction']);

%-Get current location
%-----------------------------------------------------------------------
xyzmm   = spm_results_ui('GetCoords');


%-Specify search volume
%-----------------------------------------------------------------------
SPACE   = spm_input('search volume',-1,'b','Sphere|Box|Cluster',['S','B','V']);
Q       = ones(1,size(SPM.XYZmm,2));
if     SPACE == 'S'

	D     = spm_input('radius of spherical VOI {mm}',-2);
	str   = sprintf('%0.1fmm sphere',D);
	j     = find(sum((SPM.XYZmm - xyzmm*Q).^2) <= D^2);
	D     = D./VOL.VOX;
	S     = (4/3)*pi*prod(D);

elseif SPACE == 'B'

	D     = spm_input('box dimensions [k l m] {mm}',-2);
	str   = sprintf('%0.1f x %0.1f x %0.1f mm box',D(1),D(2),D(3));
	j     = find(all(abs(SPM.XYZmm - xyzmm*Q) <= D(:)*Q/2));
	D     = D(:)./VOL.VOX;
	S     = prod(D);

elseif SPACE == 'V'
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
str = sprintf('search volume: %s at (%.0f,%.0f,%.0f)',str,...
	xyzmm(1),xyzmm(2),xyzmm(3));
TabDat = spm_list('List',SPM,VOL,4,16,str);

%-Reset title
%-----------------------------------------------------------------------
spm('FigName',['SPM{',SPM.STAT,'}: Results']);

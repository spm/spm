function spm_VOI(SPM,VOL,hReg)
% Tabular display of adjusted data
% FORMAT spm_VOI
%
% SPM  - SPM structure      {'Z' 'n' 'STAT' 'df' 'u' 'k'}
% VOL  - Spatial structure  {'R' 'FWHM' 'S' 'DIM' 'VOX' 'ORG' 'M' 'XYZ' 'QQ'}
% hReg - handle of MIP register
%
%_______________________________________________________________________
%
% spm_VOI is  called by spm_results and takes variables in SPM and VOL to
% compute a p values corrected for a specified volume of interest that is
% centred on the current voxel.
%
% This volume may be defined by the cluster in which it is embedded. Clearly
% this cluster must have been defined independently of the SPM using a
% mask based on an orthogonal contrast and u = -Inf (i.e. p = 1)
%
%_______________________________________________________________________
% %W% Karl Friston %E%

% Title
%-----------------------------------------------------------------------
Finter  = spm_figure('FindWin','Interactive');
set(Finter,'Name','Small Volume Correction');

% Identify cluster & update GUI
%-----------------------------------------------------------------------
L       = spm_XYZreg('GetCoords',hReg);
[L,i]   = spm_XYZreg('NearestXYZ',L,VOL.XYZ);

spm_XYZreg('SetCoords',L,hReg);


% Specify search volume
%-----------------------------------------------------------------------
SPACE   = spm_input('search volume',1,'b','Sphere|Box|Cluster',['S' 'B' 'V']);
Q       = ones(1,size(VOL.XYZ,2));
L       = L*Q;
if     SPACE == 'S'

	D     = spm_input('radius of spherical VOI {mm}',2);
	str   = sprintf('%0.1f mm sphere',D);
	j     = find(sum((VOL.XYZ - L).^2) <= D^2);
	D     = D./VOL.VOX;
	S     = (4/3)*pi*prod(D);

elseif SPACE == 'B'

	D     = spm_input('box dimensions [k l m] {mm}',2);
	str   = sprintf('%0.1f x %0.1f x %0.1f mm box',D(1),D(2),D(3));
	j     = find(all(abs(VOL.XYZ - L) <= D(:)*Q/2));
	D     = D(:)./VOL.VOX;
	S     = prod(D);

elseif SPACE == 'V'

	A     = spm_clusters(VOL.XYZ,VOL.VOX);
	j     = find(A == A(i));
	str   = sprintf('%0.0f voxel cluster',length(j));
	D     = VOL.XYZ(:,j); D  = D./(VOL.VOX*ones(1,size(D,2)));
	S     = length(j);
end


% Select voxels within subspace
%-----------------------------------------------------------------------
SPM.Z   = SPM.Z(j);
VOL.XYZ = VOL.XYZ(:,j);
VOL.R   = spm_resels(VOL.FWHM,D,SPACE);
VOL.S   = S;


% Tabulate p values
%-----------------------------------------------------------------------
spm_list(SPM,VOL,4,16)

title(['Search volume: ' str])

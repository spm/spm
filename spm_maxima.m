function spm_maxima(SPM,VOL,hReg)
% Tabular display of maxima in selected cluster
% FORMAT spm_maxima(SPM,VOL,hReg)
%
% SPM  - SPM structure      {'Z' 'n' 'STAT' 'df' 'u' 'k'}
% VOL  - Spatial structure  {'R' 'FWHM' 'S' 'DIM' 'VOX' 'ORG' 'M' 'XYZ' 'QQ'}
% hReg - handle of MIP register
%_______________________________________________________________________
%
% spm_maxima is called by spm_results and takes variables in the
% structures SPM and VOL to produce a table of maxima within the selected
% region.  The region and its associated maxima are characterized in
% terms of its cluster, and voxel-level p values (corrected and
% uncorrected).
%
% The maxima displayed are all at least 4mm apart. Selecting the voxel 
% coordinates of a maxima causes a green pointer to appear in the 
% appropriate place on maximum intensity projection.
%
% see also spm_list
%
%_______________________________________________________________________
% %W% Karl Friston, Andrew Holmes %E%

% Find the cluster selected
%-----------------------------------------------------------------------
L       = spm_XYZreg('GetCoords',hReg);
[L,i]   = spm_XYZreg('NearestXYZ',L,VOL.XYZ);


% Find selected cluster
%-----------------------------------------------------------------------
A       = spm_clusters(VOL.XYZ,VOL.VOX);
j       = find(A == A(i));
SPM.Z   = SPM.Z(j);
VOL.XYZ = VOL.XYZ(:,j);


% Update GUI to cluster's maximum
%-----------------------------------------------------------------------
[i j]   = max(SPM.Z);
spm_XYZreg('SetCoords',VOL.XYZ(:,j),hReg);


% Tabulate p values
%-----------------------------------------------------------------------
spm_list(SPM,VOL,4,16)

title('Maxima within cluster')

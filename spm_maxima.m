function spm_maxima(SPM,VOL,hReg)
% Tabular display of maxima in selected cluster
% FORMAT spm_maxima(SPM,VOL,hReg)
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
% .XYZ   - location of voxels {mm}
%
% VOL    - structure containing details of volume analysed
%        - required fields are:
% .S     - search Volume {voxels}
% .R     - search Volume {resels}
% .FWHM  - smoothness {voxels}     
% .iM    - mm -> voxels matrix
% .VOX   - voxel dimensions {mm}
%
% hReg   - handle of MIP XYZ registry object (see spm_XYZreg for details)
%_______________________________________________________________________
%
% spm_maxima is called by the SPM results section, and takes variables
% in the structures SPM and VOL to produce a table of maxima within the
% selected region.  The region and its associated maxima are
% characterized in terms of its cluster, and voxel-level p values
% (corrected and uncorrected).
%
% The maxima displayed are all at least 4mm apart. Selecting the voxel
% coordinates of a maxima causes a green pointer to appear in the
% appropriate place on maximum intensity projection.
%
% See also: spm_list
%
%_______________________________________________________________________
% %W% Karl Friston, Andrew Holmes %E%

%-Find the cluster selected
%-----------------------------------------------------------------------
xyz     = spm_XYZreg('GetCoords',hReg);
[xyz,i] = spm_XYZreg('NearestXYZ',xyz,SPM.XYZ);


%-Find selected cluster
%-----------------------------------------------------------------------
rcp     = VOL.iM(1:3,:)*[SPM.XYZ; ones(1,size(SPM.XYZ,2))];
A       = spm_clusters(rcp,[1,1,1]);
j       = find(A == A(i));
SPM.Z   = SPM.Z(j);
SPM.XYZ = SPM.XYZ(:,j);


%-Update GUI to cluster's maximum
%-----------------------------------------------------------------------
[i,j]   = max(SPM.Z);
spm_XYZreg('SetCoords',SPM.XYZ(:,j),hReg);


%-Tabulate p values
%-----------------------------------------------------------------------
str = 'single cluster list (corrected p-values for entire volume)';
spm_list(SPM,VOL,4,16,str)

title('Maxima within cluster')

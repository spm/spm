function spm_sections(SPM,VOL,hReg)
% rendering of regional effects [SPM{Z}] on orthogonal sections
% FORMAT spm_sections(SPM,VOL,hReg)
%
% SPM  - SPM structure      {'Z' 'n' 'STAT' 'df' 'u' 'k'}
% VOL  - Spatial structure  {'R' 'FWHM' 'S' 'DIM' 'VOX' 'ORG' 'M' 'XYZ' 'QQ'}
% hReg - handle of MIP register
%
% see spm_getSPM for details
%_______________________________________________________________________
%
% spm_sections is called by spm_results and uses variables in SPM and
% VOL to create three orthogonal sections though a background image.
% Regional foci from the selected SPM are rendered on this image.
%
%_______________________________________________________________________
% %W%	John Ashburner %E%

Fgraph = spm_figure('FindWin','Graphics');
spms   = spm_get(1,'.img','select an image for rendering');
spm_results_ui('ClearPane',Fgraph,'RNP');
spm_orthviews('Image',spms,[0.05 0.05 0.9 0.4]);
spm_orthviews('register',hReg);
spm_orthviews('addblobs',1,VOL.XYZ,SPM.Z,VOL.M);
spm_orthviews('Redraw');

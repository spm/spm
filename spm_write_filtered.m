function spm_write_filtered(SPM,VOL)
% Writes the filtered SPM as an image
% FORMAT spm_write_filtered(SPM,VOL)
%
% SPM  - SPM structure      {'Z' 'n' 'STAT' 'df' 'u' 'k'}
% VOL  - Spatial structure  {'R' 'FWHM' 'S' 'DIM' 'VOX' 'ORG' 'M' 'XYZ' 'QQ'}
%
% see spm_getSPM for details
%
%-----------------------------------------------------------------------
% This function is intended to be called from the environment set up by
% spm_results_ui.
%
%_______________________________________________________________________
% @(#)spm_write_filtered.m	1.1 FIL 96/09/10

global CWD;

%-----------------------------------------------------------------------
Q       = spm_input('Output filename',1,'s');
Finter  = spm_figure('FindWin','Interactive');
ptr     = get(Finter,'Pointer');
set(Finter,'Pointer','watch');

% Set up header information
%-----------------------------------------------------------------------
q       = max([find(Q == '/') 0]);
Q       = [CWD '/' spm_str_manip(Q((q + 1):length(Q)),'sd')];
str     = sprintf('spm{%c}-filtered: u = %5.3f, k = %d',SPM.STAT,SPM.u,SPM.k);
V       = struct(...
	'fname',	Q,...
	'dim',		[VOL.DIM' spm_type('uint8')],...
	'mat',		VOL.M,...
	'descrip', 	str);

%-Reconstruct filtered image from XYZ & SPM.Z
%-----------------------------------------------------------------------
Y      = zeros(VOL.DIM(1:3)');
IM     = inv(VOL.M);
XYZ    = round(IM(1:3,:)*[VOL.XYZ ; ones(1,size(VOL.XYZ,2))]);
OFF    = XYZ(1,:) + VOL.DIM(1)*(XYZ(2,:) + VOL.DIM(2)*XYZ(3,:));
Y(OFF)  = SPM.Z.*(SPM.Z > 0);

% Write the filtered volume
%-----------------------------------------------------------------------
spm_write_vol(V,Y);

set(Finter,'Pointer',ptr);

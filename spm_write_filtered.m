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
q       = max([find(Q == '/') 0]);
Q       = [CWD '/' spm_str_manip(Q((q + 1):length(Q)),'sd')];
str     = sprintf('spm{%c}-filtered: u = %5.3f, k = %d',SPM.STAT,SPM.u,SPM.k);
set(Finter,'Pointer','watch');

%-Reconstruct filtered image from XYZ & SPM.Z
%-----------------------------------------------------------------------
XYZ     = VOL.XYZ;
DIM     = VOL.DIM;
VOX     = VOL.VOX;
ORG     = VOL.ORG;

n       = size(XYZ,2);
rcp     = round(XYZ./meshgrid([1;1;1].*VOX,1:n)' + meshgrid(ORG,1:n)');
Dim     = cumprod([1,DIM(1:2)']);
OffSets = meshgrid([0,1,1],1:n)';
i       = Dim*(rcp - OffSets);
Z       = SPM.Z.*(SPM.Z > 0);
T       = zeros(1,prod(DIM));
T(i)    = Z;
MAX     = max(max(T));
T       = round(T*(255/MAX));

%-Write out to analyze file
%-----------------------------------------------------------------------
fid     = fopen([Q,'.img'],'w');
if (fid == -1)
	set(Finter,'Pointer',ptr);
	error(['Failed to open ' Q '.img']);
end
if fwrite(fid,T,spm_type(2)) ~= prod(size(T))
	fclose(fid);
	set(Finter,'Pointer',ptr);
	error(['Failed to write ' Q '.img']);
end
fclose(fid);
spm_hwrite([Q,'.hdr'],DIM,VOX,MAX/255,2,0,ORG,str);

%-Finished
%-----------------------------------------------------------------------
set(Finter,'Pointer',ptr);

function spm_write_filtered(Z,XYZ,u,k,V,SPMZ)
% Writes out the filtered Z or F images
% FORMAT spm_write_filtered(Z,XYZ,u,k,V,SPMZ)
%  Z      - Z or F values after filtering on height and size
%           thresholds
%  XYZ    - location in mm
%  u      - selected height threshold
%  k      - selected extent threshold {voxels}
%  V      - contains image dimensions, voxel sizes and origin
%  SPMZ   - 1 for Z image, 0 for F image
%-----------------------------------------------------------------------
% This function is intended to be called from environment set up by
% spm_results_ui. 
%_______________________________________________________________________
% %W% FIL %E%

global CWD;
Q = spm_input('Output filename',1,'s');
w = spm_figure('FindWin','Interactive');
ptr = get(w,'Pointer');
set(w,'Pointer','watch');

q = max([find(Q == '/') 0]);
Q = [CWD '/' spm_str_manip(Q((q+1):length(Q)),'sd')];

if SPMZ == 1
	tmp = 'Z';
else
	tmp = 'F';
end
str = sprintf('spm{%c}-filtered: u = %5.3f, k = %d',tmp,u,k);

%-Reconstruct filtered image from XYZ & t
%---------------------------------------------------------------
n       = size(XYZ,2);
rcp     = round(XYZ./meshgrid([1;1;1].*V(4:6),1:n)' + ...
	  meshgrid(V(7:9),1:n)');
Dim     = cumprod([1,V(1:2)']);
OffSets = meshgrid([0,1,1],1:n)';
e       = ((rcp - OffSets)'*Dim')';
Z       = Z.*(Z > 0);
T       = zeros(1,prod(V(1:3)));
T(e)    = Z;
mx = max(max(T));
T = round(T*(255/mx));

%-Write out to analyze file
%---------------------------------------------------------------
fid     = fopen([Q,'.img'],'w');
if (fid == -1)
	set(w,'Pointer',ptr);
	error(['Failed to open ' Q '.img']);
end
if fwrite(fid,T,spm_type(2)) ~= prod(size(T))
	fclose(fid);
	set(w,'Pointer',ptr);
	error(['Failed to write ' Q '.img']);
end
fclose(fid);
spm_hwrite([Q,'.hdr'],V(1:3),V(4:6),mx/255,2,0,V(7:9),str);
set(w,'Pointer',ptr);

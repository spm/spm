function spm_bias_apply(P,T)
% Apply bias field to images.
%
% FORMAT spm_bias_apply(P,T)
%   P - filename of image
%   T - DCT of bias field or filename containing this
%
% Bias corrected images are written to disk, prefixed by 'm'.
%
%_______________________________________________________________________
% %W% John Ashburner %E%

V = spm_vol(P);
if ischar(T),
	s = load(T);
	T = s.T;
end;
nbas           = [size(T) 1];
nbas           = nbas(1:3);
B1             = spm_dctmtx(V(1).dim(1),nbas(1));
B2             = spm_dctmtx(V(1).dim(2),nbas(2));
B3             = spm_dctmtx(V(1).dim(3),nbas(3));

VO             = V;
VO.dim(4)      = spm_type('float');
VO.pinfo       = [1 0 0]';
[pth,nm,xt,vr] = fileparts(deblank(P));
VO.fname       = fullfile(pth,['m' nm xt vr]);
VO             = spm_create_vol(VO);

for p=1:V.dim(3),
	M   = spm_matrix([0 0 p]);
	img = spm_slice_vol(V, M, V.dim(1:2), 1);
	t   = reshape(T,  nbas(1)*nbas(2), nbas(3));
	t   = reshape(t*B3(p,:)', nbas(1), nbas(2));
	img = img.*exp(B1*t*B2');
	VO  = spm_write_plane(VO,img,p);
end;
VO = spm_close_vol(VO);
return;
%=======================================================================

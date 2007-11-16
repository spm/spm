function [Centre,c_mm] = spm_eeg_inv_CtrBin(P)

%=======================================================================
% FORMAT [Centre,c_mm]=spm_eeg_inv_CtrBin(P)
% Determine the "centre of mass" (vx and mm) of a binary image.
%=======================================================================

if nargin == 0
	P  = spm_select(1,'*.img','Binary image');
end

Vp     = spm_vol(P) ;

Vp.dat = spm_loaduint8(Vp);
X      = (1:Vp.dim(1))'*ones(1,Vp.dim(2)); X =X(:)';
Y      = ones(Vp.dim(1),1)*(1:Vp.dim(2)); Y = Y(:)';
Z      = zeros(Vp.dim(1),Vp.dim(2)); Z = Z(:)'; 
Unit   = ones(size(Z));
	
c_mm   = zeros(1,3); SI = 0;
for pp=1:Vp.dim(3)
     I_pl = double(Vp.dat(:,:,pp)>128); I_pl = I_pl(:)*ones(1,3);
     XYZ  = Vp.mat*[X ; Y ; Z+pp ; Unit];
     c_mm = c_mm + sum(XYZ(1:3,:)'.*I_pl);
     SI   = SI+sum(I_pl(:,1));
end
c_mm   = c_mm/SI;
Centre = spm_eeg_inv_mm2vx(c_mm,Vp.mat)';

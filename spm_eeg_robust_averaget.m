function [B,Wf]=spm_eeg_robust_averaget(data,ks,FS);

% function to apply robust averaging routine to data sets and return the
% ERP (B) and the weights (Wf)
% ks is the offest of the weighting function the default is 3.


% James Kilner
% $Id: spm_eeg_robust_averaget.m 340 2005-11-30 17:25:26Z james $
if nargin==1
	ks=3;
end
data=data';
Wf=zeros(size(data(:)));
ndata=reshape(data',size(data,1)*size(data,2),1);
Xs=sparse(repmat(speye(size(data,2)),[size(data,1),1]));
s1=size(data,1);
s2=size(data,2);
%clear data;
Wis=speye(length(ndata));
h=1./((1-(diag(Xs*(Xs'*Xs)^-1*Xs'))).^0.5);
ores=1;
nres=10;
n=0;
while abs(ores-nres)>sqrt(1E-8)
	abs(ores-nres);
	ores=nres;
	n=n+1;
	B=((Xs'*Wis*Xs)^-1)*Xs'*Wis*ndata;
	if sum(isnan(B))>0
		break
	end
	if n>100
		break
	end
	res=ndata-(Xs*B);
	

	mad=median(abs(res-median(res)));	
	res=(res)./mad;
	res=res.*h;	
 	sm=gausswin(FS);
 	sm=sm/sum(sm);
 	res=conv(sm,res);
	res=res(floor(FS/2):end-ceil(FS/2));
	res=abs(res)-ks;
	res(res<0)=0;	
	nres=(sum(res.^2));
	Wf=(((abs(res)<1) .* (1 - res.^2).^2));
	clear res;
	Wis=spdiags(Wf,0,Wis);
end





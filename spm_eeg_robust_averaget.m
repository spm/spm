function [B,Wf]=robust_averaget(data);

% James Kilner
% $Id: spm_eeg_robust_averaget.m 112 2005-05-04 18:20:52Z john $

data=data';
%figure(1)
%plot(mean(data));
% clear D
% clear d1

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
	if n>25
		break
	end
	res=ndata-Xs*B;
	mad=median(abs(res-median(res)));
	res=(res)./(mad);
	res=res.*h;	
	

	res=abs(res)-1;
	res(res<0)=0;	
	nres=(sum(res.^2));
	Wf=(((abs(res)<1) .* (1 - res.^2).^2));

	clear res;
	Wis=spdiags(Wf,0,Wis);
end





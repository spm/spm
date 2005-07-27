function [B,Wf]=robust_average(data);

% James Kilner
% $Id: spm_eeg_robust_average.m 204 2005-07-27 08:51:55Z james $

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
	if n>500
		break
	end
	res=ndata-Xs*B;
	mad=median(abs(res-median(res)));
	res=(res)./mad;
	res=res.*h;	
	res=reshape(res,s2,s1);
	res=repmat((mean(abs(res)).*sign(sum(res))),[s2,1]);
	res=reshape(res,s1*s2,1);
	res=abs(res)-2;
	res(res<0)=0;
	nres=(sum(res.^2));
	Wf=(((abs(res)<1) .* (1 - res.^2).^2));
	clear res;
	Wis=spdiags(Wf,0,Wis);


end



% 
% res=ndata-Xs*B;
% res=reshape(res,s2,s1);
% res=repmat(max(abs(res)),[s2,1]);
% res=reshape(res,s1*s2,1);
% mad=median(abs(res-median(res)));
% res=res./mad;
% res=res.*h;	
% res=abs(res)-1;
% res(res<0)=0;
% nres=(sum(res.^2));
% n=0;
% while abs(ores-nres)>sqrt(1E-8)
% 	n=n+1;
% 	figure(10)
% 	plot(n,abs(ores-nres),'.')
% 	hold on
% 	ores=nres;
% 	
% 	
% 	if n>50
% 		break
% 	end
% 	
% 	
% 	B=((Xs'*Wis*Xs)^-1)*Xs'*Wis*ndata;
% 	res=ndata-Xs*B;
% 	mad=median(abs(res-median(res)));
% 	res=(res)./mad;
% 	res=res.*h;	
% 	res=reshape(res,s2,s1);
% 	[df,inds]=(max(abs(res)));
% 	y=1:length(inds);
% 	res=repmat(diag(squeeze(res(inds,y)))',[s2,1]);
% 	res=reshape(res,s1*s2,1);
% 
% 	res=abs(res)-1;
% 	res(res<0)=0;
% 	nres=(sum(res.^2));
% 	Wf=(((abs(res)<1) .* (1 - res.^2).^2));
% 	clear res;
% 	Wis=spdiags(Wf,0,Wis);
% 
% 
% end
% 

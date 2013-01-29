function [B,Wf] = spm_eeg_robust_averaget(data,ks)
% Apply robust averaging routine to data sets 
% FORMAT [B,Wf] = spm_eeg_robust_averaget(data,ks)
% data   - data matrix to be averaged
% ks     - offest of the weighting function (default: 3)
%
% Wf     - estimated weights
%__________________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% James Kilner
% $Id: spm_eeg_robust_averaget.m 5219 2013-01-29 17:07:07Z spm $

if nargin==1
    ks = 3;
end

%Wf=ones(size(data));
try
    Xs = sparse(repmat(speye(size(data,2)),[size(data,1),1]));
    h = 1./((1-(diag(Xs*(Xs'*Xs)^-1*Xs'))).^0.5);
    h = h(1);
catch
    h = 1;
    warning('File too large: assuming h=1.');
end
ores = 1;
nres = 10;
n    = 0;
B    = zeros(1,size(data,1));
while abs(ores-nres)>sqrt(1E-8)
    abs(ores-nres);
    ores=nres;
    n=n+1;

    % New method
    for t=1:size(data,1)
        if n==1
            B(t)=median(data(t,:));
        else
            B(t)=sum(Wf(t,:).*data(t,:))/sum(Wf(t,:));
        end
    end

%     sm=gausswin(FS);
%     sm=sm/sum(sm);
%     mB=mean(B);
%     B=conv(sm,B-mean(B));
%     B=B(floor(FS/2):end-ceil(FS/2));
%     B=B+mB;
   
    if sum(isnan(B))>0
        break
    end
    if n>100
        break
    end
    res  = data-repmat(B',[1,size(data,2)]);
    
    mad  = median(abs(res-repmat(median(res')',1,size(res,2)))')'; 
    res  = (res)./repmat(mad,1,size(res,2));
    res  = res.*h; 
    res  = abs(res)-ks;
    res(res<0) = 0;   
    nres = (sum(res(:).^2));
    Wf   = (((abs(res)<1) .* (1 - res.^2).^2));
    clear res;

end

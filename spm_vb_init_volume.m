 function [slice] = spm_vb_init_volume(X,p)
% Initialise generic aspects of slice structure for VB-GLM-AR models
% FORMAT [slice] = spm_vb_init_volume(X,p)
%
% X      - design matrix
% p      - AR model order
%
% slice  - data structure (see spm_vb_glmar)
%__________________________________________________________________________
% Copyright (C) 2005-2014 Wellcome Trust Centre for Neuroimaging

% Will Penny and Nelson Trujillo-Barreto
% $Id: spm_vb_init_volume.m 6079 2014-06-30 18:25:37Z spm $

% disp('Initialising volume');
% disp(' ');


slice.X   = X;
slice.p   = p;
slice.k   = size(X,2);

T         = size(X,1);
slice.T   = T;

slice.XT  = X';
slice.XTX = slice.XT*X;
    
for t=p+1:T
    slice.dX(:,:,t-p) = X(t-1:-1:t-p,:);
end

[ux,dx,vx] = svd(X);
if size(dx,2) > 1
    ddx=diag(dx);
else
    % design matrix with a single column
    ddx = dx;
end
svd_tol  = max(ddx)*eps*slice.k;
rank_X   = sum(ddx > svd_tol);
ddxm     = diag(ones(rank_X,1)./ddx(1:rank_X));
ddxm2    = diag(ones(rank_X,1)./(ddx(1:rank_X).^2));
slice.Xp = vx(:,1:rank_X)*ddxm*ux(:,1:rank_X)';
slice.X2 = vx(:,1:rank_X)*ddxm2*vx(:,1:rank_X)';

k = slice.k;
% Get lagged input covariances for updates
slice.I.Gx   = slice.X(p+1:T,:)'*slice.X(p+1:T,:);
slice.I.xtx  = slice.X(p+1:T,:)'*slice.X(p+1:T,:);
for ki=1:k
    for kj=1:k
        dXki = squeeze(slice.dX(:,ki,:));
        dXkj = squeeze(slice.dX(:,kj,:));
        if slice.p==1
            % singleton dimension already transposed
            dXki = dXki';
            dXkj = dXkj';
        end
        Stmp = dXki*dXkj';
        slice.I.S((kj-1)*k+ki,:) = Stmp(:)';
        
        R1tmp = dXkj*slice.X(p+1:T,ki);
        slice.I.R1((kj-1)*k+ki,:) = R1tmp';
    end
end

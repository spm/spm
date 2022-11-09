function spm_ness_flows(P,x,M)
% Nonequilibrium steady-state under a Helmholtz decomposition
% FORMAT spm_ness_flows(M,x)
%--------------------------------------------------------------------------
% M   - model specification structure
% Required fields:
%    M.X  - sample points
%    M.W  - (n x n) - precision matrix of random fluctuations
%    M.K  - order of polynomial expansion
%__________________________________________________________________________

% Karl Friston
% Copyright (C) 2008-2022 Wellcome Centre for Human Neuroimaging

                
% evaluate flows at sample points 
%--------------------------------------------------------------------------
X     = spm_ndgrid(x);                   % get sample points
X     = X';
n     = numel(x);
f1    = X;
f2    = X;
f3    = X;
for i = 1:size(X,2)
    M.X     = X(:,i)';
    [F,S,Q,L,H,D] = spm_NESS_gen(P,M);
    Q       = reshape(cat(1,Q{:}),n,n);
    D       = cat(1,D{:});
    G       = diag(diag(Q));
    Q       = Q - G;
    
    f1(:,i) = G*D;
    f2(:,i) = Q*D;
    f3(:,i) = -L';
end

quiver3(X(1,:),X(2,:),X(3,:),f1(1,:),f1(2,:),f1(3,:)), hold on
quiver3(X(1,:),X(2,:),X(3,:),f2(1,:),f2(2,:),f2(3,:)), hold on
quiver3(X(1,:),X(2,:),X(3,:),f3(1,:),f3(2,:),f3(3,:)), hold on
xlabel('1st state'), ylabel('2nd state'), zlabel('3rd state')
title('Flow decomposition','Fontsize',16), axis square
legend({'Gradient flow','Solenoidal flow','Correction'})

function spm_fp_display_density(M,x)
% Quiver plot of flow and equilibrium density
% FORMAT spm_fp_display_density(M,x)
%
% M   - model specifying flow; M(1).f;
% x   - cell array of domain or support
%__________________________________________________________________________
% Copyright (C) 2005 Wellcome Trust Centre for Neuroimaging
 
% Karl Friston
% $Id: spm_fp_display_density.m 2030 2008-09-02 18:28:40Z karl $
 
% evaluation points and equilibria
%--------------------------------------------------------------------------
nx      = length(x{1});
X       = spm_ndgrid(x);
[M0,q0] = spm_fp(M,x);
 
% flow fields
%--------------------------------------------------------------------------
for i = 1:size(X,1)
    f(i,:) = feval(M(1).f,X(i,:)',0,M(1).pE)';
end
for i = 1:size(X,2)
    f(:,i) = f(:,i)/max(abs(f(:,i)));
end
 
% flow and density
%--------------------------------------------------------------------------
k     = kron(ones(1,nx/4),[1 0 0 0]);
k     = find(kron(k,k));
imagesc(x{1},x{2},1 - q0'), hold on
quiver(X(k,1),X(k,2),f(k,1),f(k,2),'r'),  hold off
axis square xy
drawnow

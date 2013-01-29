function [F,X] = spm_fp_display_density(M,x)
% Quiver plot of flow and equilibrium density
% FORMAT [F,X] = spm_fp_display_density(M,x)
%
% M   - model specifying flow; M(1).f;
% x   - cell array of domain or support
%
% F   - flow
% X   - evaluation points
%__________________________________________________________________________
% Copyright (C) 2005-2013 Wellcome Trust Centre for Neuroimaging
 
% Karl Friston
% $Id: spm_fp_display_density.m 5219 2013-01-29 17:07:07Z spm $
 
% evaluation points and equilibria
%--------------------------------------------------------------------------
n             = length(x);
[M0,q0,X,x,F] = spm_fp(M,x);
 
% flow fields
%--------------------------------------------------------------------------
for i = 1:n
    f(i,:) = F(i,:)/max(eps + abs(F(i,:)));
end


% flow and density
%==========================================================================
f  = f';

% eliminate first state if 3-D
%--------------------------------------------------------------------------
if n == 3
    q     = q0;
    q     = squeeze(sum(q,1));
    q     = squeeze(sum(q,1));
    [m,j] = max(q);
    q0    = squeeze(sum(q0,3));
    k     = find(X(:,3) == x{3}(j));
    X     = X(k,[1 2]);
    f     = f(k,[1 2]);
    x     =   x([1 2]);
end

% thin out arrows for quiver
%--------------------------------------------------------------------------
k     = 1;
for i = 1:2
    nx = length(x{i});
    d  = fix(nx/16);
    k  = kron(k,kron(ones(1,nx/d),sparse(1,1,1,1,d)));
end
k     =  find(k);


% flow and density
%--------------------------------------------------------------------------
imagesc(x{1},x{2},1 - q0'), hold on

quiver(X(k,1),X(k,2),f(k,1),f(k,2),'r'),  hold off
axis square xy
drawnow

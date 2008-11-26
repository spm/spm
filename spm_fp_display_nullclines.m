function f = spm_fp_display_nullclines(M,x)
% Nullcline plot of flow and sample trajectory
% FORMAT spm_fp_display_nullclines(M,x)
%
% M   - model specifying flow; M(1).f;
% x   - cell array of domain or support
%
% f   - derivative of x(2)
%__________________________________________________________________________
% Copyright (C) 2005 Wellcome Trust Centre for Neuroimaging

% Karl Friston
% $Id: spm_fp_display_nullclines.m 2494 2008-11-26 20:08:15Z karl $

% evaluation points and equilibria
%--------------------------------------------------------------------------
for i = 1:length(x)
    nx(i) = length(x{i});
end
X     = spm_ndgrid(x);

% flow fields
%--------------------------------------------------------------------------
for i = 1:size(X,1)
    try
        f(i,:) = feval(M(1).f,X(i,:)',0,M(1).pE)';
    catch
        f(i,:) = feval(M(1).f,X(i,:)',0,[],M(1).pE)';
    end
end

% exemplar trajectory; from M(1).x
%--------------------------------------------------------------------------
U.u   = sparse(512,M(1).m);
t     = spm_int_J(M(1).pE,M,U);
t     = [M(1).x(:)'; t];

% trajectory and nullclines
%--------------------------------------------------------------------------
f     = reshape(f(:,2),nx(1),nx(2));
image(x{1},x{2},64 - (f > 0)'*8),   hold on
plot(x{1},0*x{1},'k',t(:,1),t(:,2),'r',t(1,1),t(1,2),'or'), hold off
axis square xy
drawnow

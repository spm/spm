function spm_elipsoid(X)
% renders a square matrix onto a sphere
% FORMAT spm_elipsoid(X)
% X  - square matrix of voxel values
%
%__________________________________________________________________________
% %W% %E%

%---------------------------------------------------------------------------
n   = length(X);
r   = n/3;
y   = ones(n,1)*([1:n] - mean([1:n]));
x   = y';

d   = sqrt(x.^2 + y.^2) + eps;
a   = atan(d/r);
k   = r*sin(a)./d;
x   = x.*k;
y   = y.*k;
z   = r*cos(a) - r;

surf(x,y,-z,X.*k)
view(2)
set(gca,'AspectRatio',[1 128/174])
a   = axis;
text(a(1)*0.05,a(3)*1.1,'V','FontSize',24)
axis off


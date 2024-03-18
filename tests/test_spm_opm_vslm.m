function tests = test_spm_opm_vslm
% Unit Tests for spm_opm_vslm
%__________________________________________________________________________

% Copyright (C) 2023 Wellcome Centre for Human Neuroimaging


tests = functiontests(localfunctions);


function test_spm_opm_vslm_1(testCase)
% sphere to evaluate harmonics on
sp = spm_mesh_sphere(5);
sp.vertices= sp.vertices*50;
nv = spm_mesh_normals(sp);

% general solution
S=[];
S.li=3;
S.reg=1;
S.v=sp.vertices;
S.o=nv;
S.scale=1;
H2= spm_opm_vslm(S);

% radial orienation and position ar same on unit sphere
ex = nv(:,1);
ey = nv(:,2);
ez = nv(:,3);

% position on scaled sphere(double required for precision on non spherical)
x = double(sp.vertices(:,1));
y = double(sp.vertices(:,2));
z = double(sp.vertices(:,3));


% brute force solution 
H=[];
L1 = [ey, ez,ex ];
H= [H,L1];
L2 =[y.*ex+x.*ey, z.*ey+y.*ez, ...
    -x.*ex-y.*ey+2*z.*ez, ...
    z.*ex+x.*ez, ...
    x.*ex-y.*ey];
H= [H,L2];
L3 = [(2.*x.*y).*ex+(x.^2-y.^2).*ey,...
    (y.*z).*ex+(x.*z).*ey + (x.*y).*ez,...
    (2*x.*y).*ex+(x.^2+3*y.^2-4*z.^2).*ey+(-8*z.*y).*ez,...
    (-6*x.*z).*ex+(-6*y.*z).*ey+(6*z.^2-3*x.^2-3*y.^2).*ez,...
    (3*x.^2+y.^2-4*z.^2).*ex+(2*x.*y).*ey+(-8*x.*z).*ez,...
    (2*x.*z).*ex+(-2*y.*z).*ey+(x.^2-y.^2).*ez,...
    (3*y.^2-3*x.^2).*ex+(6*x.*y).*ey];
H= [H,L3];

% correlate general and brute force solution 
 cs = zeros(size(H2,2),1);
 for i=1:15
   Cmat = corrcoef(H2(:,i),H(:,i));
     cs(i)=abs(Cmat(1,2));
 end
 
act = sum(cs);

testCase.verifyTrue(((15-act)/15) < 1e-5);

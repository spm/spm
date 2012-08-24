function [Phi,dim1,dim2] = spm_def2sparse(PY,PI)
% Generate a sparse matrix encoding a deformation
% [Phi,dim1,dim2] = spm_def2sparse(PY,PI)
% PY - Filename of deformation field
% PI - Filename of image defining field of view etc
%_______________________________________________________________________
% Copyright (C) Wellcome Trust Centre for Neuroimaging (2009)

% John Ashburner
% $Id: spm_def2sparse.m 4861 2012-08-24 15:56:39Z john $

NY=nifti(PY);
NI=nifti(PI);

dim1 = [size(NY.dat) 1 1];
dim1 = dim1(1:3);

dim2 = [size(NI.dat) 1 1];
dim2 = dim2(1:3);

X1=NY.dat(:,:,:,1,1);
X2=NY.dat(:,:,:,1,2);
X3=NY.dat(:,:,:,1,3);

M = inv(NI.mat);

Y1=M(1,1)*X1+M(1,2)*X2+M(1,3)*X3+M(1,4);
Y2=M(2,1)*X1+M(2,2)*X2+M(2,3)*X3+M(2,4);
Y3=M(3,1)*X1+M(3,2)*X2+M(3,3)*X3+M(3,4);

if false
    % Nearest Neighbour
    fY1 = (round(Y1));
    fY2 = (round(Y2));
    fY3 = (round(Y3));
    I   = find(fY1>=1 & fY1<=dim2(1) & fY2>=1 & fY2<=dim2(2) & fY3>=1 & fY3<=dim2(3));
    J   = fY1(I) + dim2(1)*(fY2(I)-1 + dim2(2)*(fY3(I)-1));
    Phi = sparse(I,J,1,prod(dim1),prod(dim2));
else
    % Trilinear
    fY1 = (floor(Y1));
    fY2 = (floor(Y2));
    fY3 = (floor(Y3));

    I   = find(fY1>=1 & fY1<=dim2(1) & fY2>=1 & fY2<=dim2(2) & fY3>=1 & fY3<=dim2(3));
    J   = fY1(I) + dim2(1)*(fY2(I)-1 + dim2(2)*(fY3(I)-1));
    S   = (1-Y1+fY1).*(1-Y2+fY2).*(1-Y3+fY3);
    S   = S(I);
    Phi = sparse(I,J,S,prod(dim1),prod(dim2));

    I   = find(fY1>=0 & fY1<=dim2(1)-1 & fY2>=1 & fY2<=dim2(2) & fY3>=1 & fY3<=dim2(3));
    J   = fY1(I)+1 + dim2(1)*(fY2(I)-1 + dim2(2)*(fY3(I)-1));
    S   = (Y1-fY1).*(1-Y2+fY2).*(1-Y3+fY3);
    S   = S(I);
    Phi = Phi + sparse(I,J,S,prod(dim1),prod(dim2));

    I   = find(fY1>=1 & fY1<=dim2(1) & fY2>=0 & fY2<=dim2(2)-1 & fY3>=1 & fY3<=dim2(3));
    J   = fY1(I) + dim2(1)*(fY2(I) + dim2(2)*(fY3(I)-1));
    S   = (1-Y1+fY1).*(Y2-fY2).*(1-Y3+fY3);
    S   = S(I);
    Phi = Phi + sparse(I,J,S,prod(dim1),prod(dim2));

    I   = find(fY1>=0 & fY1<=dim2(1)-1 & fY2>=0 & fY2<=dim2(2)-1 & fY3>=1 & fY3<=dim2(3));
    J   = fY1(I)+1 + dim2(1)*(fY2(I) + dim2(2)*(fY3(I)-1));
    S   = (Y1-fY1).*(Y2-fY2).*(1-Y3+fY3);
    S   = S(I);
    Phi = Phi + sparse(I,J,S,prod(dim1),prod(dim2));

    I   = find(fY1>=1 & fY1<=dim2(1) & fY2>=1 & fY2<=dim2(2) & fY3>=0 & fY3<=dim2(3)-1);
    J   = fY1(I) + dim2(1)*(fY2(I)-1 + dim2(2)*fY3(I));
    S   = (1-Y1+fY1).*(1-Y2+fY2).*(Y3-fY3);
    S   = S(I);
    Phi = Phi + sparse(I,J,S,prod(dim1),prod(dim2));

    I   = find(fY1>=0 & fY1<=dim2(1)-1 & fY2>=1 & fY2<=dim2(2) & fY3>=0 & fY3<=dim2(3)-1);
    J   = fY1(I)+1 + dim2(1)*(fY2(I)-1 + dim2(2)*fY3(I));
    S   = (Y1-fY1).*(1-Y2+fY2).*(Y3-fY3);
    S   = S(I);
    Phi = Phi + sparse(I,J,S,prod(dim1),prod(dim2));

    I   = find(fY1>=1 & fY1<=dim2(1) & fY2>=0 & fY2<=dim2(2)-1 & fY3>=0 & fY3<=dim2(3)-1);
    J   = fY1(I) + dim2(1)*(fY2(I) + dim2(2)*fY3(I));
    S   = (1-Y1+fY1).*(Y2-fY2).*(Y3-fY3);
    S   = S(I);
    Phi = Phi + sparse(I,J,S,prod(dim1),prod(dim2));

    I   = find(fY1>=0 & fY1<=dim2(1)-1 & fY2>=0 & fY2<=dim2(2)-1 & fY3>=0 & fY3<=dim2(3)-1);
    J   = fY1(I)+1 + dim2(1)*(fY2(I) + dim2(2)*fY3(I));
    S   = (Y1-fY1).*(Y2-fY2).*(Y3-fY3);
    S   = S(I);
    Phi = Phi + sparse(I,J,S,prod(dim1),prod(dim2));
end


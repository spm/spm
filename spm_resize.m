function [O] = spm_resize(I,x,y)
% resizes a matrix
% FORMAT [O] = spm_resize(I,X,Y);
% I - image
% X - rows
% Y - columns
% O - resized {X x Y} image
%______________________________________________________________________
%
% spm_resize resizes a matrix (I) in working memory to x x y using
% bilinear interpolation
%
%__________________________________________________________________________
% %W% %E%

%----------------------------------------------------------------------
[m n] = size(I);
O     = interp2(([1:n]-1),([1:m]-1)',I,...
		([1:y]-1)/(y-1)*(n-1),([1:x]-1)'/(x-1)*(m-1));

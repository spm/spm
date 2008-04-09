function [y] = plgndr(n,k,x)

% PLGNDR associated Legendre function
% 
% y = plgndr(n,k,x) computes the values of the associated Legendre functions
% of degree N and order K
%
% implemented as MEX file

% the original implementation was based on "Numerical Recipes in C", version 2.0
% but has been replaced with an equvalent function from GNU Scientific Library

% Copyright (C) 2002, Robert Oostenveld
%
% $Log: plgndr.m,v $
% Revision 1.3  2008/03/05 16:26:18  roboos
% updated documentation
%
% Revision 1.2  2003/03/11 14:45:37  roberto
% updated help and copyrights
%

error(sprintf('could not locate MEX file for %s', mfilename))


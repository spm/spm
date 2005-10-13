function [y] = plgndr(n,k,x)

% PLGNDR associated Legendre function
% 
% y = plgndr(n,k,x) computes the values of the associated Legendre functions
% of degree N and order K
%
% implemented as MEX file

% the actual implementation is based on "Numerical Recipes in C", version 2.0

% Copyright (C) 2002, Robert Oostenveld
%
% $Log: plgndr.m,v $
% Revision 1.2  2003/03/11 14:45:37  roberto
% updated help and copyrights
%

error(sprintf('could not locate MEX file for %s', mfilename))


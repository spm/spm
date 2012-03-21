
% Script to compile C code for DARTEL registration

% $Id: make.m 4698 2012-03-21 14:00:44Z john $
% John Ashburner
mex dartel3.c optimizer3d.c diffeo3d.c -O
mex optimN.c optimN_mex.c -O
mex optimN.c optimN_mex.c -O -DNEUMANN -output optimNn


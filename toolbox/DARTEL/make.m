
% Script to compile C code for DARTEL registration

mex dartel2.c optimizer2d.c diffeo2d.c -O
mex dartel3.c optimizer3d.c diffeo3d.c -O
mex optimizer2d.c optimizer2d_mex.c -O
mex optimizer3d.c optimizer3d_mex.c -O
mex optimizer2d.c optimizer2d_mex.c -O -DNEUMANN -output optimizer2dn
mex optimizer3d.c optimizer3d_mex.c -O -DNEUMANN -output optimizer3dn
mex optim1.c optim1_mex.c -O
mex optimN.c optimN_mex.c -O
mex optimN.c optimN_mex.c -O -DNEUMANN -output optimNn


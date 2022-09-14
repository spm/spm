% Script to compile C code for EEG-BEM

% Copyright (C) 1999-2022  Christophe Phillips

mex -O bem_Cii_cog.c
mex -O bem_Cii_cst.c
mex -O bem_Cii_lin.c
mex -O bem_Cij_cog.c
mex -O bem_Cij_cst.c
mex -O bem_Cij_lin.c
mex -O bem_Gi_cog.c
mex -O bem_Gi_vert.c

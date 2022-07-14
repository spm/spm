%
% Script to compile C MEX-files for FieldMap
%
% See also: mex, Makefile
%__________________________________________________________________________

% Chloe Hutton
% Copyright (C) 2008-2022 Wellcome Centre for Human Neuroimaging


mex -O pm_invert_phasemap_dtj.c
mex -O pm_merge_regions.c
mex -O pm_create_connectogram_dtj.c
mex -O pm_pad.c
mex -O pm_estimate_ramp.c
mex -O pm_restore_ramp.c
mex -O pm_ff_unwrap.c

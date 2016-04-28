%
% Script to compile C MEX-files for FieldMap
%
% See also: mex, Makefile
%__________________________________________________________________________
% Copyright (C) 2008-2011 Wellcome Trust Centre for Neuroimaging

% Chloe Hutton
% $Id: make_FieldMap.m 6791 2016-04-28 14:47:20Z john $

mex -O pm_invert_phasemap_dtj.c            
mex -O pm_merge_regions.c
mex -O pm_create_connectogram_dtj.c  
mex -O pm_pad.c
mex -O pm_estimate_ramp.c
mex -O pm_restore_ramp.c
mex -O pm_ff_unwrap.c                

% make.m
%
% This is a matlab makefile which simply runs through all the 
% necessary .c files compiling them individually using mex. 
% It may be necessary to use any special options.
% if there are problems, in matlab check:
% >help mex.
%--------------------------------------------------------------------------
% Chloe Hutton

mex ip_bwlabel.c 
mex pm_invert_phasemap_dtj.c
mex ip_dilate_erode.c             
mex pm_merge_regions.c
mex pm_create_connectogram_dtj.c  
mex pm_pad.c
mex pm_estimate_ramp.c
mex pm_restore_ramp.c
mex pm_ff_unwrap.c                
mex pm_smooth_phasemap_dtj.c

% If the above commands do not work,it may be necessary to specify another 
% location for the options file. 
% In this case, comment out the above lines and uncomment the lines below
% replacing /usr/local/matlab6.5 with the location of your local version of 
% matlab.

%mex -f /local/matlab6.5/bin/gccopts.sh ip_bwlabel.c 
%mex -f /local/matlab6.5/bin/gccopts.sh pm_invert_phasemap_dtj.c
%mex -f /local/matlab6.5/bin/gccopts.sh ip_dilate_erode.c             
%mex -f /local/matlab6.5/bin/gccopts.sh pm_merge_regions.c
%mex -f /local/matlab6.5/bin/gccopts.sh pm_create_connectogram_dtj.c  
%mex -f /local/matlab6.5/bin/gccopts.sh pm_pad.c
%mex -f /local/matlab6.5/bin/gccopts.sh pm_estimate_ramp.c
%mex -f /local/matlab6.5/bin/gccopts.sh pm_restore_ramp.c
%mex -f /local/matlab6.5/bin/gccopts.sh pm_ff_unwrap.c                
%mex -f /local/matlab6.5/bin/gccopts.sh pm_smooth_phasemap_dtj.c

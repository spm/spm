
% Rendering of regional effects [SPM{Z}] on drawings of the cortical surface
% FORMAT spm_picture
%____________________________________________________________________________
%
% spm_picture has been replaced by the SPM rendering utility - spm_render.m
%
%                           ----------------
%
% spm_picture is called after spm_results_ui, and uses variables in working
% memory to load images/pictures of the cortical surface.  Regional
% foci from the selected SPM{Z} are rendered on these drawings.
%
% These are crude renders usually used for didactic purposes only.  The
% voxels are assigned to four parasaggistal blocks (x < -24, 0 > x > -24
% 0 < x  x< 24 and x > 24) and are displayed (as a MIP) on the appropriate
% surface.
%__________________________________________________________________________
% %W% Karl Friston %E%

%-Print warning of obsolescence
%-----------------------------------------------------------------------
warning('spm_picture is obsolete: use spm_render instead')


%-NB: Right & left hemisphere cortex & lateral surface drawings are
% stored in the MAT files spm_rctx.mat & spm_lctx.mat

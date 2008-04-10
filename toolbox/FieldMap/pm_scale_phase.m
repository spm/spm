% Script to scale a phase map so that max = pi and min =-pi radians.
% Writes out scaled image prepended with 'sc'.
% Chloe Hutton 25/02/04
% SPM Update - 13/11/06
%_______________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% Chloe Hutton
% $Id: pm_scale_phase.m 1358 2008-04-10 11:20:26Z guillaume $

spm('defaults','FMRI');

V=spm_vol(spm_select(1,'image','Select phase image to scale'));
vol=spm_read_vols(V);

mn=min(vol(:));
mx=max(vol(:));
svol=-pi+(vol-mn)*2*pi/(mx-mn);

oV=V;
name='sc';
oV.fname=[spm_str_manip(V.fname, 'h') ['/' name] deblank(spm_str_manip(V.fname,'t'))];
spm_write_vol(oV,svol);

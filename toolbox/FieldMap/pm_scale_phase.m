% Script to scale a phase map so that max = pi and min =-pi radians.
% Writes out scaled image prepended with 'sc'.
%__________________________________________________________________________

% Chloe Hutton
% Copyright (C) 2004-2022 Wellcome Centre for Human Neuroimaging


V   = spm_vol(spm_select(1,'image','Select phase image to scale'));
vol = spm_read_vols(V);

mn   = min(vol(:));
mx   = max(vol(:));
svol = -pi+(vol-mn)*2*pi/(mx-mn);

% Output image struct
oV = struct(...
    'fname',   spm_file(V.fname,'prefix','sc'),...
    'dim',     V.dim(1:3),...
    'dt',      [4 spm_platform('bigend')],...
    'mat',     V.mat,...
    'descrip', 'Scaled phase');

spm_write_vol(oV,svol);

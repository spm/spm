
% Little script to generate the BEM image where the 3 volumes (i/o-skull &
% scalp) are coded by intensity in the image, as well as their outside 
% "border".
%_______________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% Christophe Phillip
% $Id: BuildBEMimg.m 1930 2008-07-18 11:49:51Z christophe $


P = spm_select(3,'image','Select i/o-skull & scalp');
V = spm_vol(P);

e_iskull = spm_erode(V(1).private.dat(:,:,:));
e_oskull = spm_erode(V(2).private.dat(:,:,:));
e_scalp = spm_erode(V(3).private.dat(:,:,:));

% Keeping Robert's intensity convention:
% - brain volume,           14 = 12+2
% - brain volume border,    94 = 14+80
% - skull volume,           12 = 8+4
% - skull volume border,    92 = 12+80
% - scalp volume,           8
% - scalp volume border,    24 = 8+16
% - background,             0

bem = V(3).private.dat(:,:,:)*8 ...                 % scalp volume
      + (V(3).private.dat(:,:,:)-e_scalp)*16 ...    % scalp vol. border
      + V(2).private.dat(:,:,:)*4 ...               % skull volume
      + (V(2).private.dat(:,:,:)-e_oskull)*80 ...   % skull vol. border
      + V(1).private.dat(:,:,:)*2 ...               % brain volume
      + (V(1).private.dat(:,:,:)-e_iskull)*80 ;     % brain vol. border
  
Vbem = V(1);
Vbem.fname = [spm_str_manip(Vbem.fname,'h'),filesep,'kk_bem.nii'];
Vbem.pinfo(1) = 1;
Vbem = spm_write_vol(Vbem,bem);

function VO = spm_add2(VI,VO,flg)
% add a series of images - a compiled routine
% FORMAT VO = spm_add2(VI,VO)
% VI    - Vector of volume handles (from spm_vol).
% VO    - Description of output volume that gets passed to
%         spm_write_vol.m
% flg - Flags can be:
%       'm' - masks the mean to zero or NaN wherever
%             a zero occurs in the input images.
%_______________________________________________________________________
%
% spm_add computes a sum of a set of image volumes to produce an
% integral image that is written to a named file (VO.fname).
%
% A mean can be effected by modifying the scalefactors (and offsets) of
% VI. A weighted sum can be effected by using different weightings for
% image scalefactors.
%
% This function is intended to replace spm_add.
%_______________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% John Ashburner
% $Id: spm_add2.m 1143 2008-02-07 19:33:33Z spm $

if nargin>2 && any(flg=='m'), msk = true;
else msk = false; end;

dat = zeros(VI(1).dim(1:3));
for j=1:numel(VI),
    for i=1:VI(1).dim(3),
        slice = spm_slice_vol(VI(j),spm_matrix([0 0 i]),VI(1).dim(1:2),0);
        if msk, slice(slice==0) = NaN; end;
        dat(:,:,i) = dat(:,:,i) + slice;
    end;
end;
VO = spm_write_vol(VO,dat);


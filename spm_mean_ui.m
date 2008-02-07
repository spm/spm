function spm_mean_ui
% promts for a series of images and averages them
% FORMAT spm_mean_ui
%_______________________________________________________________________
%
% spm_mean_ui simply averages a set of images to produce a mean image
% that is written as type int16 to "mean.img" (in the current directory).
%
% The images must have the same dimensions, orientations (as defined by
% the Origin header field or any associated *.mat files), and the same
% voxel sizes.
%
% This is not a "softmean" - zero voxels are treated as zero.
%_______________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% John Ashburner, Andrew Holmes
% $Id: spm_mean_ui.m 1143 2008-02-07 19:33:33Z spm $

SCCSid = '$Rev: 1143 $';


%-Say hello
%-----------------------------------------------------------------------
SPMid = spm('FnBanner',mfilename,SCCSid);


%-Select images & check dimensions, orientations and voxel sizes
%-----------------------------------------------------------------------
fprintf('\t...select files')
P = spm_select(Inf,'image','Select images to be averaged');
fprintf(' ...mapping & checking files')
Vi = spm_vol(P);

n  = prod(size(Vi));
if n==0, fprintf('\t%s : no images selected\n\n',mfilename), return, end

spm_check_orientations(Vi);


%-Compute mean and write headers etc.
%-----------------------------------------------------------------------
fprintf(' ...computing')
Vo = struct(    'fname',    'mean.img',...
        'dim',      Vi(1).dim(1:3),...
        'dt',           [4, spm_platform('bigend')],...
        'mat',      Vi(1).mat,...
        'pinfo',    [1.0,0,0]',...
        'descrip',  'spm - mean image');

%-Adjust scalefactors by 1/n to effect mean by summing
for i=1:prod(size(Vi))
    Vi(i).pinfo(1:2,:) = Vi(i).pinfo(1:2,:)/n; end;

Vo            = spm_create_vol(Vo);
Vo.pinfo(1,1) = spm_add(Vi,Vo);
Vo            = spm_create_vol(Vo);


%-End - report back
%-----------------------------------------------------------------------
fprintf(' ...done\n')
fprintf('\tMean image written to file ''%s'' in current directory\n\n',Vo.fname)

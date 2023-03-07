function spm_mean(P)
% Compute a mean image from a set
% FORMAT spm_mean(P)
% P   - list of images to average [Default: GUI]
%__________________________________________________________________________
%
% spm_mean_ui simply averages a set of images to produce a mean image
% that is written as type int16 to "mean.img" (in the current directory).
%
% The images must have the same dimensions, orientations and the same
% voxel sizes.
%
% This is not a "softmean" - zero voxels are treated as zero.
%__________________________________________________________________________

% John Ashburner, Andrew Holmes
% Copyright (C) 1998-2022 Wellcome Centre for Human Neuroimaging


persistent runonce
if isempty(runonce)
    warning('ImCalc should be preferred to spm_mean.');
    runonce = 1;
end

%-Say hello
%--------------------------------------------------------------------------
spm('FnBanner',mfilename);

%-Select images & check dimensions, orientations and voxel sizes
%--------------------------------------------------------------------------
if ~nargin
    [P, sts] = spm_select([2 Inf],'image','Select images to be averaged');
    if ~sts, return; end
end

Vi = spm_vol(P);

spm_check_orientations(Vi);

%-Compute mean and write image
%--------------------------------------------------------------------------
Vo = struct('fname',    ['mean' spm_file_ext],...
            'dim',      Vi(1).dim(1:3),...
            'dt',       [spm_type('int16') spm_platform('bigend')],...
            'mat',      Vi(1).mat,...
            'pinfo',    [1 0 0]',...
            'descrip',  'spm - mean image');

%-Adjust scalefactors by 1/n to effect mean by summing
n  = numel(Vi);
for i=1:n
    Vi(i).pinfo(1:2,:) = Vi(i).pinfo(1:2,:)/n;
end

Vo            = spm_create_vol(Vo);
Vo.pinfo(1,1) = spm_add(Vi,Vo);
Vo            = spm_create_vol(Vo);

%-End - report back
%--------------------------------------------------------------------------
fprintf('Mean image written to file ''%s'' in current directory\n',Vo.fname);

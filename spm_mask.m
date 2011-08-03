function spm_mask(P1, P2, thresh)
% Mask images
% FORMAT spm_mask(P1, P2, thresh)
% P1     - matrix of input image filenames from which
%          to compute the mask.
% P2     - matrix of input image filenames on which
%          to apply the mask.
% thresh - optional threshold(s) for defining the mask.
% The masked images are prepended with the prefix `m'.
%
% If any voxel in the series of images is zero (for data types without
% a floating point representation) or does not have a finite value (for
% floating point and double precision images), then that voxel is set to
% NaN or zero in all the images.  If a threshold, or vector of
% thresholds is passed, then the masking is based on voxels whos
% values are above all the thresholds.
%
% Images sampled in different orientations and positions can be passed
% to the routine.
%__________________________________________________________________________
% Copyright (C) 1999-2011 Wellcome Trust Centre for Neuroimaging

% John Ashburner
% $Id: spm_mask.m 4419 2011-08-03 18:42:35Z guillaume $


persistent runonce
if isempty(runonce)
    warning('ImCalc should be preferred to spm_mask whenever possible.');
    runonce = 1;
end

SVNid = '$Rev: 4419 $';

%-Say hello
%--------------------------------------------------------------------------
SPMid = spm('FnBanner',mfilename,SVNid);

%-Parameters & Arguments
%--------------------------------------------------------------------------
if ~nargin
    [P1, sts] = spm_select(Inf,'image','Images to compute mask from');
    if ~sts, return; end
    [P2, sts] = spm_select(Inf,'image','Images to apply mask to');
    if ~sts, return; end
end

if nargin==1
    P2 = P1;
end

V1 = spm_vol(P1);
V2 = spm_vol(P2);

if nargin==3
    if numel(thresh)==1
        thresh = repmat(thresh,numel(V1),1);
    else
        if numel(V1) ~= numel(thresh)
            error('Input argument ''thresh'' has wrong size.');
        end
        thresh = thresh(:);
    end
end

m1 = numel(V1);
m2 = numel(V2);

%-Create headers
%--------------------------------------------------------------------------
VO = V2;
for i=1:m2
    [pth,nm,ext,num] = spm_fileparts(deblank(VO(i).fname));
    VO(i).fname      = fullfile(pth,['m', nm, ext, num]);
    VO(i).descrip    = 'Masked';
    VO(i).mat        = VO(1).mat;
    VO(i).dim(1:3)   = VO(1).dim(1:3);
end
VO  = spm_create_vol(VO);
M   = VO(1).mat;
dim = VO(1).dim(1:3);

%-Compute masked images
%--------------------------------------------------------------------------
spm_progress_bar('Init',VO(1).dim(3),'Masking','planes completed');

for j=1:dim(3)

    msk = true(dim(1:2));
    Mi  = spm_matrix([0 0 j]);

    % Load slice j from all images
    for i=1:m1
        M1  = M\V1(i).mat\Mi;
        %if sum((M1(:)-Mi(:)).^2<eps) M1 = Mi; end;

        img = spm_slice_vol(V1(i),M1,dim(1:2),[0 NaN]);
        
        msk = msk & isfinite(img);
        
        if ~spm_type(V1(i).dt(1),'nanrep')
            msk = msk & (img ~= 0);
        end
        
        if nargin == 3
            msk = msk & (img >= thresh(i));
        end
    end
    
    % Write the images
    for i=1:m2
        M1        = M\V2(i).mat\Mi;
        img       = spm_slice_vol(V2(i),M1,dim(1:2),[1 0]);
        img(~msk) = NaN;
        VO(i)     = spm_write_plane(VO(i),img,j);
    end

    spm_progress_bar('Set',j);
end

spm_progress_bar('Clear');

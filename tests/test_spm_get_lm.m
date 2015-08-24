function tests = test_spm_get_lm
% Unit Tests for spm_get_lm
%__________________________________________________________________________
% Copyright (C) 2015 Wellcome Trust Centre for Neuroimaging

% $Id: test_spm_get_lm.m 6534 2015-08-24 16:02:56Z guillaume $

tests = functiontests(localfunctions);


function test_spm_get_lm_2D(testCase)
dim     = [5 5];
vol     = rand(dim(1),dim(2));
subvol  = vol(2:4,2:4);
vol(3,3) = max(subvol(:)) + 0.1;
[X,Y] = ndgrid(1:dim(1),1:dim(2));
gm      = sub2ind(size(vol),3,3);
idx     = spm_get_lm(vol,[X(:)';Y(:)']);
testCase.verifyTrue(any(idx == gm));


function test_spm_get_lm_3D(testCase)
dim     = [5 5 5];
vol     = rand(dim(1),dim(2),dim(3));
subvol  = vol(2:4,2:4,2:4);
vol(3,3,3) = max(subvol(:)) + 0.1;
[X,Y,Z] = ndgrid(1:dim(1),1:dim(2),1:dim(3));
gm      = sub2ind(size(vol),3,3,3);
idx     = spm_get_lm(vol,[X(:)';Y(:)';Z(:)']);
testCase.verifyTrue(any(idx == gm));

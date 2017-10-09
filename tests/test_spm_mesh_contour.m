function tests = test_spm_mesh_contour
% Unit Tests for spm_mesh_contour
%__________________________________________________________________________
% Copyright (C) 2017 Wellcome Trust Centre for Neuroimaging

% $Id: test_spm_mesh_contour.m 7183 2017-10-09 15:26:47Z guillaume $

tests = functiontests(localfunctions);


function test_spm_mesh_contour_(testCase)
M = gifti(fullfile(spm('Dir'),'canonical','cortex_20484.surf.gii'));
S = spm_mesh_contour(M,mean(M.vertices(:,3)));

% at least two hemispheres
exp = 2;
act = numel(S);
testCase.verifyGreaterThanOrEqual(act, exp);

% all surfaces are closed
testCase.verifyFalse(any([S.isopen]));

function tests = test_spm_fileparts
% Unit Tests for spm_fileparts
%__________________________________________________________________________

% Copyright (C) 2022 Wellcome Centre for Human Neuroimaging


tests = functiontests(localfunctions);


function test_spm_fileparts_char(testCase)

[pth,nam,ext,num] = spm_fileparts('');
testCase.verifyEqual(pth, '');
testCase.verifyEqual(nam, '');
testCase.verifyEqual(ext, '');
testCase.verifyEqual(num, '');

[pth,nam,ext,num] = spm_fileparts(spm('Dir'));
[p,n] = fileparts(spm('Dir'));
testCase.verifyEqual(pth, p);
testCase.verifyEqual(nam, n);
testCase.verifyEqual(ext, '');
testCase.verifyEqual(num, '');

[pth,nam,ext,num] = spm_fileparts([spm('Dir') filesep]);
testCase.verifyEqual(pth, spm('Dir'));
testCase.verifyEqual(nam, reshape('',1,0));
testCase.verifyEqual(ext, '');
testCase.verifyEqual(num, '');

[pth,nam,ext,num] = spm_fileparts(fullfile(spm('Dir'),'spm.m'));
testCase.verifyEqual(pth, spm('Dir'));
testCase.verifyEqual(nam, 'spm');
testCase.verifyEqual(ext, '.m');
testCase.verifyEqual(num, '');

[pth,nam,ext,num] = spm_fileparts(fullfile(spm('Dir'),'image.nii,1'));
testCase.verifyEqual(pth, spm('Dir'));
testCase.verifyEqual(nam, 'image');
testCase.verifyEqual(ext, '.nii');
testCase.verifyEqual(num, ',1');

[pth,nam,ext,num] = spm_fileparts(fullfile(spm('Dir'),'image.nii,2,5'));
testCase.verifyEqual(pth, spm('Dir'));
testCase.verifyEqual(nam, 'image');
testCase.verifyEqual(ext, '.nii');
testCase.verifyEqual(num, ',2,5');

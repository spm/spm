classdef test_spm_bids < matlab.unittest.TestCase
% Unit Tests for spm_bids
%   spm_tests('class','unit','test','test_spm_bids')
%__________________________________________________________________________

% Copyright (C) 2025 Wellcome Centre for Human Neuroimaging

properties
   bidsdir % temporary bids root folder
end

methods (TestMethodSetup)
   function setup(testCase)
      import matlab.unittest.fixtures.TemporaryFolderFixture
      % per test bids root. matlab removes when test finishes
      tmp = testCase.applyFixture(TemporaryFolderFixture);
      testCase.bidsdir = tmp.Folder;
      % Note: R2023 introduced 'testCase.createTemporaryFolder();'
      %       more ergonomic than instead of applyFixture().Folder

      % required json at root of BIDS/
      desc_f = fullfile(testCase.bidsdir,'dataset_description.json');
      fid = fopen(desc_f,'w');
      fprintf(fid, '{"Name": "Test", "BIDSVersion": "1.7.0", "Authors": ["None"]}');
      fclose(fid);
   end
end % TestMethodSetup


methods (Test)
function test_spm_bids_dwi(testCase)
import matlab.unittest.constraints.*


% create files:
%  /tmp/bids/sub-1/dwi/sub-1_dir-AP.{nii,json,bvec,bval}
dwidir = fullfile(testCase.bidsdir, 'sub-1','dwi');
mkdir(dwidir);
prefix = 'sub-1_dir-AP_dwi';
% Note: nii just needs to exist for this test. value/header is not read
for ext = {'nii','bvec','bval'};
   bids_file = fullfile(dwidir,[prefix '.' ext{1}]);
   fid = fopen(bids_file,'w');
   fprintf(fid,'1'); % any value will do?
   fclose(fid);
end
% json handled individually
fid = fopen(fullfile(dwidir,[prefix '.json']),'w');
fprintf(fid,'{}');
fclose(fid);

% run BIDS
BIDS = spm_BIDS(testCase.bidsdir);

% confirm we were able to read single dwi
testCase.verifyEqual(fieldnames(BIDS), {'dir','description','sessions','scans','sess','participants','subjects'}');
testCase.assertTrue(ismember('dwi', fieldnames(BIDS.subjects(1))));
testCase.verifyEqual(BIDS.subjects(1).dwi(1).bval, 1);
end % dwi test

end % method
end % classdef

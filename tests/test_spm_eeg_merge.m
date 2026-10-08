classdef test_spm_eeg_merge < matlab.unittest.TestCase
% Unit Tests for spm_eeg_merge
%__________________________________________________________________________

% Copyright (C) 2023 Wellcome Centre for Human Neuroimaging



methods (TestMethodSetup)
    function useTemporaryFolder(testCase)
        % Each test runs in a new temporary folder, which receives the output
        testCase.applyFixture(matlab.unittest.fixtures.WorkingFolderFixture);
    end
end % methods (TestMethodSetup)

methods (Test)


function test_spm_eeg_merge_1(testCase)

spm('defaults','eeg');

fname = fullfile(spm('Dir'),'tests','data','OPM','test_opm.mat');
tname = fullfile(spm('Dir'),'tests','data','OPM','test_opm2.dat');

D     = spm_eeg_load(fname);
D = chantype(D,1:110,'MEG');
D.save();

D2 = clone(D,tname);
D2.save()


S=[];
S.D={D,D2};
Dout = spm_eeg_merge(S);
siD = size(Dout);

Dout.delete();
D2.delete();

testCase.verifyTrue(all(siD == [110,1000,2]));
end
end % methods (Test)

end % classdef
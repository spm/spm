classdef test_spm_eeg_crop < matlab.unittest.TestCase
% Unit Tests for spm_eeg_crop
%__________________________________________________________________________

% Copyright (C) 2023 Wellcome Centre for Human Neuroimaging

methods (TestClassSetup)
    function setupSPM(testCase)
        spm('defaults', 'eeg');
    end
end % methods (TestClassSetup)

methods (Test)


function test_spm_eeg_crop_1(testCase)

fname = fullfile(spm('Dir'),'tests','data','OPM','test_opm.mat');
D     = spm_eeg_load(fname);
D = chantype(D,1:110,'MEG');
D.save();

labs = chanlabels(D);

S=[];
S.D=D;
S.timewin = [1 500];
S.channels = {labs{1:30}};
cD = spm_eeg_crop(S);
siD = size(cD);

delete(cD);

testCase.verifyTrue(all(siD==[30, 500, 1]));
end

end % methods (Test)

end % classdef
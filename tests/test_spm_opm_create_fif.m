classdef test_spm_opm_create_fif < matlab.unittest.TestCase
% Unit Tests for spm_opm_hfc
%__________________________________________________________________________

% Copyright (C) 2023 Wellcome Centre for Human Neuroimaging



methods (Test)


function test_spm_opm_create_fif_1(testCase)

    spm('defaults','eeg');
    
    data = fullfile(spm('Dir'),'tests','data','OPM','OPM_meg_fif_part_pos.fif');
    
    S = [];
    S.data = data;
    D = spm_opm_create(S);
    
    testCase.verifyTrue(size(D.sensors('MEG').coilpos,1)==75);
end

function test_spm_opm_create_fif_2(testCase)

    spm('defaults','eeg');
    
    data = fullfile(spm('Dir'),'tests','data','OPM','OPM_meg_fif_no_pos.fif');
 
    S = [];
    S.data = data;
    D = spm_opm_create(S);
    
    testCase.verifyTrue(~isfield(D, 'sensors'));
end

end % methods (Test)

end % classdef
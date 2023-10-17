function [coverage, tocover] = spm_eeg_test_coverage
% Return number of M/EEG functions and number of associated tests
% FORMAT [coverage, tocover] = spm_eeg_test_coverage
%
% Output:
%   coverage - number of M/EEG tested functions
%   tocover  - number of M/EEG functions 
%__________________________________________________________________________

% Tim Tierney
% Copyright (C) 2023 Wellcome Centre for Human Neuroimaging


fpmfiles   = spm_select('FPListRec',spm('Dir'),'^spm_(eeg|opm).*');
mfiles     = spm_file(cellstr(fpmfiles),'basename');
tocover    = numel(mfiles);

coverage   = 0;
testfolder = fullfile(spm('Dir'),'tests');
for i=1:tocover
    if spm_existfile(fullfile(testfolder,['test_' mfiles{i},'.m']))
        coverage = coverage + 1;
    end
end

fprintf('spm M/EEG test coverage is %d out of %d functions.\n',...
    coverage, tocover);

function [coverage, tocover] = spm_eeg_test_coverage()
% Returns number of M/EEG functions and number of associated tests
% FORMAT  [coverage, tocover] = spm_eeg_test_coverage()
%
% Output:
%   coverage          - number of M/EEG tests
%   tocover           - number of M/EEG funcitons 
%__________________________________________________________________________
% Copyright  Tim Tierney



[fpmfiles,~] = spm_select('FPListRec',spm('dir'),'^spm_(eeg|opm).*');
[~,mfiles]=fileparts(string(fpmfiles));
coverage = 0;
testfolder = fullfile(spm('dir'),'tests');
tocover = size(mfiles,1);

for i = 1:tocover
[testfiles,~] = spm_select('List',testfolder,[[mfiles{i}] '.*m$']);
if ~isempty(testfiles)
  coverage=coverage+1; 
end
end

msg = ['spm M/EEG test coverage is ' num2str(coverage) ....
  ' out of ' num2str(tocover) ' functions'];
disp(msg);
end




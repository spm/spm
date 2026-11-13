classdef test_checkcode < matlab.unittest.TestCase
% Test for possible problems in all of MATLAB code files
%__________________________________________________________________________

% Copyright (C) 2017-2022 Wellcome Centre for Human Neuroimaging


methods (Test)

function test_display_matlab_version(testCase)
% Output MATLAB version for reporting to Mathworks if needed
ver -support
end

function test_checkcode_all(testCase)
testCase.assumeTrue(strcmp(spm_check_version,'matlab'));
testCase.assumeFalse(isdeployed);

f = cellstr(spm_select('FPlistRec',spm('Dir'),'\.m$')); % discard external?
s = checkcode(f,'-id');
checked = false(size(s));
for i=1:numel(s)
    for j=1:numel(s{i})
        if isequal(s{i}(j).id,'SYNER') % Syntax error
            checked(i) = true;
            fprintf('%s (line %d): %s\n', f{i}, s{i}(j).line, s{i}(j).message)
        end
    end
end
testCase.verifyFalse(any(checked));
if any(checked)
    fprintf('Problematic files are:\n');
    disp(f(checked));
end
end

end % methods (Test)


methods (Test, TestTags = {'Disabled'})

function test_analyzeCodeCompatibility(testCase) % disabled for now
testCase.assumeFalse(~strcmp(spm_check_version,'matlab') | ...
    spm_check_version('matlab','9.3') < 0);
    
import matlab.unittest.constraints.*
r = analyzeCodeCompatibility(spm('Dir'),'IncludeSubfolders',true); % discard external?
idx = find(r.Recommendations.Severity == 'Error');
testCase.verifyThat(numel(idx), IsEqualTo(0));
for i=1:numel(idx)
    fprintf('%s:\n',r.Recommendations.File(idx(i)));
    fprintf('  %s\n',r.Recommendations.Description(idx(i)));
end

end

end % methods (Test)

end % classdef
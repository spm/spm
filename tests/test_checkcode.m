function tests = test_checkcode
% Test for possible problems in all of MATLAB code files
%__________________________________________________________________________

% Copyright (C) 2017-2022 Wellcome Centre for Human Neuroimaging


tests = functiontests(localfunctions);

function test_display_matlab_version(testCase)
% Output MATLAB version for reporting to Mathworks if needed
ver -support

function test_checkcode_all(testCase)
testCase.assumeTrue(strcmp(spm_check_version,'matlab'));

f = cellstr(spm_select('FPlistRec',spm('Dir'),'\.m$')); % discard external?
s = checkcode(f,'-id');
checked = false(size(s));
for i=1:numel(s)
    for j=1:numel(s{i})
        if isequal(s{i}(j).id,'SYNER') % Syntax error
            checked(i) = true;
        end
    end
end
testCase.verifyFalse(any(checked));
if any(checked)
    fprintf('Problematic files are:\n');
    disp(f(checked));
end


function t_e_s_t_analyzeCodeCompatibility(testCase) % disabled for now
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

function tests = test_checkcode
% Test for possible problems in all of MATLAB code files
%__________________________________________________________________________
% Copyright (C) 2017 Wellcome Trust Centre for Neuroimaging

% $Id: test_checkcode.m 7107 2017-06-14 09:31:40Z guillaume $

tests = functiontests(localfunctions);

function test_checkcode_all(testCase)

f = cellstr(spm_select('FPlistRec',spm('Dir'),'\.m$')); % discard external?
s = checkcode(f,'-id');
check = false(size(s));
for i=1:numel(s)
    for j=1:numel(s{i})
        if isequal(s{i}(j).id,'SYNER'); % Syntax error
            check(i) = true;
        end
    end
end
testCase.verifyFalse(any(check));
if any(check)
    fprintf('Problematic files are:\n');
    disp(f(check));
end

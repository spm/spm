function spm_eeg_hist2script(S)
% function that generates a script from the history of an SPM for M/EEG
% file. % FORMAT D = spm_eeg_hist2script(S)
%
% S  - filename or input struct (optional)
% (optional) fields of S:
% history         - history of M/EEG object (D.history)
% fname           - filename of to be generated script
%__________________________________________________________________________
%
% In SPM for M/EEG, each preprocessing step enters its call and input
% arguments into an internal history. The sequence of function calls that
% led to a given file can be read by the history method (i.e. call
% 'D.history'). From this history this function generates a script (m-file)
% which can be run without user interaction and will repeat, if run, the
% exact sequence on the preprocessing steps stored in the history. Of
% course, the generated script can also be used as a template for a
% slightly different analysis or for different subjects.
%__________________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% Stefan Kiebel
% $Id: spm_eeg_hist2script.m 2861 2009-03-11 18:41:03Z guillaume $

try
    h = S.history;
catch
    [D, sts] = spm_select(1, 'mat', 'Select M/EEG mat file');
    if ~sts, return; end
    D = spm_eeg_load(D);
    h = D.history;
end

try
    S.fname;
catch
    [filename, pathname] = uiputfile('*.m', 'Select to be generated script file');
    S.fname = fullfile(pathname, filename);
end

Nh = length(h);
fp = fopen(S.fname, 'w+');

for i = 1:Nh

    fprintf(fp, '%s\n', 'S = [];');
    s = gencode(h(i).args(1), 'S');

    for j = 1:length(s)
        fprintf(fp, '%s\n', s{j});
    end
    
    fprintf(fp, '%s\n\n\n', ['D = ' eval('h(i).fun') '(S);']);
    
end
fclose(fp);



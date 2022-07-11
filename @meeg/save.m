function this = save(this)
% Save an meeg object into a file
% FORMAT this = save(this)
%
% Converts an meeg object to struct and saves it.
%__________________________________________________________________________

% Vladimir Litvak
% Copyright (C) 2008-2022 Wellcome Centre for Human Neuroimaging


D = struct(this);
try
    save(fullfile(this), 'D', spm_get_defaults('mat.format'));
catch
    err = lasterror;
    if strcmp(err.identifier, 'MATLAB:MatFile:UnsupportedCharacters')
        warning('MATLAB:MatFile:UnsupportedCharacters',...
        ['Found characters the default encoding is unable to represent.\n',...
         'Trying again using default version for MAT-files.']);
        save(fullfile(this), 'D');
    else
        rethrow(err);
    end
end

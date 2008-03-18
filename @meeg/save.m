function res = save(this)
% save an meeg object into a file
% FORMAT res = save(this)
%_______________________________________________________________________
%
% Converts an meeg object to struct and saves it.
%_______________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% $Id: save.m 1224 2008-03-18 12:19:15Z vladimir $


D = struct(this);
D = rmfield(D, 'cache');

res = 1;

try
    save(fullfile(D.path, D.fname), 'D');
catch
    [filename, pathname] = uiputfile('*.mat', 'Select a file to save');
    try
        save(fullfile(pathname, filename), 'D');
    catch
        res = 0;
    end
end
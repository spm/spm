function res = save(this)
% save an meeg object into a file
% FORMAT res = save(this)
%
% Converts an meeg object to struct and saves it.
%_______________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% Vladimir Latvik
% $Id: save.m 1373 2008-04-11 14:24:03Z spm $


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

function res = path(obj)
% Method for getting path
% FORMAT res = path(obj)
% _______________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% Stefan Kiebel

try
    res = obj.path;
catch
    res = '';
end
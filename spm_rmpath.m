function varargout = spm_rmpath(d)
% Recursively removes SPM paths from the MATLAB path
%   SPM_RMPATH checks if the file spm.m is found and removes the
%   path to that file and any subdirectories below it from the MATLAB
%   path.
%
%   P = SPM_RMPATH performs the same function as above and returns the
%   cleaned path string in P.
%
%   SPM_RMPATH(D) strips the path string D from the MATLAB path.
%
%   P = SPM_RMPATH(D) strips the path string D from the MATLAB path and
%   returns the cleaned path string in P.
%
%   See also PATH, ADDPATH, RMPATH, GENPATH, PATHTOOL, SAVEPATH.
%__________________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% Darren Gitelman & Guillaume Flandin
% $Id: spm_rmpath.m 4025 2010-07-29 11:10:15Z guillaume $ 

varargout = {};
if ~nargin
    % Get the actual SPM directory
    try, d = spm('dir');
    catch, return; end
end

% Recursively remove directories in the MATLAB path
p = textscan(path,'%s','delimiter',pathsep); p = p{1};
i = strncmp(d,p,length(d)); P = p(i); p(i) = [];
if ~nargin && ~isempty(P)
    fprintf('Removed %s paths starting from base path: "%s"\n',spm('ver','',1),d);
elseif ~isempty(P)
    fprintf('Removed paths starting from base path: "%s" from:\n',d);
else
    fprintf('No matching path strings found to remove\n')
end
if numel(P), fprintf('\t%s\n',P{:}); end

% Set the new MATLAB path
p = strcat(p,pathsep);
path(strcat(p{:}));

% Return the cleaned path if requested
if nargout
    varargout{1} = p;
end

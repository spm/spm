function pathparts = cfg_getfile_pathparts(inpaths)
% Split paths in cellstr inpaths into cellstr lists of directories
% FORMAT pathparts = cfg_getfile_pathparts(inpaths)
% inpaths          - cellstr list of pathnames
% pathparts        - cell array of cellstr list, containing the directory
%                    names of each inpath path
% returns cell array of path component cellstr arrays
% For PC (WIN) targets, both '\' and '/' are accepted as filesep, similar
% to MATLAB fileparts
%____________________________________________________________________________
% Copyright (C) 2005 Wellcome Department of Imaging Neuroscience
%
% This code is part of a batch job configuration system for MATLAB. See 
%      help matlabbatch
% for a general overview.
%_______________________________________________________________________
% Copyright (C) 2007 Freiburg Brain Imaging

% John Ashburner and Volkmar Glauche
% $Id: cfg_getfile.m 7067 2017-04-20 16:49:53Z guillaume $

if ispc
    fs = '\\/';
else
    fs = filesep;
end
pathparts = cellfun(@(p1)textscan(p1,'%s','delimiter',fs,'MultipleDelimsAsOne',1),inpaths);
% Canonicalise first entry for Windows drive letters and UNC paths
if ispc
    for k = 1:numel(pathparts)
        if ~isempty(regexp(pathparts{k}{1}, '^[a-zA-Z]:$', 'once'))
            pathparts{k}{1} = strcat(pathparts{k}{1}, filesep);
        elseif ~isempty(regexp(inpaths{k}, '^\\\\', 'once'))
            pathparts{k}{1} = strcat(filesep, filesep, pathparts{k}{1});
        end
    end
end

function cpaths = cfg_getfile_cpath(inpaths, basedir)
% Canonicalise paths
% FORMAT cpaths = cfg_getfile_cpath(inpaths, basedir)
% inpaths       - cell array of (relative or absolute) paths
% basedir       - single cell containing the base path of relative input paths
% cpaths        - cell array of canonical paths, one for each element of
%                 inpaths
% This function canonicalises paths by
% - removing /./ and /../ constructs
% - prepending basedir to relative pathnames
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

if ispc % valid absolute paths
    % Allow drive letter or UNC path
    mch = '^([a-zA-Z]:)|(\\\\[^\\]*)';
else
    mch = '^/';
end
if (nargin<2)||isempty(basedir), basedir = {pwd}; end
% Find partial paths, prepend them with d
ppsel    = cellfun(@isempty, regexp(inpaths,mch,'once'));
inpaths(ppsel) = cellfun(@(t1)fullfile(basedir{1},t1),inpaths(ppsel),'UniformOutput',false);
% Break paths into cell lists of folder names
pt = cfg_getfile_pathparts(inpaths);
% Remove single '.' folder names
sd = cellfun(@(pt1)strcmp(pt1,'.'),pt,'UniformOutput',false);
for cp = 1:numel(pt)
    pt{cp} = pt{cp}(~sd{cp});
end
% Go up one level for '..' folders, don't remove drive letter/server name
% from PC path
if ispc
    ptstart = 2;
else
    ptstart = 1;
end
for cp = 1:numel(pt)
    tmppt = {};
    for cdir = ptstart:numel(pt{cp})
        if strcmp(pt{cp}{cdir},'..')
            tmppt = tmppt(1:end-1);
        else
            tmppt{end+1} = pt{cp}{cdir};
        end
    end
    if ispc
        pt{cp} = [pt{cp}(1) tmppt];
    else
        pt{cp} = tmppt;
    end
end
% Assemble paths
if ispc
    cpaths = cellfun(@(pt1)fullfile(pt1{:}),pt,'UniformOutput',false);
else
    cpaths = cellfun(@(pt1)fullfile(filesep,pt1{:}),pt,'UniformOutput',false);
end

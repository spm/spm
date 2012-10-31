function p = fullfile(this)
% Returns full path to the meeg mat file
% FORMAT p = fullfile(this)
% _______________________________________________________________________
% Copyright (C) 2008-2011 Wellcome Trust Centre for Neuroimaging

% Vladimir Litvak
% $Id: fullfile.m 5025 2012-10-31 14:44:13Z vladimir $


p = fullfile(path(this), fname(this));
function p = fullfile(this)
% Returns full path to the meeg mat file
% FORMAT p = fullfile(this)
%__________________________________________________________________________

% Vladimir Litvak
% Copyright (C) 2008-2022 Wellcome Centre for Human Neuroimaging


p = fullfile(path(this), fname(this));

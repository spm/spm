function this = reload(this)
% Reload the file from disk
% FORMAT this = reload(this)
%
% Useful to update the object e.g. after running a batch.
%__________________________________________________________________________

% Vladimir Litvak
% Copyright (C) 2008-2022 Wellcome Centre for Human Neuroimaging


this = meeg(getfield(load(fullfile(this)), 'D'));

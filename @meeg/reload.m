function this = reload(this)
% Reload the file from disk
% FORMAT this = reload(this)
%
% Useful to update the object e.g. after running a batch.
%__________________________________________________________________________
% Copyright (C) 2013 Wellcome Trust Centre for Neuroimaging

% Vladimir Litvak
% $Id: reload.m 5385 2013-04-03 16:24:05Z vladimir $

this = meeg(getfield(load(fullfile(this)), 'D'));

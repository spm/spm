function this = save(this)
% Save an meeg object into a file
% FORMAT this = save(this)
%
% Converts an meeg object to struct and saves it.
%__________________________________________________________________________
% Copyright (C) 2008-2012 Wellcome Trust Centre for Neuroimaging

% Vladimir Litvak
% $Id: save.m 5078 2012-11-25 15:08:05Z vladimir $

D = struct(this);
save(fullfile(this), 'D', spm_get_defaults('mat.format'));

function save(this)
% Save an meeg object into a file
% FORMAT save(this)
%
% Converts an meeg object to struct and saves it.
%__________________________________________________________________________
% Copyright (C) 2008-2012 Wellcome Trust Centre for Neuroimaging

% Vladimir Litvak
% $Id: save.m 5025 2012-10-31 14:44:13Z vladimir $

D = struct(this);
save(fullfile(this), 'D');

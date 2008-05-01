function [ok, this] = check(this, option)
% Method that performs integrity checks of the meeg object
% and its readiness for particular purposes.
% FORMAT  this = check(this, option)
% IN
% option - 'basic' (default) - just check the essential fields
%          'sensfid' - also checks sensor and fiducial definitions
% OUT
% ok - 1 - OK, 0- failed
% _______________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% Vladimir Litvak
% $Id: check.m 1531 2008-05-01 14:17:54Z vladimir $

[ok, this] = checkmeeg(struct(this), option);

this = meeg(this);
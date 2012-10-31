function res = isequal(this, that)
% Method to check if 2 MEEG objects are the same
% FORMAT res = isequal(this, that)
% _______________________________________________________________________
% Copyright (C) 2011 Wellcome Trust Centre for Neuroimaging

% Christophe Phillips
% $Id: isequal.m 5025 2012-10-31 14:44:13Z vladimir $

res = isequal(struct(this), struct(that));
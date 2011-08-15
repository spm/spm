function res = isequal(this, that)
% Method to check if 2 MEEG objects are the same
% FORMAT res = isequal(this, that)
% _______________________________________________________________________
% Copyright (C) 2011 Wellcome Trust Centre for Neuroimaging

% Christophe Phillips
% $Id: $

this = struct(this);
that = struct(that);
res = isequal(this,that);
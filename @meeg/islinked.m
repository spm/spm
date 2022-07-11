function out = islinked(this)
% True if the object is linked to data file
% FORMAT out = islinked(this)
%__________________________________________________________________________

% Vladimir Litvak
% Copyright (C) 2008-2022 Wellcome Centre for Human Neuroimaging


out  = isa(this.data, 'file_array');

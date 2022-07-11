function res = getfield(this, varargin)
% Returns  fields in .other
% FORMAT res = getfield(this, varargin)
%
% An overloaded function...
%__________________________________________________________________________

% Vladimir Litvak
% Copyright (C) 2008-2022 Wellcome Centre for Human Neuroimaging


res = getfield(this.other, varargin{:});

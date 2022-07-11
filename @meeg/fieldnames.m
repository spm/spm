function res = fieldnames(this, varargin)
% Returns names of the fields in .other
% FORMAT res = fieldnames(this)
%
% An overloaded function...
%__________________________________________________________________________

% Vladimir Litvak
% Copyright (C) 2008-2022 Wellcome Centre for Human Neuroimaging


res = fieldnames(this.other, varargin{:});

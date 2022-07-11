function res = isfield(this, varargin)
% Returns true if the string fieldname is the name of a field in the 
% substructure 'other' in the meeg object 'this'.
% FORMAT res = isfield(this,fieldname)
%
% An overloaded function...
%__________________________________________________________________________

% Vladimir Litvak
% Copyright (C) 2008-2022 Wellcome Centre for Human Neuroimaging


res = isfield(this.other, varargin{:});

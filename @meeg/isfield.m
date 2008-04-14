function res = isfield(this,fieldname)
% Returns true if the string fieldname is the name of a field in the 
% substructure 'other' in the meeg object 'this'.
% FORMAT res = isfield(this,fieldname)
%
% An overloaded function...
% _______________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% Christophe Phillips
% $$

if isfield(this.other,fieldname)
    res = true;
else
    res = false;
end
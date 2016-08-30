function json = spm_jsonread(filename)
% JSON (JavaScript Object Notation) parser - a compiled routine
% FORMAT json = spm_jsonread(filename)
% filename - name of a JSON file or JSON string
% json     - JSON structure
% 
% References:
%   JSON Standard: http://www.json.org/
%   JSMN C parser: http://zserge.com/jsmn.html
%   jsondecode: http://www.mathworks.com/help/matlab/ref/jsondecode.html
%__________________________________________________________________________
% Copyright (C) 2015-2016 Wellcome Trust Centre for Neuroimaging

% Guillaume Flandin
% $Id: spm_jsonread.m 6863 2016-08-30 14:56:27Z guillaume $


%-This is merely the help file for the compiled routine
error('spm_jsonread.c not compiled - see Makefile')

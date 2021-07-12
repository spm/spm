function json = spm_jsonread(filename, opts)
% JSON (JavaScript Object Notation) parser - a compiled routine
% FORMAT json = spm_jsonread(filename, opts)
% filename - name of a JSON file or JSON string
% json     - JSON structure
% opts     - structure or list of name/value pairs of optional parameters:
%              ReplacementStyle: string to control how non-alphanumeric
%                characters are replaced {'underscore','hex','delete','nop'}
%                [Default: 'underscore']
%              Prefix: string to prepend when first character of a field is
%                not alphabetical [Default: 'x']
% 
% References:
%   JSON Standard: https://www.json.org/
%   JSMN C parser: https://zserge.com/jsmn/
%   jsondecode: https://www.mathworks.com/help/matlab/ref/jsondecode.html
%__________________________________________________________________________
% Copyright (C) 2015-2021 Wellcome Trust Centre for Neuroimaging

% Guillaume Flandin
% $Id: spm_jsonread.m 8124 2021-07-12 16:26:10Z guillaume $


%-This is merely the help file for the compiled routine
error('spm_jsonread.c not compiled - see Makefile')

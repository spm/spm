function str = gencode_substruct(subs, name)

% GENCODE_SUBSTRUCT  String representation of subscript structure.
% Generate MATLAB code equivalent to subscript structure subs. See help
% on SUBSTRUCT, SUBSASGN and SUBSREF for details how subscript structures
% are used.
%
% str = gencode_substruct(subs, name)
% Input arguments:
%  subs - a subscript structure
%  name - optional: name of variable to be dereferenced
% Output arguments:
%  str  - a one-line cellstr containing a string representation of the
%         subscript structure
% If name is given, it is prepended to the string.
% For '()' and '{}' also pseudo subscripts are allowed: if subs.subs{...}
% is a string, it will be printed literally, even if it is not equal to
% ':'. This way, it is possible create code snippets that contain
% e.g. references to a loop variable by name.
%
% See also GENCODE, GENCODE_RVALUE, GENCODE_SUBSTRUCTCODE.
%
% This code has been developed as part of a batch job configuration
% system for MATLAB. See  
%      http://sourceforge.net/projects/matlabbatch
% for details about the original project.
%_______________________________________________________________________
% Copyright (C) 2007 Freiburg Brain Imaging

% Volkmar Glauche
% $Id: gencode_substruct.m 3355 2009-09-04 09:37:35Z volkmar $

rev = '$Rev: 3355 $'; %#ok

ind = 1;
if nargin < 2
    name = '';
end

if ~isstruct(subs) || ~all(isfield(subs, {'type','subs'}))
    if any(exist('cfg_message') == 2:6)
        cfg_message('matlabbatch:usage', 'Item is not a substruct.');
    else
        warning('gencode_substruct:usage', 'Item is not a substruct.');
    end
else
    str = {name};
    for k = 1:numel(subs)
        switch subs(k).type
            case '.',
                str{1} = sprintf('%s.%s', str{1}, subs(k).subs);
            case {'()','{}'},
                str{1} = sprintf('%s%s', str{1}, subs(k).type(1));
                for l = 1:numel(subs(k).subs)
                    if ischar(subs(k).subs{l})
                        substr = subs(k).subs{l};
                    else
                        substr = sprintf('%d ', subs(k).subs{l});
                        if numel(subs(k).subs{l}) > 1
                            substr = sprintf('[%s]', substr(1:end-1));
                        else
                            substr = substr(1:end-1);
                        end
                    end
                    str{1} = sprintf('%s%s, ', str{1}, substr);
                end
                str{1} = str{1}(1:end-2);
                str{1} = sprintf('%s%s', str{1}, subs(k).type(2));
        end
    end
end

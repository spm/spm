function str = gencode_substructcode(subs, name)

% GENCODE_SUBSTRUCTCODE  Create code for a subscript structure
% Generate MATLAB code (using SUBSTRUCT) to create subscript structure
% subs. See help on SUBSTRUCT, SUBSASGN and SUBSREF for details how
% subscript structures are used.
%
% str = gencode_substructcode(subs, name)
% Input arguments:
%  subs - a subscript structure
%  name - optional: name of variable
% Output arguments:
%  str  - a one-line cellstr containing a call to SUBSTRUCT that returns
%         an substruct equivalent to subs.
% If name is supplied as input argument, the generated code will assign
% the output of SUBSTRUCT to the variable 'name'.
% then only the rhs of the expression will be returned.
% For '()' and '{}' also pseudo subscripts are allowed: if subs.subs{...}
% is a string, it will be printed literally, even if it is not equal to
% ':'. This way, one can create code snippets that contain e.g. references
% to a loop variable by name.
%
% See also GENCODE, GENCODE_RVALUE, GENCODE_SUBSTRUCT.
%
% This code has been developed as part of a batch job configuration
% system for MATLAB. See  
%      http://sourceforge.net/projects/matlabbatch
% for details about the original project.
%_______________________________________________________________________
% Copyright (C) 2007 Freiburg Brain Imaging

% Volkmar Glauche
% $Id: gencode_substructcode.m 3355 2009-09-04 09:37:35Z volkmar $

rev = '$Rev: 3355 $'; %#ok

ind = 1;
if nargin < 2
    name = '';
end

if ~isstruct(subs) || ~all(isfield(subs, {'type','subs'}))
    if any(exist('cfg_message') == 2:6)
        cfg_message('matlabbatch:usage', 'Item is not a substruct.');
    else
        warning('gencode_substructcode:usage', ...
                'Item is not a substruct.');
    end
else
    if isempty(subs)
        str = {'struct(''type'',{},''subs'',{})'};
    else
        str = {'substruct('};
        for k = 1:numel(subs)
            str{1} = sprintf('%s''%s'',', str{1}, subs(k).type);
            switch subs(k).type
                case '.',
                    substr = sprintf('''%s''', subs(k).subs);
                case {'()','{}'},
                    substr = '{';
                    for l = 1:numel(subs(k).subs)
                        if ischar(subs(k).subs{l})
                            if strcmp(subs(k).subs{l},':')
                                substr = sprintf('%s'':'', ', substr);
                            else
                                substr = sprintf('%s%s, ', substr, subs(k).subs{l});
                            end
                        else
                            substr1 = sprintf('%d ', subs(k).subs{l});
                            if numel(subs(k).subs{l}) > 1
                                substr1 = sprintf('[%s]', substr1(1:end-1));
                            else
                                substr1 = substr1(1:end-1);
                            end
                            substr = sprintf('%s%s, ', substr, substr1);
                        end
                    end
                    substr = sprintf('%s}', substr(1:end-2));
            end
            str{1} = sprintf('%s%s, ', str{1}, substr);
        end
        str{1} = sprintf('%s)', str{1}(1:end-2));
    end
    if ~isempty(name)
        str{1} = sprintf('%s = %s;', name, str{1});
    end
end

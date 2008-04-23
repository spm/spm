function [str, name, ind] = gencode_substructcode(subs, name)

% Generate MATLAB code that recreates a given subscript structure using
% a substruct call. 
% str is a 1-line cellstr, name is the used name and ind is 1 (to conform
% to general gencode calling syntax). If name is not supplied or empty,
% then only the rhs of the expression will be returned.
% For '()' and '{}' also pseudo subscripts are allowed: if subs.subs{...}
% is a string, it will be printed literally, even if it is not equal to
% ':'. This way, one can create code snippets that contain e.g. references
% to a loop variable by name.
%_______________________________________________________________________
% Copyright (C) 2007 Freiburg Brain Imaging

% Volkmar Glauche
% $Id: gencode_substructcode.m 1472 2008-04-23 12:02:44Z volkmar $

rev = '$Rev: 1472 $';

ind = 1;
if nargin < 2
    name = '';
end;

if ~isstruct(subs) || ~all(isfield(subs, {'type','subs'}))
    error('gencode:substruct:notsubs', 'Item is not a substruct.');
end;
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
                        end;
                    else
                        substr1 = sprintf('%d ', subs(k).subs{l});
                        if numel(subs(k).subs{l}) > 1
                            substr1 = sprintf('[%s]', substr1(1:end-1));
                        else
                            substr1 = substr1(1:end-1);
                        end;
                        substr = sprintf('%s%s, ', substr, substr1);
                    end;
                end
                substr = sprintf('%s}', substr(1:end-2));
        end;
        str{1} = sprintf('%s%s, ', str{1}, substr);
    end;
    str{1} = sprintf('%s)', str{1}(1:end-2));
end;
if ~isempty(name)
    str{1} = sprintf('%s = %s;', name, str{1});
end;
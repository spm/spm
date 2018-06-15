function [str, sts] = gencode_rvalue(item, cflag)

% GENCODE_RVALUE  Code for right hand side of MATLAB assignment
% Generate the right hand side for a valid MATLAB variable
% assignment. This function is a helper to GENCODE, but can be used on
% its own to generate code for the following types of variables:
% * scalar, 1D or 2D numeric, logical or char arrays
% * scalar or 1D cell arrays, where each item can be one of the supported
%   array types (i.e. nested cells are allowed)
%
% function [str, sts] = gencode_rvalue(item, cflag)
% Input argument:
%  item  - value to generate code for
%  cflag - (optional) if true, try to generate 1-line code also for 2D
%          arrays. This may reduce readability of the generated code.
%          Defaults to false.
% Output arguments:
%  str - cellstr with generated code, line per line
%  sts - true, if successful, false if code could not be generated
%
% See also GENCODE, GENCODE_SUBSTRUCT, GENCODE_SUBSTRUCTCODE.
%
% This code has been developed as part of a batch job configuration
% system for MATLAB. See  
%      http://sourceforge.net/projects/matlabbatch
% for details about the original project.
%_______________________________________________________________________
% Copyright (C) 2007 Freiburg Brain Imaging

% Volkmar Glauche
% $Id: gencode_rvalue.m 7344 2018-06-15 12:44:46Z volkmar $

rev = '$Rev: 7344 $'; %#ok

if nargin < 2
    cflag = false;
end

str = {};
sts = true;
switch class(item)
    case 'char'
        if ndims(item) == 2 %#ok<ISMAT>
            if isempty(item)
                if all(size(item)==0)
                    str1 = {''''''};
                else
                    str1 = {sprintf('char(zeros%s)', char(gencode_substruct(substruct('()', num2cell(size(item))))));};
                end
            else
                % Create cell string, keep white space padding
                cstr = cell(1,size(item,1));
                for k = 1:size(item,1)
                    cstr{k} = item(k,:);
                end
                str1 = genstrarray(cstr);
            end
            if numel(str1) == 1
                % One string, do not print brackets
                str = str1;
            else
                % String array, print brackets and concatenate str1
                if cflag
                    str = {['[' sprintf('%s;', str1{:}) ']']};
                else
                    str = [ {'['} str1(:)' {']'} ];
                end
            end
        else
            % not an rvalue
            sts = false;
        end
    case 'cell'
        if isempty(item)
            if all(size(item)==0)
                str = {'{}'};
            else
                str = {sprintf('cell%s', char(gencode_substruct(substruct('()', num2cell(size(item))))))};
            end
        elseif ndims(item) == 2 && any(size(item) == 1) %#ok<ISMAT>
            str1 = {};
            for k = 1:numel(item)
                [str2, sts] = gencode_rvalue(item{k}, cflag);
                if ~sts
                    break;
                end
                str1 = [str1(:)' str2(:)'];
            end
            if sts
                if numel(str1) == 1
                    % One item, print as one line
                    str{1} = sprintf('{%s}', str1{1});
                else
                    % Cell vector, print braces and concatenate str1
                    if size(item,1) == 1
                        endstr = '}''';
                    else
                        endstr = '}';
                    end
                    if cflag
                        str = {['{' sprintf('%s;', str1{:}) endstr]};
                    else
                        str = [{'{'} str1(:)' {endstr}];
                    end
                end
            end
        else
            sts = false;
        end
    case 'function_handle'
        fstr = func2str(item);
        % sometimes (e.g. for anonymous functions) '@' is already included
        % in func2str output
        if fstr(1) == '@'
            str{1} = fstr;
        else
            str{1} = sprintf('@%s', fstr);
        end
    otherwise
        if isobject(item) || ~(isnumeric(item) || islogical(item)) || issparse(item) || ndims(item) > 2 %#ok<ISMAT>
            sts = false;
        else
            % treat item as numeric or logical, don't create 'class'(...)
            % classifier code for double
            clsitem = class(item);
            if isempty(item)
                if all(size(item)==0)
                    tmpstr = '[]';
                else
                    tmpstr = sprintf('zeros%s', char(gencode_substruct(substruct('()', num2cell(size(item))))));
                end
                if strcmp(clsitem, 'double') %#ok<STISA>
                    str{1} = tmpstr;
                else
                    str{1} = sprintf('%s(%s)', clsitem, tmpstr);
                end
            else
                % Use mat2str with standard precision 15
                if any(strcmp(clsitem, {'double', 'logical'}))
                    sitem = mat2str(item);
                else
                    sitem = mat2str(item,'class');
                end
                if cflag
                    str = {sitem};
                else
                    try
                        if ~verLessThan('matlab', '8.4')
                            bszopt = {};
                        else
                            error('Need bufsize option');
                        end
                    catch
                        bsz   = max(numel(sitem)+2,100); % bsz needs to be > 100 and larger than string length
                        bszopt = {'bufsize', bsz};
                    end
                    str1 = textscan(sitem, '%s', 'delimiter',';', bszopt{:});
                    if numel(str1{1}) > 1
                        str = str1{1};
                    else
                        str{1} = str1{1}{1};
                    end
                end
            end
        end
end
function str = genstrarray(stritem)
% generate a cell string of properly quoted strings suitable for code
% generation.
str = strrep(stritem, '''', '''''');
for k = 1:numel(str)
    if ~any(str{k} == char(0)) &&  ~any(str{k} == char(9)) && ~any(str{k} == char(10))
        str{k} = sprintf('''%s''', str{k});
    else
        % first, quote sprintf special chars % and \
        % second, replace special characters by sprintf equivalents
        replacements = {'%', '%%'; ...
            '\', '\\'; ...
            char(0), '\0'; ...
            char(9), '\t'; ...
            char(10), '\n'; ...
            char(13), '\r'};
        for cr = 1:size(replacements, 1)
            str{k} = strrep(str{k}, replacements{cr,:});
        end
        str{k} = sprintf('sprintf(''%s'')', str{k});
    end
end


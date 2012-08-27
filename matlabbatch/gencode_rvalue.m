function [str, sts] = gencode_rvalue(item)

% GENCODE_RVALUE  Code for right hand side of MATLAB assignment
% Generate the right hand side for a valid MATLAB variable
% assignment. This function is a helper to GENCODE, but can be used on
% its own to generate code for the following types of variables:
% * scalar, 1D or 2D numeric, logical or char arrays
% * scalar or 1D cell arrays, where each item can be one of the supported
%   array types (i.e. nested cells are allowed)
%
% function [str, sts] = gencode_rvalue(item)
% Input argument:
%  item - value to generate code for
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
% $Id: gencode_rvalue.m 4864 2012-08-27 13:57:31Z volkmar $

rev = '$Rev: 4864 $'; %#ok

str = {};
sts = true;
switch class(item)
    case 'char'
        if ndims(item) == 2 %#ok<ISMAT>
            cstr = {''};
            % Create cell string, keep white space padding
            for k = 1:size(item,1)
                cstr{k} = item(k,:);
            end
            str1 = genstrarray(cstr);
            if numel(str1) == 1
                % One string, do not print brackets
                str = str1;
            else
                % String array, print brackets and concatenate str1
                str = [ {'['} str1(:)' {']'} ];
            end
        else
            % not an rvalue
            sts = false;
        end
    case 'cell'
        if isempty(item)
            str = {'{}'};
        elseif ndims(item) == 2 && any(size(item) == 1) %#ok<ISMAT>
            str1 = {};
            for k = 1:numel(item)
                [str2 sts] = gencode_rvalue(item{k});
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
                        endstr = {'}'''};
                    else
                        endstr = {'}'};
                    end
                    str = [{'{'} str1(:)' endstr];
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
                if strcmp(clsitem, 'double') %#ok<STISA>
                    str{1} = '[]';
                else
                    str{1} = sprintf('%s([])', clsitem);
                end
            else
                % Use mat2str with standard precision 15
                if any(strcmp(clsitem, {'double', 'logical'}))
                    sitem = mat2str(item);
                else
                    sitem = mat2str(item,'class');
                end
                bsz   = max(numel(sitem)+2,100); % bsz needs to be > 100 and larger than string length
                str1 = textscan(sitem, '%s', 'delimiter',';', 'bufsize',bsz); 
                if numel(str1{1}) > 1
                    str = str1{1};
                else
                    str{1} = str1{1}{1};
                end
            end
        end
end
function str = genstrarray(stritem)
% generate a cell string of properly quoted strings suitable for code
% generation.
str = strrep(stritem, '''', '''''');
for k = 1:numel(str)
    str{k} = sprintf('''%s''', str{k});
end


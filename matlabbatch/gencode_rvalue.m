function [str, sts] = gencode_rvalue(item)
str = {};
sts = true;
switch class(item)
    case 'char'
        if ndims(item) == 2
            cstr = {''};
            % Create cell string, keep white space padding
            for k = 1:size(item,1)
                cstr{k} = item(k,:);
            end;
            str1 = genstrarray(cstr);
            if numel(str1) == 1
                % One string, do not print brackets
                str = str1;
            else
                % String array, print brackets and concatenate str1
                str = {'[' str1{:} ']'};
            end;
        else
            % not an rvalue
            sts = false;
        end;
    case 'cell'
        if isempty(item)
            str = {'{}'};
        elseif ndims(item) == 2 && any(size(item) == 1)
            str1 = {};
            for k = 1:numel(item)
                [str2 sts] = gencode_rvalue(item{k});
                if ~sts
                    break;
                end;
                str1 = {str1{:} str2{:}};
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
                    end;
                    str = {'{' str1{:} endstr};
                end;
            end;
        else
            sts = false;
        end;
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
        if isobject(item) || ~(isnumeric(item) || islogical(item)) || ndims(item) > 2
            sts = false;
        else
            % treat item as numeric or logical, don't create 'class'(...)
            % classifier code for double
            clsitem = class(item);
            if isempty(item)
                if strcmp(clsitem, 'double')
                    str{1} = '[]';
                else
                    str{1} = sprintf('%s([])', clsitem);
                end
            else
                % Use mat2str with standard precision 15
                if any(strcmp(clsitem, {'double', 'logical'}))
                    str1 = textscan(mat2str(item), '%s', 'delimiter', ';');
                else
                    str1 = textscan(mat2str(item,'class'), '%s', 'delimiter', ';');
                end;
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


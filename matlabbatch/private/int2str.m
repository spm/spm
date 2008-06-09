function s = int2str(x)
%INT2STR Convert integer to string.
%   S = INT2STR(X) rounds the elements of the matrix X to
%   integers and converts the result into a string matrix.
%   Return NaN and Inf elements as strings 'NaN' and 'Inf', respectively.
%
%   Modified by Volkmar Glauche to return 'true' and 'false' instead of 0
%   and 1 for logical arrays.
%
%   See also NUM2STR, SPRINTF, FPRINTF, MAT2STR.

%   Copyright 1984-2006 The MathWorks, Inc.
%   $Revision: 284 $  $Date: 2006/11/11 22:45:05 $

if ~islogical(x)
    x = round(real(x));
end;
if length(x) == 1
    if islogical(x) 
        % handle logical
        if x
            s = 'true';
        else
            s = 'false';
        end;
    elseif isfinite(x)
        % handle special case of single infinite or NaN element
        s = sprintf('%.1f',x); % .1 to avoid precision loss on hpux
        s = s(1:length(s)-2);  % remove .d decimal digit
    else
        s = sprintf('%.0f',x);
    end
else

    t = '';
    [m,n] = size(x);
    if islogical(x)
        scell = cell(m,n);
        [scell{x}] = deal('true  ');
        [scell{~x}] = deal('false ');
        for k = 1:m
            s(k,:) = [scell{k,:}];
        end;
    else
        % Determine elements of x that are finite.
        xfinite = x(isfinite(x));
        % determine the variable text field width quantity
        
        xmax = double(max(abs(xfinite(:))));
        if xmax == 0
            d = 3;
        else
            d = floor(log10(xmax)) + 3;
        end

    
        % delimit string array with one space between all-NaN or all-Inf columns
        if numel(xfinite) ~= numel(x)
            d = max([d;5]);
        end
    
        clear('xfinite')
    
        % create cell arrays for storage
        scell = cell(1,m);
        empties = true(1,m);

        precision_string = sprintf('%%%d.0f', d);
        % walk through numbers array and convert elements to strings
        for i = 1:m
            % use vectorized version of sprintf
            t = sprintf(precision_string,x(i,:));
            if ~isempty(t)
                scell{i} = t;
                empties(i) = false;
            end
        end
        if m > 1
            s = char(scell{~empties});
        else
            s = t;
        end
        % trim leading spaces from string array within constraints of rectangularity.
        if ~isempty(s)
            while all(s(:,1) == ' ')
                s(:,1) = [];
            end
        end
    end;
end

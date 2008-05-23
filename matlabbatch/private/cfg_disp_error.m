function disp_error(l)

% function disp_error(errstruct)
%
% Display a condensed version of a MATLAB error without rethrowing it.
%
% This code is part of a batch job configuration system for MATLAB. See 
%      help matlabbatch
% for a general overview.
%_______________________________________________________________________
% Copyright (C) 2007 Freiburg Brain Imaging

% Volkmar Glauche
% $Id: copyright_cfg.m 218 2008-04-17 10:34:11Z glauche $

rev = '$Rev: 218 $'; %#ok

disp(l.message);
if isfield(l,'stack'), % Does not always exist
    for m = 1:numel(l.stack),
        try
            fp  = fopen(l.stack(m).file,'r');
            str = fread(fp,Inf,'*uchar');
            fclose(fp);
            str = char(str(:)');
            re  = regexp(str,'\$Id: \w+\.\w+ ([0-9]+) [0-9][0-9][0-9][0-9].*\$','tokens');
            if numel(re)>0 && numel(re{1})>0,
                id = [' (v', re{1}{1}, ')'];
            else
                id = ' (???)';
            end
        catch
            id = '';
        end
        fprintf('In file "%s"%s, function "%s" at line %d.\n', ...
                l.stack(m).file, id, l.stack(m).name, l.stack(m).line);
    end
end;

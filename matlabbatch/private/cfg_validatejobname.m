function nname = cfg_validatejobname(name, cont)
% CFG_VALIDATEJOBNAME - validate the name of a job file
% 
% See also GENVARNAME
%
% This code is part of a batch job configuration system for MATLAB. See 
%      help matlabbatch
% for a general overview.
%_______________________________________________________________________
% Copyright (C) 2007 Freiburg Brain Imaging

% Volkmar Glauche
% $Id: cfg_validatejobname.m 2673 2009-01-30 13:34:53Z volkmar $

rev = '$Rev: 2673 $'; %#ok

nname = genvarname(name);
if ~strcmp(nname, name)
    if cont
        cfg_message('matlabbatch:validatejobname:soft',...
                    '''%s'' is not a valid name for a MATLAB script file, using ''%s'' instead.',...
                    name, nname);
    else
        cfg_message('matlabbatch:validatejobname:hard',...
                    ['''%s'' is not a valid name for a MATLAB script file.\n',...
                     'Allowed characters are\n',...
                     '* letters,\n',...
                     '* numbers and\n',...
                     '* underscores.\n',...
                     'The first character must be a letter.\n',...
                     'See MATLAB documentation for details.'],...
                    name);
    end
end
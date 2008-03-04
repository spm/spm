function sts = subsasgn_check_num(val)

% Check, whether a num value is a numeric 2-vector
%
% This code is part of a batch job configuration system for MATLAB. See 
%      help matlabbatch
% for a general overview.
%_______________________________________________________________________
% Copyright (C) 2007 Freiburg Brain Imaging

% Volkmar Glauche
% $Id: subsasgn_check_num.m 1184 2008-03-04 16:27:57Z volkmar $

rev = '$Rev: 1184 $';
if numel(val)==2 && isnumeric(val)
    sts = true;
else
    sts = false;
    warning('matlabbatch:subsasgn_check_num:notdim', 'Value of field ''num'' must be a 2-vector.');
    disp(val);
end;

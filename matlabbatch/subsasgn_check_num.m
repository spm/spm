function sts = subsasgn_check_num(val)

% function sts = subsasgn_check_num(val)
% Check, whether a num value is a numeric 2-vector, denoting a
% minimum/maximum number of elements. val(1) must be >= 0 and 
% val(2) >= val(1). 
% This function is called for all num fields, except those in cfg_entry
% objects. 
%
% This code is part of a batch job configuration system for MATLAB. See 
%      help matlabbatch
% for a general overview.
%_______________________________________________________________________
% Copyright (C) 2007 Freiburg Brain Imaging

% Volkmar Glauche
% $Id: subsasgn_check_num.m 1862 2008-06-30 14:12:49Z volkmar $

rev = '$Rev: 1862 $'; %#ok

sts = numel(val)==2 && isnumeric(val) && val(1)>=0 && val(2)>=val(1);
if ~sts
    cfg_message('matlabbatch:check:num', ...
            ['Value of field ''num'' must be a 2-vector with non-negative ' ...
             'equal or ascending elements.\n%s'], evalc('disp(val)'));
end

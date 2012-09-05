function njob = cfg_basicio_rewrite(ojob)
% Rewrite job to conform to new submenu structure of BasicIO
%
% This code is part of a batch job configuration system for MATLAB. See 
%      help matlabbatch
% for a general overview.
%_______________________________________________________________________
% Copyright (C) 2007 Freiburg Brain Imaging

% Volkmar Glauche
% $Id: cfg_basicio_rewrite.m 4899 2012-09-05 13:44:17Z volkmar $

rev = '$Rev: 4899 $'; 

% new job layout has fields file_dir, var_ops, run_ops
switch char(fieldnames(ojob))
    case {'file_dir', 'var_ops', 'run_ops'}
        njob = ojob;
    case {'file_move', 'cfg_named_file', 'file_fplist', 'file_filter', 'cfg_file_split'}
        njob.file_dir.file_ops = ojob;
    case {'cfg_cd' 'cfg_mkdir' 'cfg_named_dir'}
        njob.file_dir.dir_ops = ojob;
    case {'cfg_fileparts'}
        njob.file_dir = ojob;
    case {'cfg_named_input', 'load_vars', 'cfg_save_vars', 'subsrefvar', 'cfg_assignin'}
        njob.var_ops = ojob;
    case {'runjobs', 'call_matlab'}
        njob.run_ops = ojob;
end
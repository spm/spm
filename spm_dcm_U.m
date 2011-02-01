function spm_dcm_U(DCM_filename,SPM_filename,session,input_nos)
% Insert new inputs into a DCM model
% FORMAT spm_dcm_U(DCM_filename,SPM_filename,session,input_nos)
%
% DCM_filename  - Name of DCM file  (char array)
% SPM_filename  - Name of SPM file  (char array)
% session       - Session number    (integer)
% input_nos     - Inputs to include (cell array)
%
% Examples of specification of parameter 'input_nos':
% * without parametric modulations:
%   {1, 0, 1} includes inputs 1 and 3.
% * with parametric modulations:
%   {1,0,[0 0 1],[0 1]} includes the non-modulated first input, the second
%   PM of the third input and the first PM of the fourth input.
% Note that this cell array only has to be specified up to the last input
% that is replaced.
%
% This function can be used, for example, to replace subject X's inputs by
% subject Y's. The model can then be re-estimated without having to go
% through model specification again.
%__________________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% Will Penny & Klaas Enno Stephan
% $Id: spm_dcm_U.m 4185 2011-02-01 18:46:18Z guillaume $


%-Load DCM and SPM files
%--------------------------------------------------------------------------
load(DCM_filename);
load(SPM_filename);


%-Get session
%--------------------------------------------------------------------------
try
    Sess = SPM.Sess(session);
catch
    error('SPM file does not have a session %d.',session);
end


%-Check numbers of inputs
%--------------------------------------------------------------------------
if size(DCM.c,2) ~= sum(cellfun(@nnz,input_nos))
    error('Number of specified inputs does not match DCM.');
end
if length(input_nos) > length(Sess.U)
    error('More inputs specified than exist in SPM.mat.');
end


%-Replace inputs
%--------------------------------------------------------------------------
U.name = {};
U.u    = [];
U.dt   = DCM.U.dt;
for i  = 1:length(input_nos)
    if any(input_nos{i})
        mo = find(input_nos{i});
        if (length(mo)-1) > length(Sess.U(i).P)
            error(['More parametric modulations specified than exist ' ...
                'for input %s in SPM.mat.'],Sess.U(i).name{1});
        end
        for j=mo
            U.u             = [U.u Sess.U(i).u(33:end,j)];
            U.name{end + 1} = Sess.U(i).name{j};
        end
    end
end
DCM.U = U;


%-Check inputs and outputs match up (to the nearest DCM.U.dt)
%--------------------------------------------------------------------------
DCM.U.dt      = Sess.U(1).dt;
DCM.Y.dt      = SPM.xY.RT;

num_inputs    = size(DCM.U.u,1);
input_period  = DCM.U.dt*num_inputs;
output_period = DCM.v*DCM.Y.dt;
if round(DCM.v*DCM.Y.dt/DCM.U.dt) ~= num_inputs
    error(sprintf(['Input period and output period do not match.\n'...
      sprintf('Number of inputs=%d, input dt=%1.2f, input period=%1.2f\n',...
        num_inputs,DCM.U.dt,input_period) ...
      sprintf('Number of outputs=%d, output dt=%1.2f, output period=%1.2f\n',...
        DCM.v,DCM.Y.dt,output_period)]));
end


% Save DCM with replaced inputs
%--------------------------------------------------------------------------
if spm_check_version('matlab','7') >= 0
    save(DCM_filename, 'DCM','-V6');
else
    save(DCM_filename, 'DCM');
end

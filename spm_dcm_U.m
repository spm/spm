function DCM = spm_dcm_U(DCM,SPM,sess,inputs)
% Insert new inputs into a DCM
% FORMAT DCM = spm_dcm_U(DCM,SPM,sess,inputs)
%
% DCM     - DCM structure or its filename
% SPM     - SPM structure or its filename
% sess    - session index     (integer)
% inputs  - Inputs to include (cell array)
%
% Examples of specification of parameter 'inputs':
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
% Copyright (C) 2003-2014 Wellcome Trust Centre for Neuroimaging

% Will Penny & Klaas Enno Stephan
% $Id: spm_dcm_U.m 7228 2017-11-21 12:17:03Z peter $


%-Load DCM and SPM files
%--------------------------------------------------------------------------
if ~isstruct(DCM)
    DCMfile = DCM;
    load(DCM);
end
if ~isstruct(SPM)
    load(SPM);
end

%-Get session
%--------------------------------------------------------------------------
try
    Sess = SPM.Sess(sess);
catch
    error('SPM file does not have a session %d.',sess);
end


%-Check numbers of inputs
%--------------------------------------------------------------------------
if size(DCM.c,2) ~= sum(cellfun(@nnz,inputs))
    error('Number of specified inputs does not match DCM.');
end
if numel(inputs) > numel(Sess.U)
    error('More inputs specified than exist in SPM.mat.');
end


%-Replace inputs
%--------------------------------------------------------------------------
U.name = {};
U.u    = [];
U.idx  = [];
try
    U.dt   = DCM.U.dt;
catch
    U.dt   = Sess.U(1).dt;
end
for i  = 1:numel(inputs)
    if any(inputs{i})
        mo = find(inputs{i});
        num_regressors = size(Sess.U(i).u,2);
        if length(mo) > num_regressors
            error(['More regressors specified than exist ' ...
                'for input %s in SPM.mat.'],Sess.U(i).name{1});
        end
        for j=mo
            U.u             = [U.u Sess.U(i).u(33:end,j)];
            U.name{end + 1} = Sess.U(i).name{j};
            U.idx           = [U.idx; i j];   
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


% Save (overwrite) DCM with replaced inputs
%--------------------------------------------------------------------------
if exist('DCMfile','var')
    save(DCMfile, 'DCM', spm_get_defaults('mat.format'));
end

function [] = spm_dcm_U (DCM_filename,SPM_filename,input_nos)
% Insert new inputs into a DCM model
% FORMAT [] = spm_dcm_U (DCM_filename,SPM_filename,input_nos)
%
% DCM_filename      Name of DCM file
% Sess              Session field from SPM.Sess
% input_nos         Inputs to include (eg. [1 3] to include inputs 1 and 3)
%
% This function can be used, for example, to replace subject X's inputs by subject Y's.
% The model can then be re-estimated without having to go through
% model specification again.
%
% %W% Will Penny %E%

load(DCM_filename);
load(SPM_filename);

Sess   = SPM.Sess(DCM.xY(1).Sess);
U.dt   = Sess.U(1).dt;

% Number of inputs in SPM file
u      = length(Sess.U);

% Number of selected inputs from SPM file
m_sel=length(input_nos);

% Number of inputs in DCM file
m = size(DCM.C,2);

if ~(m_sel==m)
    disp(sprintf('Error in spm_dcm_U: must include %d inputs',m));
    return
end

U.name = {};
U.u    = [];
for  k = 1:m_sel;
    i=input_nos(k);
    if ~(i<=u)
        disp(sprintf('Error in spm_dcm_U: input number %d not in SPM file',i));
        keyboard
        return
    end
    U.u             = [U.u Sess.U(i).u(33:end,1)];
    U.name{end + 1} = Sess.U(i).name{1};
end
DCM.U=U;

instr=['save ',DCM_filename,' DCM'];
eval(instr);



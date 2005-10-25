 function [] = spm_dcm_U (DCM_filename,SPM_filename,session,input_nos)
% Insert new inputs into a DCM model
% FORMAT [] = spm_dcm_U (DCM_filename,SPM_filename,session,input_nos)
%
% DCM_filename      Name of DCM file
% SPM_filename      Name of SPM file (eg. 'SPM')
% session           Session number (eg. 1)
% input_nos         Inputs to include (eg. [1 3] to include inputs 1 and 3)
%
% This function can be used, for example, to replace subject X's inputs by subject Y's.
% The model can then be re-estimated without having to go through
% model specification again.
%_______________________________________________________________________
% Copyright (C) 2005 Wellcome Department of Imaging Neuroscience

% Will Penny
% $Id: spm_dcm_U.m 270 2005-10-25 17:25:49Z klaas $


load(DCM_filename,'-mat');
load(SPM_filename,'-mat');

if session>length(SPM.Sess)
    disp(sprintf('Error in spm_dcm_U: SPM file doesnt have %d sessions',session));
    return
end
Sess   = SPM.Sess(session);

% Get input sampling rate from first input
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
        return
    end
    U.u             = [U.u Sess.U(i).u(33:end,1)];
    U.name{end + 1} = Sess.U(i).name{1};
end
DCM.U=U;

% Use the TR from the SPM data structure
DCM.Y.dt=SPM.xY.RT;

% Check inputs and outputs match up (to the nearest DCM.U.dt)
num_inputs=size(DCM.U.u,1);
input_period=DCM.U.dt*num_inputs;
output_period=DCM.v*DCM.Y.dt;
if ~(round(DCM.v*DCM.Y.dt/DCM.U.dt) == num_inputs)
    disp('Error in spm_dcm_U: input period and output period do not match');
    disp(sprintf('Number of inputs=%d, input dt=%1.2f, input period=%1.2f',num_inputs,DCM.U.dt,input_period));
    disp(sprintf('Number of outputs=%d, output dt=%1.2f, output period=%1.2f',DCM.v,DCM.Y.dt,output_period));
    return
end

if str2num(version('-release'))>=14,
    save(DCM_filename, 'DCM','-V6');
else
    save(DCM_filename, 'DCM');
end;


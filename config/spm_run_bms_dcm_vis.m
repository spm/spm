function out = spm_run_bms_dcm_vis(varargin)

% API to review BMS results
%
% 
% __________________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% Chun-Chuan Chen and Maria Joao Rosa
% $Id: spm_run_bms_dcm_vis.m 3445 2009-10-06 11:22:23Z maria $

job     = varargin{1};

if isempty(job.file{1})==1
    fname = spm_select([1 1],'^BMS.mat$','select BMS.mat file');
    load(fname);
else
    load(job.file{1});
end

for loop=1:1e10

    A = spm_input('Inference method','0','b',{'Fixed |Random|Bye'},[0 1 2],3);

    switch A
        
        case(0)

            if isfield(BMS.DCM,'ffx')
                N =size(BMS.DCM.ffx.F,2);
                N=1:N;
                out=spm_api_bmc(BMS.DCM.ffx.SF,N);
            else
                msgbox('This result does not exist! Please perform BMS first!');
                spm_input('Thank you',1,'d');
                return

            end

        case (1)
            if isfield(BMS.DCM,'rfx')
                N        = size(BMS.DCM.rfx.F,2);
                N        = 1:N;
                
                if isfield(BMS.DCM.rfx,'model')
                out      = spm_api_bmc(BMS.DCM.rfx.SF,...
                    N,BMS.DCM.rfx.model.exp_r,BMS.DCM.rfx.model.xp);
                else % Older version (prior to family level)
                out      = spm_api_bmc(BMS.DCM.rfx.SF,...
                    N,BMS.DCM.rfx.exp_r,BMS.DCM.rfx.xp);
                end
            else
                msgbox('This result does not exist! Please perform BMS first!');
                spm_input('Thank you',1,'d');
                return
            end

        case(2)
 
            spm_input('Thank you',1,'d');
            return

    end
   
end

  

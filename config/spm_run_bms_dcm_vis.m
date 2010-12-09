function out = spm_run_bms_dcm_vis(job)
% Review BMS results
%__________________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% Chun-Chuan Chen and Maria Joao Rosa
% $Id: spm_run_bms_dcm_vis.m 4136 2010-12-09 22:22:28Z guillaume $

if ~exist('job','var') || isempty(job.file{1})
    fname = spm_select([1 1],'^BMS\.mat$','select BMS.mat file');
    load(fname);
else
    load(job.file{1});
end

A = 1;

while A

    A = spm_input('Inference method','0','b',{'Fixed |Random|Bye'},[1 2 0],3);

    switch A
        
        case 1

            if isfield(BMS.DCM,'ffx')
                N   = size(BMS.DCM.ffx.F,2);
                N   = 1:N;
                out = spm_api_bmc(BMS.DCM.ffx.SF,N);
            else
                spm('alert*', {'This result does not exist.',...
                               'Please perform BMS first.'},'BMS');
                A   = 0;
            end

        case 2
            
            if isfield(BMS.DCM,'rfx')
                N   = size(BMS.DCM.rfx.F,2);
                N   = 1:N;
                if isfield(BMS.DCM.rfx,'model')
                    out = spm_api_bmc(BMS.DCM.rfx.SF,...
                        N,BMS.DCM.rfx.model.exp_r,BMS.DCM.rfx.model.xp);
                else % Older version (prior to family level)
                    out = spm_api_bmc(BMS.DCM.rfx.SF,...
                        N,BMS.DCM.rfx.exp_r,BMS.DCM.rfx.xp);
                end
            else
                spm('alert*', {'This result does not exist.',...
                               'Please perform BMS first.'},'BMS');
                A   = 0;
            end

    end
    
    spm_input('Thank you',1,'d');
   
end

function out = spm_run_dcm_bms_vis(job)
% Review BMS results
%__________________________________________________________________________
% Copyright (C) 2008-2014 Wellcome Trust Centre for Neuroimaging

% CC Chen and Maria Joao Rosa
% $Id: spm_run_dcm_bms_vis.m 6004 2014-05-21 14:24:14Z guillaume $


if ~nargin || isempty(job.bmsmat{1})
    [bmsmat, sts] = spm_select(1,'^BMS\.mat$','Select BMS.mat');
    if ~sts, out = []; return; end
else
    bmsmat = job.file{1};
end

try
    load(bmsmat);
catch
    error('Cannot load file: %s', bmsmat);
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

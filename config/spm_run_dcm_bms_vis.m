function out = spm_run_dcm_bms_vis(job)
% Review BMS results
%__________________________________________________________________________
% Copyright (C) 2008-2014 Wellcome Trust Centre for Neuroimaging

% CC Chen and Maria Joao Rosa
% $Id: spm_run_dcm_bms_vis.m 6005 2014-05-21 16:46:26Z guillaume $


if ~nargin || isempty(job.bmsmat{1})
    [bmsmat, sts] = spm_select(1,'^BMS\.mat$','Select BMS.mat');
    if ~sts, out = []; return; end
else
    bmsmat = job.bmsmat{1};
end

try
    load(bmsmat);
catch
    error('Cannot load file: %s', bmsmat);
end

lab = {};
val = {};
if isfield(BMS.DCM,'ffx')
    lab = [lab 'Fixed'];
    val = [val 'ffx'];
end
if isfield(BMS.DCM,'rfx')
    lab = [lab 'Random'];
    val = [val 'rfx'];
end
if numel(val) == 2
    lab = [lab 'Bye'];
    val = [val {''}];
end

method = '?';

while ~isempty(method)
    
    method = spm_input('Inference method','0','b',lab,char(val),numel(val));
    
    switch method
        
        case 'ffx'
            
            N   = size(BMS.DCM.ffx.F,2);
            N   = 1:N;
            out = spm_api_bmc(BMS.DCM.ffx.SF,N);
            
            if ~isfield(BMS.DCM,'rfx'), method = ''; end
            
        case 'rfx'
            
            N   = size(BMS.DCM.rfx.F,2);
            N   = 1:N;
            if isfield(BMS.DCM.rfx,'model')
                out = spm_api_bmc(BMS.DCM.rfx.SF,...
                    N,BMS.DCM.rfx.model.exp_r,BMS.DCM.rfx.model.xp);
            else % Older version (prior to family level)
                out = spm_api_bmc(BMS.DCM.rfx.SF,...
                    N,BMS.DCM.rfx.exp_r,BMS.DCM.rfx.xp);
            end
            
            if ~isfield(BMS.DCM,'ffx'), method = ''; end
            
    end
    
    spm_input('Thank you',1,'d');
    
end

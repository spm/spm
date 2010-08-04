function [inverse] = spm_eeg_inv_custom_ui(D)
% GUI for parameters of inversion of forward model for EEG-MEG
% FORMAT [inverse] = spm_eeg_inv_custom_ui(D)
%
% D  - M/EEG data structure
%
% gets:
%
%     inverse.type   - 'GS' Greedy search on MSPs
%                      'ARD' ARD search on MSPs
%                      'LOR' LORETA-like model
%                      'IID' LORETA and minimum norm
%     inverse.woi    - time window of interest ([start stop] in ms)
%     inverse.Han    - switch for Hanning window
%     inverse.lpf    - band-pass filter - low frequency cut-off (Hz)
%     inverse.hpf    - band-pass filter - high frequency cut-off (Hz)
%     inverse.pQ     - any source priors (eg from fMRI) - cell array
%     inverse.xyz    - (n x 3) locations of spherical VOIs
%     inverse.rad    - radius (mm) of VOIs
%__________________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging
 
% Karl Friston
% $Id: spm_eeg_inv_custom_ui.m 4032 2010-08-04 10:15:35Z vladimir $
 
% defaults from D is specified
%==========================================================================
try
    woi = fix([D.time(1) D.time(end)]*1000);
    if (woi(end) - woi(1)) > 128
        q = 1;
    else
        q = 3;
    end
catch
    woi = [-100 200];
    q   = 1;
end
 
% get inversion parameters
%==========================================================================
inverse.type = 'GS';
if spm_input('Model','+1','b',{'Standard|Custom'},[0 1],1)
    
    % Search strategy
    %--------------------------------------------------------------------------
    type         = spm_input('Model inversion','+1','GS|COH|IID',{'GS','LOR','IID'},1);
    inverse.type = type{1};
    
    % Time window of interest
    %----------------------------------------------------------------------
    woi          = spm_input('Time-window (ms)','+1','r',woi);
    inverse.woi  = fix([min(woi) max(woi)]);
    
    % Hanning
    %----------------------------------------------------------------------
    inverse.Han  = spm_input('PST Hanning','+1','yes|no',[1 0],1);
 
    % High-pass filter
    %----------------------------------------------------------------------
    inverse.lpf  = spm_input('High-pass (Hz)','+1','0|1|8|16',[0 1 8 16],q);
    
    % Low-pass filter
    %----------------------------------------------------------------------
    inverse.hpf  = spm_input('Low-pass (Hz)','+1','48|128|256',[48 128 256],q);
    
    % Other source priors (eg from fMRI)
    %----------------------------------------------------------------------
    if spm_input('Source priors','+1','no|yes',[0 1],1)
        f = '(.*\.gii$)|(.*\.mat$)|(.*\.nii(,\d+)?$)|(.*\.img(,\d+)?$)';
        [P, sts]   = spm_select(1, f, 'Select source priors');
        if sts
            [p,f,e] = fileparts(P);
            switch lower(e)
                case '.gii'
                    g = gifti(P);
                    inverse.pQ = cell(1,size(g.cdata,2));
                    for i=1:size(g.cdata,2)
                        inverse.pQ{i} = double(g.cdata(:,i));
                    end
                case '.mat'
                    load(P);
                    inverse.pQ = pQ;
                case {'.img', '.nii'}
                    S.D = D;
                    S.fmri = P;
                    D = spm_eeg_inv_fmripriors(S);
                    inverse.fmri = D.inv{D.val}.inverse.fmri;
                    load(inverse.fmri.priors);
                    inverse.pQ = pQ;
                otherwise
                    error('Unknown file type.');
            end
        end
    else
        inverse.pQ = {};
    end
    
    % Source space restictions
    %----------------------------------------------------------------------
    if spm_input('Restrict solutions','+1','no|yes',[0 1],1)        
        [P, sts] = spm_select(1, 'mat', 'Select source (n x 3) location file');
        if sts
            xyz         = load(P);
            name        = fieldnames(xyz);
            inverse.xyz = xyz.(name{1});
            inverse.rad = spm_input('radius of VOI (mm)','+1','r',32);
        end
        
    end
end
 
return
%==========================================================================
% other GUI options
 
    % Hanning
    %----------------------------------------------------------------------
    inverse.Han = spm_input('PST Hanning','+1','yes|no',[1 0],1);
 
    % Channel modes
    %----------------------------------------------------------------------
    inverse.Nm  = spm_input('Channel modes (max)','+1','32|64|128',[32 64 128],2);
    
    % Temporal modes
    %----------------------------------------------------------------------
    inverse.Nr  = spm_input('Temporal modes (max)','+1','4|8|16',[4 8 16],2);
    
    % D.inverse.sdv    - smoothness of source priors (ms)
    %----------------------------------------------------------------------
    inverse.sdv      = spm_input('Temporal smoothness (ms)','+1','1|4|16',[1 4 16],2);
 
    % Number of sparse priors
    %----------------------------------------------------------------------
    switch inverse.type, case{'MSP','GS','ARD'}
        inverse.Np   = spm_input('MSPs per hemisphere','+1','64|128|256|512',[64 128 256 512],3);
    end
    
    % D.inverse.smooth - smoothness of source priors (mm)
    %----------------------------------------------------------------------
    switch inverse.type, case{'GS','MSP','ARD''LOR'}
        inverse.smooth = spm_input('Spatial smoothness (0-1)','+1','0.2|0.4|0.6',[0.2 0.4 0.6],3);
    end

function DCM = spm_dcm_fmri_csd_data(DCM)
% Get cross-spectral density data-features using a VAR model
% FORMAT DCM = spm_dcm_fmri_csd_data(DCM)
% DCM    -  DCM structure or fMRI
%
% sets
%
%    DCM.Y.pst     - Peristimulus Time [ms] sampled
%    DCM.Y.dt      - sampling in seconds [s] (down-sampled)
%    DCM.Y.csd     - cross spectral density over sources
%    DCM.Y.Hz      - Frequency bins
%__________________________________________________________________________
% Copyright (C) 2013 Wellcome Trust Centre for Neuroimaging
 
% Karl Friston
% $Id: spm_dcm_fmri_csd_data.m 5600 2013-08-10 20:20:49Z karl $


if ~isdeployed
    addpath(fullfile(spm('Dir'),'toolbox', 'spectral'));
end

% Time[s] of data
%--------------------------------------------------------------------------
Nb        = size(DCM.Y.y,1);              % number of bins
Nu        = size(DCM.U.u,1);              % number of bins
Nc        = size(DCM.U.u,2);              % number of inputs
DCM.Y.pst = (1:Nb)*DCM.Y.dt;              % PST

% Get frequency range
%--------------------------------------------------------------------------
try
    Hz1   = DCM.options.Fdcm(1);          % lower frequency
    Hz2   = DCM.options.Fdcm(2);          % upper frequency
catch
    Hz1   = 1/128;
    Hz2   = min(1/(2*DCM.Y.dt),1);
end

% Frequencies
%--------------------------------------------------------------------------
Nw        = 64;
DCM.Y.Hz  = linspace(Hz1,Hz2,Nw);          % Frequencies


% Cross spectral density - respones (MAR(p) model)
%==========================================================================
p         = 4;
mar       = spm_mar(DCM.Y.y,p);
mar       = spm_mar_spectra(mar,DCM.Y.Hz,1/DCM.Y.dt);
DCM.Y.csd = mar.P;

% Decimate U.u from micro-time
% -------------------------------------------------------------------------
if any(diff(DCM.U.u))
    
    Dy        = spm_dctmtx(Nb,Nb);
    Du        = spm_dctmtx(Nu,Nb);
    Dy        = Dy*sqrt(Nb/Nu);
    u         = Dy*(Du'*DCM.U.u);
    
    % Cross spectral density - inputs
    %======================================================================
    mar       = spm_mar(full(u),p);
    mar       = spm_mar_spectra(mar,DCM.Y.Hz,1/DCM.Y.dt);
    DCM.U.csd = mar.P;
    
else
    DCM.U.csd = zeros(Nw,Nc,Nc);
end

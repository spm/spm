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
%
%    DCM.U.csd     - cross spectral density of inputs
%__________________________________________________________________________
% Copyright (C) 2013 Wellcome Trust Centre for Neuroimaging
 
% Karl Friston
% $Id: spm_dcm_fmri_csd_data.m 7270 2018-03-04 13:08:10Z karl $


% check for cell array data
%--------------------------------------------------------------------------
if iscell(DCM.Y.y)
    
    % frequnecy ranges
    %----------------------------------------------------------------------
    dcm     = DCM;
    dcm.Y.y = DCM.Y.y{1};
    dcm     = spm_dcm_fmri_csd_data(dcm);
    dcm.options.Fdcm = [dcm.Y.Hz(1), dcm.Y.Hz(end)];
    
    % session-specifc spectral data features
    %----------------------------------------------------------------------
    for i = 1:numel(DCM.Y.y)
        dcm.Y.y = DCM.Y.y{i};
        dcm = spm_dcm_fmri_csd_data(dcm);
        DCM.Y.csd{i} = dcm.Y.csd;
        DCM.Y.mar{i} = dcm.Y.mar;
    end
    DCM.Y.Hz  = dcm.Y.Hz;
    DCM.Y.pst = dcm.Y.pst;
    DCM.Y.p      = dcm.Y.p;
    return
end


% add spectral toolbox
%--------------------------------------------------------------------------
if ~isdeployed
    addpath(fullfile(spm('Dir'),'toolbox', 'spectral'));
end

% Time[s] of data
%--------------------------------------------------------------------------
Nn        = size(DCM.Y.y,2);              % number of nodes
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
    Hz1   = 1/min(128,Nb*DCM.Y.dt);
    Hz2   = 1/max(8,2*DCM.Y.dt);
end

% Frequencies
%--------------------------------------------------------------------------
Nw        = 32;
DCM.Y.Hz  = linspace(Hz1,Hz2,Nw);          % Frequencies


% Cross spectral density - responses (MAR(p) model)
%==========================================================================
try
    p = DCM.Y.p;
catch
    p = 4;
end
mar       = spm_mar(DCM.Y.y,p);
mar       = spm_mar_spectra(mar,DCM.Y.Hz,1/DCM.Y.dt);
DCM.Y.csd = mar.P;
DCM.Y.p   = mar.p;

% organise MAR coefficients
%--------------------------------------------------------------------------
for i = 1:Nn
    for j = 1:Nn
        for k = 1:p
            A{i,j}(k,1) = -mar.lag(k).a(i,j);
        end
    end
end
DCM.Y.mar = spm_cat(A);


% simulated data features
%++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
% load true_parameters

% % MAR coefficients
% % -------------------------------------------------------------------------
% DCM.M.Hz  = DCM.Y.Hz;
% DCM.M.p   = DCM.Y.p;
% DCM.M.dt  = DCM.Y.dt;
% DCM.M.N   = 32;
% DCM.Y.mar = spm_csd_fmri_mar(pP,DCM.M,DCM.U);
% 
% % Cross spectral density
% % -----------------------------------------------------------------------
% DCM.M.Hz  = DCM.Y.Hz;
% DCM.M.dt  = 1;
% DCM.M.N   = 32;
% DCM.M.ns  = 1/DCM.Y.dt;
% DCM.Y.csd = spm_csd_fmri_mtf(pP,DCM.M,DCM.U);
% 
% % or generate CSD from simulated coefficients
% % -----------------------------------------------------------------------
% DCM.M.p   = 4;
% mar       = spm_csd_fmri_mar(pP,DCM.M,DCM.U);
% csd       = spm_mar2csd(mar,DCM.Y.Hz,1/DCM.M.dt);
% DCM.Y.csd = csd/sqrt(sum(abs(csd(:,1,1)).^2));

%++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++



% Cross spectral density and CCF - inputs
%==========================================================================
if any(diff(DCM.U.u))
    
    % Decimate U.u from micro-time
    % ---------------------------------------------------------------------
    Dy        = spm_dctmtx(Nb,Nb);
    Du        = spm_dctmtx(Nu,Nb);
    Dy        = Dy*sqrt(Nb/Nu);
    u         = Dy*(Du'*DCM.U.u);
     
    % Cross spectral density - inputs
    %----------------------------------------------------------------------
    mar       = spm_mar(full(u),16);
    DCM.U.csd = spm_mar2csd(mar,DCM.Y.Hz,1/DCM.Y.dt);
    
    % cross-correlation functions
    %----------------------------------------------------------------------
    Hz        = (1:64)/DCM.Y.dt/256;
    csd       = spm_mar2csd(mar,Hz,1/DCM.Y.dt);
    [ccf,pst] = spm_csd2ccf(csd,Hz,DCM.Y.dt);
    i         = 128 + (1:128);
    DCM.U.ccf = ccf(i,:,:)/max(abs(spm_vec(ccf)));
    DCM.U.pst = pst(i);

else
    
    DCM.U.csd = zeros(Nw,Nc,Nc);
    
end

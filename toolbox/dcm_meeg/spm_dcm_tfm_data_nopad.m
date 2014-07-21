function DCM = spm_dcm_tfm_data(DCM)
% gets cross-spectral density data-features using a wavelet transform
% FORMAT DCM = spm_dcm_tfm_data(DCM)
% DCM    -  DCM structure
% requires
%
%    DCM.xY.Dfile        - name of data file
%    DCM.M.U             - channel subspace
%    DCM.options.trials  - trial to evaluate
%    DCM.options.Tdcm    - time limits
%    DCM.options.Fdcm    - frequency limits
%    DCM.options.D       - Down-sampling
%
% sets
%
%    DCM.xY.pst     - Peristimulus Time [ms] sampled
%    DCM.xY.dt      - sampling in seconds [s] (down-sampled)
%    DCM.xY.U       - channel subspace
%    DCM.xY.y       - cross spectral density over channels
%    DCM.xY.csd     - cross spectral density over channels
%    DCM.xY.erp     - event-related average over channels
%    DCM.xY.It      - Indices of time bins
%    DCM.xY.Ic      - Indices of good channels
%    DCM.xY.Hz      - Frequency bins
%    DCM.xY.code    - trial codes evaluated
%    DCM.xY.scale   - scalefactor applied to data
%    DCM.xY.Rft     - Wavelet number or ratio of frequency to time
%__________________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging
 
% Karl Friston
% $Id: spm_dcm_tfm_data_nopad.m 6110 2014-07-21 09:36:13Z karl $
 
% Set defaults and Get D filename
%-------------------------------------------------------------------------

% channel indices (excluding bad channels)
%--------------------------------------------------------------------------
DCM       = spm_dcm_erp_data(DCM,0);
Ic        = DCM.xY.Ic;
It        = DCM.xY.It;
Nb        = length(It); 

% ensure spatial modes have been computed (see spm_dcm_csd)
%-------------------------------------------------------------------------
if ~isfield(DCM.M,'U')
    DCM.M.U = spm_dcm_eeg_channelmodes(DCM.M.dipfit,DCM.options.Nmodes);
end
if size(DCM.M.U,1) ~= length(Ic);
    DCM.M.U = spm_dcm_eeg_channelmodes(DCM.M.dipfit,DCM.options.Nmodes);
end
Nm          = size(DCM.M.U,2);


% options
%--------------------------------------------------------------------------
try, trial = DCM.options.trials; catch, trial = D.nconditions; end
try, Rft   = DCM.options.Rft;    catch, Rft   = 8;             end

 
% get frequency range
%--------------------------------------------------------------------------
try
    Hz1     = DCM.options.Fdcm(1);          % lower frequency
    Hz2     = DCM.options.Fdcm(2);          % upper frequency
catch
    pst     = DCM.xY.pst(end) - DCM.xY.pst(1);
    Hz1     = max(ceil(2*1000/pst),4);
    if Hz1 < 8;
        Hz2 = 48;
    else
        Hz2 = 128;
    end
end
 
 
% Frequencies
%--------------------------------------------------------------------------
Hz  = fix(Hz1:Hz2);                         % Frequencies
Nf  = length(Hz);                           % number of frequencies
Ne  = length(trial);                        % number of ERPs
 
% Cross spectral density for each trial type
%==========================================================================
erp   = cell(1,Ne);
csd   = cell(1,Ne);
for e = 1:Ne;
    
    % evoked response
    %----------------------------------------------------------------------
    Nt    = DCM.xY.nt(e);
    Y     = zeros(Nb,Nm,Nt);
    P     = zeros(Nb,Nf,Nm,Nm);
    Q     = zeros(Nb,Nf,Nm,Nm);
    for k = 1:Nt
        Y(:,:,k) = full(DCM.xY.y{e}(:,:,k)*DCM.M.U);
    end
 
    % ERP
    %----------------------------------------------------------------------
    erp{e} = mean(Y,3);
    
    % induced response
    %----------------------------------------------------------------------
    for k = 1:Nt
        
        fprintf('\nevaluating condition %i (trial %i)',e,k)
        G     = spm_morlet(Y(:,:,k) - erp{e},Hz*DCM.xY.dt,Rft);
        for i = 1:Nm
            for j = 1:Nm
                P(:,:,i,j) = (G(:,:,i).*conj(G(:,:,j)));
            end
        end
        Q = Q + P;
    end
    
    % normalise induced responses in relation to evoked responses
    %--------------------------------------------------------------------------
    Vm    = mean(mean(squeeze(var(Y,[],3))));
    Vs    = mean(diag(squeeze(mean(squeeze(mean(Q))))));
    Q     = Vm*Q/Vs;
    
    % store
    %----------------------------------------------------------------------
    csd{e} = Q;
    
end
 
% place cross-spectral density in xY.y
%==========================================================================
try
    scale = DCM.xY.scale;
catch
    scale = 1;
end

DCM.xY.erp   = spm_unvec(spm_vec(erp)*(scale^1),erp);
DCM.xY.csd   = spm_unvec(spm_vec(csd)*(scale^2),csd);
DCM.xY.y     = DCM.xY.csd;
DCM.xY.U     = DCM.M.U;
DCM.xY.scale = scale;
DCM.xY.Hz    = Hz;
DCM.xY.Rft   = Rft;

 
return
 
% plot responses
%--------------------------------------------------------------------------
spm_dcm_tfm_response(DCM.xY,DCM.xY.pst,DCM.xY.Hz);

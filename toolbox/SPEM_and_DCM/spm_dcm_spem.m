function DCM = spm_dcm_spem(DCM)
% Estimate parameters of a DCM of smooth pursuit eye movements
% FORMAT DCM = spm_dcm_spem(DCM)
%
% DCM
%    DCM.name: name string
%    DCM.xY: data   {1 x nc struct}
%
%        xY.Y{i}  - eye-gaze position for i-th condition
%        xY.C{i}  - target position for i-th condition
%        xY.DT    - time bin (ms)
%        xY.occ   - occlusion function occ(x) = {0,1}: -1 > x > 1
%        xY.C{i}  - target position for i-th condition
%
%    DCM.xU: design [nu x nc array]
%    DCM.pE: prior expectation
%    DCM.pC: prior covariance
%
% This routine checks the data and inverts a meta-model of observed slow
% pursuit eye movements using the standard variational Laplacian scheme
%
% See also: spm_SEM_gen; spm_dcm_spem_data; spm_dcm_spem_results
%__________________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% Karl Friston
% $Id: spm_dcm_spem.m 6014 2014-05-23 15:00:35Z guillaume $

% name
%--------------------------------------------------------------------------
if ~isfield(DCM,'name')
    DCM.name =  ['DCM-' date '.mat'];
end

% check data
%--------------------------------------------------------------------------
DCM.xY = spm_dcm_spem_data(DCM.xY);

% check occluder function
%--------------------------------------------------------------------------
if isfield(DCM.xY,'occ')
    M.occ = DCM.xY.occ;
else
    M.occ = @(x) x - x + 1;
end


% model specification
%--------------------------------------------------------------------------
M.IS    = 'spm_SEM_gen';
M.pE    = DCM.pE;
M.pC    = DCM.pC;
M.hE    = (8 + 2);
M.hC    = exp(-8);
M.ns    = length(DCM.xY.y{1});
M.u     = DCM.xY.u;
M.w     = 2*pi/M.ns;
M.x     = DCM.xY.x;
M.Nmax  = 64;


% Variational Laplace: model inversion: invert active inference model
%==========================================================================
[Ep,Cp,Eh,F] = spm_nlsi_GN(M,DCM.xU,DCM.xY);

% posterior predictions
%--------------------------------------------------------------------------
[Y,DEM] = spm_SEM_gen(Ep,M,DCM.xU);

% report solution for the first condition
%--------------------------------------------------------------------------
spm_figure('GetWin','DEM');
spm_DEM_qU(DEM{1}.qU,DEM{1}.pU)


% store estimates in DCM
%--------------------------------------------------------------------------
DCM.M   = M;                    % (Meta) model
DCM.Ep  = Ep;                   % conditional expectation
DCM.Cp  = Cp;                   % conditional covariance
DCM.Ce  = exp(-Eh);             % ReML error covariance
DCM.F   = F;                    % Laplace log evidence
DCM.Y   = Y;                    % predicted responses
DCM.DEM = DEM;                  % inference model

% save
%--------------------------------------------------------------------------
save(DCM.name, 'DCM', spm_get_defaults('mat.format'));

% and report predictions over conditions
%--------------------------------------------------------------------------
spm_figure('GetWin','DCM for SPEM');
spm_dcm_spem_results(DCM)

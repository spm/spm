function DCM = spm_dcm_erp_dipfit(DCM, save_vol_sens)
% Prepare structures for ECD forward model (EEG, MEG and LFP)
% FORMAT DCM = spm_dcm_erp_dipfit(DCM, save_vol_sens)
% DCM           - DCM structure
% save_vol_sens - optional argument indicating whether to perform
%                 the time consuming step required for actually using
%                 the forward model to compute lead fields (1, default)
%                 or skip it if the function is only called for
%                 verification of the input (0).
%
% Input DCM structure requires:
%       DCM.xY.Dfile
%       DCM.xY.Ic
%       DCM.Lpos
%       DCM.options.spatial - 'ERP', 'LFP' or 'IMG'
%
% fills in:
%
%       DCM.M.dipfit
%
%    dipfit.location - 0 or 1 for source location priors
%    dipfit.symmetry - 0 or 1 for symmetry constraints on sources
%    dipfit.modality - 'EEG', 'MEG', 'MEGPLANAR' or 'LFP'
%    dipfit.type     - 'ECD', 'LFP' or 'IMG''
%    dipfit.symm     - distance (mm) for symmetry constraints (ECD)
%    dipfit.Lpos     - x,y,z source positions (mm)            (ECD)
%    dipfit.Nm       - number of modes                        (Imaging)
%    dipfit.Ns       - number of sources
%    dipfit.Nc       - number of channels
%
%    dipfit.vol      - volume structure (for M/EEG)
%    dipfit.datareg  - registration structure (for M/EEG)
%__________________________________________________________________________
% Copyright (C) 2007-2015 Wellcome Trust Centre for Neuroimaging
 
% Karl Friston
% $Id: spm_dcm_erp_dipfit.m 6360 2015-03-04 19:24:56Z spm $
 
% Get data filename and good channels
%--------------------------------------------------------------------------
try
    DCM.xY.Dfile;
    DCM.M.dipfit.Ic = DCM.xY.Ic;
catch
    spm('alert*','Please specify data');
    error('Data not specified.')
end
 
if nargin == 1
  save_vol_sens = 0;
end
    

% D - SPM data structure
%--------------------------------------------------------------------------
try
    D  = spm_eeg_load(DCM.xY.Dfile);
catch
    D  = spm_eeg_load(spm_file(DCM.xY.Dfile,'filename'));
end
 
 
% set options in dipfit
%--------------------------------------------------------------------------
try, spatial  = DCM.options.spatial;  catch, spatial  = 'IMG'; end
try, location = DCM.options.location; catch, location = 0;     end
try, symmetry = DCM.options.symmetry; catch, symmetry = 0;     end

DCM.M.dipfit.type     = spatial;
DCM.M.dipfit.location = location;
DCM.M.dipfit.symmetry = symmetry;


% Get source locations if MEG or EEG
%--------------------------------------------------------------------------
switch DCM.xY.modality
 
    % get source priors for EEG or MEG
    %----------------------------------------------------------------------
    case{'EEG','MEG','MEGPLANAR'}
 
        [D, ok] = check(D, '3d');
        if ~ok
            spm('alert*',{'File not ready for source reconstruction.',...
                'Use Prepare to specify sensors and fiducials.'});
        end
 
        try
            DCM.M.dipfit.Lpos = DCM.Lpos;
        catch
            spm('alert*',{'Please specify source locations','in DCM.Lpos'})
        end
 
        DCM.M.dipfit.modality = DCM.xY.modality;
        DCM.M.dipfit.Ns       = length(DCM.Sname);
        DCM.M.dipfit.Nc       = length(DCM.xY.Ic);
 
        % otherwise assume LFP
        %------------------------------------------------------------------
    otherwise
 
        DCM.M.dipfit.modality = 'LFP';
        DCM.M.dipfit.Ns       = length(DCM.Sname);
        DCM.M.dipfit.Nc       = length(DCM.xY.Ic);
        return
end

% Detect silent sources for non-LFP type from GUI source name text
%--------------------------------------------------------------------------
DCM.M.dipfit.silent_source  = strfind(DCM.Sname,'silent');

% If not LFP, get electromagnetic forward model
%==========================================================================
if ~isfield(D, 'val'), D.val = 1; end

% detects old version of the struct 
%--------------------------------------------------------------------------
if ~isfield(D, 'inv') || ~iscell(D.inv) ||...
        ~(isfield(D.inv{D.val}, 'forward') && isfield(D.inv{D.val}, 'datareg')) ||...
        ~isa(D.inv{D.val}.mesh.tess_ctx, 'char') 
    D = spm_eeg_inv_mesh_ui(D, D.val);
    D = spm_eeg_inv_datareg_ui(D, D.val);
    D = spm_eeg_inv_forward_ui(D, D.val);
    save(D);
end

% fill in dipfit
%--------------------------------------------------------------------------
for m = 1:numel(D.inv{D.val}.forward)
    if strncmp(DCM.xY.modality, D.inv{D.val}.forward(m).modality, 3)
        DCM.M.dipfit.vol      = D.inv{D.val}.forward(m).vol;
        DCM.M.dipfit.datareg  = D.inv{D.val}.datareg(m);
        if isfield(D.inv{D.val}.forward(m), 'siunits') && D.inv{D.val}.forward(m).siunits
            DCM.M.dipfit.sens     = D.inv{D.val}.forward(m).sensors;
            DCM.M.dipfit.siunits  = true;
        else            
            DCM.M.dipfit.sens     = DCM.M.dipfit.datareg.sensors;
            DCM.M.dipfit.siunits  = false;
        end
    end
end


% channels
%--------------------------------------------------------------------------
if save_vol_sens
    if ischar(DCM.M.dipfit.vol)
        DCM.M.dipfit.vol = ft_read_vol(DCM.M.dipfit.vol);
    end

    [DCM.M.dipfit.vol, DCM.M.dipfit.sens] = ft_prepare_vol_sens(DCM.M.dipfit.vol, ...
        DCM.M.dipfit.sens, 'channel', D.chanlabels(DCM.xY.Ic));
end

switch DCM.options.spatial
 
    % Imaging (distributed source reconstruction)
    %----------------------------------------------------------------------
    case{'IMG'}
 
        % Load Gain or Lead field matrix
        %------------------------------------------------------------------
        DCM.val = D.val;
        [L,D]   = spm_eeg_lgainmat(D, [], D.chanlabels(DCM.xY.Ic));
        
        % centers
        %------------------------------------------------------------------
        xyz = DCM.M.dipfit.Lpos;
        Np  = size(xyz,2);
 
        % parameters
        %==================================================================
 
        % defaults: Nm = 6; number of modes per region
        %------------------------------------------------------------------
        try, rad  = DCM.M.dipfit.radius; catch, rad  = 16;    end
        try, Nm   = DCM.M.dipfit.Nm;     catch, Nm   = 6;     end
        
 
        % Compute spatial basis (eigenmodes of lead field)
        %==================================================================
 
        % create MSP spatial basis set in source space
        %------------------------------------------------------------------
        vert   = D.inv{D.val}.mesh.tess_mni.vert;
        for i  = 1:Np
            Dp = sum([vert(:,1) - xyz(1,i), ...
                      vert(:,2) - xyz(2,i), ...
                      vert(:,3) - xyz(3,i)].^2,2);
 
            % nearest mesh points
            %--------------------------------------------------------------
            Ip = find(Dp < rad^2);
            if length(Ip) < Nm;
                [y,Ip] = sort(Dp);
                Ip     = Ip(1:Nm);
            end
 
            % left hemisphere
            %--------------------------------------------------------------
            U                  = spm_svd(L(:,Ip)',0);
            U                  = U(:,1:Nm);
            DCM.M.dipfit.G{i}  = L(:,Ip)*U;
            DCM.M.dipfit.U{i}  = U;
            DCM.M.dipfit.Ip{i} = Ip;
        end
 
        % Save results
        %==================================================================
        DCM.M.dipfit.radius  = rad;                   % VOI (XYZ, Radius)
        DCM.M.dipfit.Nm      = Nm;                    % modes per region
        DCM.M.dipfit.Nd      = length(vert);          % number of dipoles
        DCM.M.dipfit.gainmat = D.inv{D.val}.gainmat;  % Lead field filename
 
    otherwise
 
end

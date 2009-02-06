function DCM = spm_dcm_erp_dipfit(DCM, save_vol_sens)
% prepares structures for ECD forward model (EEG, MEG and LFP)
% FORMAT DCM = spm_dcm_erp_dipfit(DCM, save_vol_sens)
% save_vol_sens - optional argument indicating whether to perform
%                 the time consuming step required for actually using
%                 the forward model to compute lead fields (1)
%                 or skip it if the function is only called for
%                 verification of the input (0, default).
% requires:
%
% needs:
%       DCM.xY.Dfile
%       DCM.xY.Ic
%       DCM.Lpos
%       DCM.options.spatial     - 'ERP', 'LFP' or 'IMG'
%       DCM.M.dipfit.sensorfile - (for ECD models)
%
% fills in:
%
%       DCM.M.dipfit
%
%    dipfit.location - 0 or 1 for source location priors
%    dipfit.symmetry - 0 or 1 for symmetry constraints on sources
%    dipfit.modality - 'EEG', 'MEG' or 'LFP'
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
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% Karl Friston
% $Id: spm_dcm_erp_dipfit.m 2709 2009-02-06 19:56:19Z karl $

% Get data filename and good channels
%--------------------------------------------------------------------------
try
    DCM.xY.Dfile;
    DCM.M.dipfit.Ic = DCM.xY.Ic;
catch
    errordlg('Please specify data');
    error('')
end

if nargin == 1
    save_vol_sens = 0;
end


% D - SPM data structure
%--------------------------------------------------------------------------
D    = spm_eeg_load(DCM.xY.Dfile);


% set options in dipfit
%--------------------------------------------------------------------------
try, spatial  = DCM.options.spatial;  catch, spatial  = 'IMG'; end
try, model    = DCM.options.model;    catch, model    = 'NMM'; end
try, location = DCM.options.location; catch, location = 0;     end
try, symmetry = DCM.options.symmetry; catch, symmetry = 0;     end

DCM.M.dipfit.type     = spatial;
DCM.M.dipfit.model    = model;
DCM.M.dipfit.location = location;
DCM.M.dipfit.symmetry = symmetry;


% Get source locations if MEG or EEG
%--------------------------------------------------------------------------
switch DCM.xY.modality

    % get source priors for EEG or MEG
    %----------------------------------------------------------------------
    case{'EEG','MEG'}

        [ok, D] = check(D, 'sensfid');
        if ~ok
            if check(D, 'basic')
                errordlg(['File not ready for source reconstruction.'...
                    'Use prep to specify sensors and fiducials.']);
            else
                errordlg('The meeg file is corrupt or incomplete');
            end
        end

        try
            DCM.M.dipfit.Lpos = DCM.Lpos;
        catch
            errordlg({'Please specify source locations','in DCM.Lpos'})
        end

        DCM.M.dipfit.modality  = DCM.xY.modality;
        DCM.M.dipfit.Ns        = length(DCM.Sname);
        DCM.M.dipfit.Nc        = length(DCM.xY.Ic);

        % otherwise assume LFP
        %------------------------------------------------------------------
    otherwise

        DCM.M.dipfit.modality  = 'LFP';
        DCM.M.dipfit.Ns        = length(DCM.Sname);
        DCM.M.dipfit.Nc        = length(DCM.xY.Ic);
        return
end

% If not LFP, get electromagnetic forward model
%==========================================================================
if ~isfield(D, 'val'), D.val = 1; end

try
    DCM.M.dipfit.vol     = D.inv{D.val}.forward.vol;
    DCM.M.dipfit.datareg = D.inv{D.val}.datareg;
catch
    D = spm_eeg_inv_mesh_ui(D, D.val, [], 4);
    D = spm_eeg_inv_datareg_ui(D, D.val);
    DCM.M.dipfit.datareg = D.inv{D.val}.datareg;
    save(D);
end

% channels
%--------------------------------------------------------------------------
chanind = strmatch(DCM.xY.modality, D.chantype, 'exact');
chanind = setdiff(chanind, D.badchannels);
if ~isempty(chanind)
    channels = D.chanlabels(chanind);
else
    error(['No good ' DCM.xY.modality ' channels were found.']);
end

if save_vol_sens
    [DCM.M.dipfit.vol, DCM.M.dipfit.sens] = forwinv_prepare_vol_sens(D.inv{D.val}.forward.vol, ...
        D.inv{D.val}.datareg.sensors, 'channel', channels);
end

switch DCM.options.spatial

    % Imaging (distributed source reconstruction)
    %----------------------------------------------------------------------
    case{'IMG'}

        % Load Gain or Lead field matrix
        %------------------------------------------------------------------
        DCM.val = D.val;
        [L D]   = spm_eeg_lgainmat(D);

        % centers
        %------------------------------------------------------------------
        xyz = DCM.M.dipfit.Lpos;
        Np  = size(xyz,2);

        % parameters
        %==================================================================

        % defaults: Nm = 8; number of modes per region
        %------------------------------------------------------------------
        try, rad  = DCM.M.dipfit.radius; catch, rad  = 16;    end
        try, Nm   = DCM.M.dipfit.Nm;     catch, Nm   = 8;     end

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
        DCM.M.dipfit.radius  = rad;                           % VOI (XYZ, Radius)
        DCM.M.dipfit.Nm      = Nm;                            % modes per region
        DCM.M.dipfit.Nd      = length(vert);                  % number of dipoles
        DCM.M.dipfit.gainmat = D.inv{D.val}.forward.gainmat;  % Lead field filename

    otherwise

end

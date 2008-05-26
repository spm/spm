function DCM = spm_dcm_erp_dipfit(DCM)
% prepares structures for ECD forward model (both EEG and MEG)
% FORMAT DCM = spm_dcm_erp_dipfit(DCM)
% requires:
%
% needs:
%       DCM.xY.Dfile
%       DCM.xY.Ic
%       DCM.Lpos
%       DCM.options.type        - 'ECD' (1) or 'imaging' (2)
%       DCM.M.dipfit.sensorfile - (for ECD models)
% fills in:
%
%       DCM.M.dipfit
%
%    dipfit.modality - 'EEG', 'MEG'     or 'LFP'
%    dipfit.type     - 'ECD', 'Imaging' or 'LFP'
%    dipfit.symm     - distance (mm) for symmetry constraints (ECD)
%    dipfit.Lpos     - x,y,z source positions (mm)            (ECD)
%    dipfit.Nm       - number of modes                        (Imaging)
%    dipfit.Ns       - number of sources
%__________________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% Karl Friston
% $Id: spm_dcm_erp_dipfit.m 1726 2008-05-26 16:45:55Z vladimir $

% Get data filename and good channels
%--------------------------------------------------------------------------
try
    DCM.xY.Dfile;
    DCM.M.dipfit.Ic = DCM.xY.Ic;
catch
    errordlg('Please specify data');
    error('')
end

% D - SPM data structure
%--------------------------------------------------------------------------
D    = spm_eeg_load(DCM.xY.Dfile);

% Get source locations if MEG or EEG
%--------------------------------------------------------------------------
switch DCM.xY.modality

    % get source priors for EEG or MEG
    %----------------------------------------------------------------------
    case{'EEG','MEG'}

        [ok, D] = check(D, 'sensfid');

        if ~ok
            if check(D, 'basic')
                errordlg(['The requested file is not ready for source reconstruction.'...
                    'Use prep to specify sensors and fiducials.']);
            else
                errordlg('The meeg file is corrupt or incomplete');
            end
        end

        try
            DCM.M.dipfit.Lpos = DCM.Lpos;
        catch
            try
                DCM.Lpos = DCM.M.dipfit.Lpos;
            catch
                errordlg({'Please specify source locations','in DCM.Lpos'})
            end
        end

        DCM.M.dipfit.modality  = DCM.xY.modality;
        DCM.M.dipfit.Ns        = size(DCM.Lpos,2);

        % othrerwise assume LFP
        %----------------------------------------------------------------------
    otherwise
        DCM.M.dipfit.modality  = 'LFP';
        DCM.M.dipfit.type      = 'LFP';
        DCM.M.dipfit.Ns        = size(DCM.C,1);
        return
end

if ~isfield(D, 'val')
    D.val = 1;
end

try
    DCM.M.dipfit.vol = D.inv{D.val}.forward.vol;
    DCM.M.dipfit.datareg = D.inv{D.val}.datareg;
catch
    D = spm_eeg_inv_mesh_ui(D, D.val, [], 4);
    D = spm_eeg_inv_datareg_ui(D, D.val);
    DCM.M.dipfit.datareg = D.inv{D.val}.datareg;
end

[DCM.M.dipfit.vol, DCM.M.dipfit.sens] = forwinv_prepare_vol_sens(D.inv{D.val}.forward.vol, ...
    D.inv{D.val}.datareg.sensors, 'channel', D.inv{D.val}.forward.channels);

switch DCM.options.type
    case 1 %ECD
        DCM.M.dipfit.type = 'ECD';
    case 2  % Imaging (distributed source reconstruction)
        DCM.M.dipfit.type = 'Imaging';

        % Load Gain or Lead field matrix
        %------------------------------------------------------------------
        DCM.val = D.val;
        try
            gainmat = D.inv{D.val}.forward.gainmat;
        catch
            D = spm_eeg_inv_forward_ui(D, D.val);
            gainmat = D.inv{D.val}.forward.gainmat;
        end

        try
            L       = load(gainmat);
        catch
            [p,f]   = fileparts(gainmat);
            L       = load(f);
            gainmat = fullfile(pwd,f);
        end

        name    = fieldnames(L);
        L       = sparse(getfield(L, name{1}));
        L       = spm_cond_units(L);

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
        DCM.M.dipfit.radius  = rad;                  % VOI (XYZ, Radius)
        DCM.M.dipfit.Nm      = Nm;                   % modes per region
        DCM.M.dipfit.Nd      = length(vert);         % number of dipoles
        DCM.M.dipfit.gainmat = gainmat;              % Lead field filename
        DCM.M.dipfit.type    = 'Imaging';

        save(D);
end
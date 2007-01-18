function DCM = spm_dcm_erp_dipfit(DCM,sensorfile,Dfile)
% prepares structures for ECD forward model (both EEG and MEG)
% FORMAT DCM = spm_dcm_erp_dipfit(DCM,Dfile,sensorfile)
% Stefan Kiebel
% $Id: spm_dcm_erp_prepareSpatial.m 668 2006-10-26 16:35:28Z karl $

% Put data file in DCM.xY
%--------------------------------------------------------------------------
try
    DCM.xY.Dfile = Dfile;
catch
    try
        Dfile = DCM.xY.Dfile;
    catch
        try
            Dfile = DCM.M.Dfile;
            DCM.xY.Dfile = Dfile;
        catch
            [f p] = uigetfile({'*.mat'},'Please select data file');
            DCM.xY.Dfile = fullfile(p,f);
            Dfile = DCM.xY.Dfile;
        end
    end
end

% Polhemus file skip
%==========================================================================
try
    Pol_skip = DCM.options.Pol_skip;
catch
    Pol_skip = 4;
    DCM.options.Pol_skip = Pol_skip;
end
%==========================================================================

% EEG
%--------------------------------------------------------------------------
if DCM.options.Spatial_type == 1

    % Get Polyhmus file
    %----------------------------------------------------------------------
    try
        DCM.M.dipfit.sensorfile = sensorfile;
    catch
        try
            sensorfile = DCM.M.dipfit.sensorfile;
        catch
            try
                sensorfile = DCM.xY.sensorfile;
            catch
                [f p] = uigetfile({'*.mat'},'Please select sensor location file');
                DCM.M.dipfit.sensorfile = fullfile(p,f);
                sensorfile = DCM.xY.dipfit.sensorfile;
            end
        end
    end

    % Polyhmus file
    %----------------------------------------------------------------------
    if strcmpi('pol', spm_str_manip(sensorfile, 'e'))
        try
            % returns coordinates of fiducials and sensors (filenames)
            % coordinates are in CTF coordinate system
            [fid, sens] = spm_eeg_inv_ReadPolhemus(sensorfile, Pol_skip);
        catch
            try
                DCM.M.dipfit.elc;
                return
            catch
                warndlg('Could not read sensor file')
                return
            end
        end
    else
        try

            % mat-file must contain 2 matrices with fid and sensor
            % coordinates
            tmp = load(sensorfile, '-mat');
            fid = tmp.fid;
            sens = tmp.sens;
            
        catch
            try
                DCM.M.dipfit.elc;
                return
            catch
                warndlg('Could not read sensor file')
                return
            end
        end
    end

    % Use approximated MRI fiducials (MNI space)... for now
    %----------------------------------------------------------------------
    mni_fid = [[0   86 -42];...
              [-84 -18 -55];...
              [ 84 -18 -55]];
          
          
    % EEG coordinates in MNI space
    %----------------------------------------------------------------------
    % coregistration between EEG and MRI fiducials
    [eeg2mri,sensreg,fid_reg] = spm_eeg_inv_datareg(sens, fid, mni_fid);

    try
        chansel = DCM.M.dipfit.chansel;
    catch
        D       = load(DCM.xY.Dfile);
        chansel = setdiff(D.D.channels.eeg, D.D.channels.Bad);
        DCM.M.dipfit.chansel = chansel;
    end
    dipfit  = DCM.M.dipfit;
    elc     = sensreg(chansel, :);

    % fit origin of outer sphere, with fixed radius, to electrodes
    %----------------------------------------------------------------------
    dipfit.vol.r = [71 72 79 85];
    dipfit.vol.c = [0.3300 1 0.0042 0.3300];
    dipfit.vol.o = [0 0 0]; % origin is always zero (forward model can't handle non-zero origin)

    % origin of best fitting sphere (in MNI space)
    %----------------------------------------------------------------------
    dipfit.vol.o_sphere = fminsearch(@(x) d_elc_sphere(x, sensreg, dipfit.vol.r(4)),[0 0 0]);

    % compute POL->MNI coordinate transformation, following definition of
    % CTF origin, and axes
    %----------------------------------------------------------------------

    % origin in MNI space
    %----------------------------------------------------------------------
    POLorigin = mean(fid_reg(2:3, :));
    POLvx = (fid_reg(1, :) - POLorigin)/norm((fid_reg(1, :) - POLorigin));
    POLvz = cross(POLvx, fid_reg(2, :) - fid_reg(3, :));
    POLvz = POLvz/norm(POLvz);  % Vector normal to the NLR plane, pointing upwards
    POLvy = cross(POLvz,POLvx); POLvy = POLvy/norm(POLvy);

    % transformation POL->MNI. Must be the same as computed by registration
    % (file RT).
    %----------------------------------------------------------------------
    dipfit.Mpol2mni = [[POLvx' POLvy' POLvz' POLorigin']; [0 0 0 1]];
    dipfit.MNIelc   = sensreg(chansel, :); % save

    % original coordinates
    elc = inv(dipfit.Mpol2mni)*[dipfit.MNIelc'; ones(1, size(dipfit.MNIelc, 1))];
    elc = elc(1:3,:)';

    % centre of sphere in POL space
    o   = inv(dipfit.Mpol2mni)*[dipfit.vol.o_sphere'; 1];
    elc = elc - kron(ones(size(elc, 1), 1), o(1:3)');

    % transformation matrix from MNI-space to sphere in POL space
    Msphere = eye(4);
    Msphere(1:3,4) = -o(1:3);
    dipfit.Mpol2sphere = Msphere; % in pol-space: translation of origin to zero
    dipfit.Mmni2polsphere = Msphere*inv(dipfit.Mpol2mni); % from MNI-space to pol space (plus translation to origin)

    % projecting channels to outer sphere (remove origin, necessary for forward model)
    dist = sqrt(sum((elc).^2,2));
    dipfit.elc  = dipfit.vol.r(4) * elc ./[dist dist dist];

    % save in M
    %----------------------------------------------------------------------
    DCM.M.dipfit = dipfit;


%MEG
%--------------------------------------------------------------------------
elseif DCM.options.Spatial_type == 2

    % not much to do, because the sensors location/orientations were
    % already read at the time of conversion.

    DCM.M.dipfit.vol.r = [85];
    DCM.M.dipfit.vol.o = [0 0 0];
end

% used for fitting sphere to sensors
%--------------------------------------------------------------------------
function x = d_elc_sphere(x, elc, r)

% returns sum of squares of distances between sensors elc and sphere.
Nelc = size(elc, 1);
d    = sqrt(sum((elc - repmat(x, Nelc, 1)).^2, 2));
x    = sum((d - r).^2);



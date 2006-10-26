function DCM = spm_dcm_erp_prepareSpatial(DCM,Dfile,sensorfile)
% prepares structures for ECD forward model (both EEG and MEG)

% Stefan Kiebel
% $Id: spm_dcm_erp_prepareSpatial.m 668 2006-10-26 16:35:28Z karl $

% Put data file in DCM.Y
%--------------------------------------------------------------------------
try
    DCM.Y.Dfile = Dfile;
catch
    try
        Dfile = DCM.Y.Dfile;
    catch
        [f p] = uigetfile({'*.mat'},'Please data file');
        DCM.Y.Dfile = fullfile(p,f);
        Dfile = DCM.Y.Dfile;
    end
end

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
            [f p] = uigetfile({'*.mat'},'Please select sensor location file');
            DCM.M.dipfit.sensorfile = fullfile(p,f);
            sensorfile = DCM.M.dipfit.sensorfile;
        end
    end
    
    % Polyhmus file
    %----------------------------------------------------------------------
    if strcmpi('pol', spm_str_manip(sensorfile, 'e'))
        try
            % returns coordinates of fiducials and sensors (filenames)
            % coordinates are in CTF coordinate system
            [Fname_fid, Fname_sens] = spm_eeg_inv_ReadPolhemus(sensorfile, 0);
        catch
            try
                DCM.M.dipfit.vol;
                return
            catch
                warndlg('Could not read sensor file')
                return
            end
        end
    else
        try
            % must be in CTF-space
            tmp = load(sensorfile, '-mat');
            xyz = fieldnames(tmp);
            eval(['xyz = tmp.' xyz{1}]);
            x   = xyz(:, 1);
            y   = xyz(:, 2);
            z   = xyz(:, 3);
        catch
            try
                DCM.M.dipfit.vol;
                return
            catch
                warndlg('Could not read sensor file')
                return
            end
        end
    end

% Use approximated MRI fiducials (MNI space)... for now
%----------------------------------------------------------------------
MNI_FID = [[0.0 86 -39];...
          [-85 -17 -35];...
          [ 84 -14 -35]];

save MNI_FID MNI_FID
mni_fid = 'MNI_FID';

% coregistration between EEG and MRI fiducials, and transforms EEG
% coordinates to MNI space
%----------------------------------------------------------------------
[RT, sensors_reg, fid_reg] = spm_eeg_inv_datareg(1,Fname_sens, Fname_fid, mni_fid);

% EEG coordinates in MNI space
%----------------------------------------------------------------------
load(sensors_reg);

% hard-coded for now (Olivier measured as first 2 channels ear-lobes ->
% non-standard), and first two channels are (always) CMS and DRL)
%----------------------------------------------------------------------
try
    chansel = DCM.M.dipfit.chansel;
catch
    try
        D = load(DCM.Y.Dfile);
    catch
        D = load(DCM.M.Dfile);
    end
    chansel = setdiff(D.D.channels.eeg, D.D.channels.Bad);
    DCM.M.dipfit.chansel = chansel;
end
dipfit  = DCM.M.dipfit;
sensreg = sensreg(5:end, :);
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
load(fid_reg);

% origin in MNI space
%----------------------------------------------------------------------
POLorigin = mean(fideegreg(2:3, :));
POLvx = (fideegreg(1, :) - POLorigin)/norm((fideegreg(1, :) - POLorigin));
POLvz = cross(POLvx, fideegreg(2, :) - fideegreg(3, :));
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
elc  = dipfit.vol.r(4) * elc ./[dist dist dist];

dipfit.elc   = elc;
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



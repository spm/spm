<<<<<<< .mine
function DCM = spm_dcm_erp_dipfit(DCM)
% prepares structures for ECD forward model (both EEG and MEG)
% FORMAT DCM = spm_dcm_erp_dipfit(DCM)
% requires:
% 
% needs:
%       DCM.xY.Dfile
%       DCM.xY.Ic
%       DCM.Lpos
%       DCM.options.type
%       DCM.M.dipfit.sensorfile - 'ECD'
%
% fills in:
%       
%       DCM.M.dipfit.type       - 'ECD (EEG)','ECD (MEG)','Imaging','LFP'
%       DCM.M.dipfit. ...
%__________________________________________________________________________
% Copyright (C) 2005 Wellcome Trust Centre for Neuroimaging
 
% Karl Friston
% $Id: spm_dcm_erp_dipfit.m 1040 2007-12-21 20:28:30Z karl $

% Get data filename and good channels
%--------------------------------------------------------------------------
try
    DCM.xY.Dfile;
    DCM.M.dipfit.Ic = DCM.xY.Ic;
catch
    errordlg('Please specify data');
    error('')
end

% Get source locations
%--------------------------------------------------------------------------
try
    DCM.M.dipfit.L.pos = DCM.Lpos;
catch
    try
        DCM.Lpos = DCM.M.dipfit.L.pos;
    catch
        errordlg({'Please specify source locations',...
                  'in DCM.Lpos'})
    end
end


% Polhemus file skip
%==========================================================================
try
    Pol_skip = DCM.options.Pol_skip;
catch
    Pol_skip = 2;
    DCM.options.Pol_skip = Pol_skip;
end
%==========================================================================


% D - SPM data structure
%==========================================================================
D       = spm_eeg_ldata(DCM.xY.Dfile);

switch DCM.options.type

    % EEG - ECD
    %----------------------------------------------------------------------
    case 1
    if strcmp(D.modality,'EEG')

        % try to get sen_reg and fid_reg from D.inv
        %------------------------------------------------------------------
        try
            sen_reg = D.inv{1}.datareg.sens_coreg;
            fid_reg = D.inv{1}.datareg.fid_coreg;
        catch

            % Get sensor file
            %--------------------------------------------------------------
            try
                sensorfile = DCM.M.dipfit.sensorfile;
            catch
                [f,p]      = uigetfile({'*.mat;*.pol'},'select sensor file');
                sensorfile = fullfile(p,f);
                DCM.M.dipfit.sensorfile = sensorfile;
            end

            % Polyhmus file
            %--------------------------------------------------------------
            if strcmpi('pol', spm_str_manip(sensorfile, 'e'))
                
                % returns coordinates of fiducials and sensors (filenames)
                % coordinates are in CTF coordinate system
                %----------------------------------------------------------
                [fid, sens] = spm_eeg_inv_ReadPolhemus(sensorfile, Pol_skip);
                
            else
 
                % mat-file must contain 2 matrices with fid and sensor coordinates
                %----------------------------------------------------------
                tmp  = load(sensorfile, '-mat');
                fid  = tmp.fid;
                sens = tmp.sens;
            end

            % Use approximated MRI fiducials (MNI space)... for now
            %--------------------------------------------------------------
            mni_fid = [[0   86 -42];...
                      [-84 -18 -55];...
                      [ 84 -18 -55]];


            % EEG coordinates in MNI space
            % coregistration between EEG and MRI fiducials
            %--------------------------------------------------------------
            [eeg2mri,sen_reg,fid_reg] = spm_eeg_inv_datareg(sens, fid, mni_fid);

            try
                Ic = DCM.M.dipfit.Ic;
            catch
                D  = spm_eeg_ldata(DCM.xY.Dfile);
                Ic = setdiff(D.channels.eeg, D.channels.Bad);
                DCM.M.dipfit.Ic = Ic;
            end
            sen_reg  = sen_reg(Ic, :);

        end

        % evaluate dipfit and save in M
        %------------------------------------------------------------------
        DCM.M.dipfit       = spm_dipfit(DCM.M.dipfit,sen_reg,fid_reg);
        DCM.M.dipfit.type  = 'ECD (EEG)';
    


    % MEG - ECD
    %----------------------------------------------------------------------
    else

        % not much to do because the sensors location/orientations were
        % already read at the time of conversion.
        %------------------------------------------------------------------
        DCM.M.dipfit.vol.r = [85];
        DCM.M.dipfit.vol.o = [0 0 0];
        DCM.M.grad         = D.channels.grad;
        DCM.M.dipfit.type  = 'ECD (MEG)';
        
    end
        
    % Imaging (distributed source reconstruction)
    %----------------------------------------------------------------------
    case 2

        % Load Gain or Lead field matrix
        %------------------------------------------------------------------
        try
            DCM.val = D.val;
            gainmat = D.inv{D.val}.forward.gainmat;
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
        catch
            errordlg({'Please create and save a foward model',...
                      'using spm_eeg_inv_imag_api'})
            error('')
        end

        % centers
        %------------------------------------------------------------------
        xyz = DCM.M.dipfit.L.pos;
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
end

return


function dipfit = spm_dipfit(dipfit,sen_reg,fid_reg)
%==========================================================================
% computes dipfit model for ECD, DCM

% fit origin of outer sphere, with fixed radius, to electrodes
%--------------------------------------------------------------------------
dipfit.vol.r = [71 72 79 85];
dipfit.vol.c = [0.3300 1 0.0042 0.3300];
dipfit.vol.o = [0 0 0]; % origin is always zero (forward model can't handle non-zero origin)

% origin of best fitting sphere (in MNI space)
%--------------------------------------------------------------------------
dipfit.vol.o_sphere = fminsearch(@(x) d_elc_sphere(x, sen_reg, dipfit.vol.r(4)),[0 0 0]);

% compute POL->MNI coordinate transformation, following definition of
% CTF origin, and axes
%--------------------------------------------------------------------------

% origin in MNI space
%--------------------------------------------------------------------------
POLorigin = mean(fid_reg(2:3, :));
POLvx = (fid_reg(1, :) - POLorigin)/norm((fid_reg(1, :) - POLorigin));
POLvz = cross(POLvx, fid_reg(2, :) - fid_reg(3, :));
POLvz = POLvz/norm(POLvz);  % Vector normal to the NLR plane, pointing upwards
POLvy = cross(POLvz,POLvx); POLvy = POLvy/norm(POLvy);

% transformation POL->MNI. Must be the same as computed by registration
% (file RT).
%--------------------------------------------------------------------------
dipfit.Mpol2mni = [[POLvx' POLvy' POLvz' POLorigin']; [0 0 0 1]];
dipfit.MNIelc   = sen_reg;

% original coordinates
%--------------------------------------------------------------------------
elc = inv(dipfit.Mpol2mni)*[dipfit.MNIelc'; ones(1, size(dipfit.MNIelc, 1))];
elc = elc(1:3,:)';

% centre of sphere in POL space
%--------------------------------------------------------------------------
o   = inv(dipfit.Mpol2mni)*[dipfit.vol.o_sphere'; 1];
elc = elc - kron(ones(size(elc, 1), 1), o(1:3)');

% transformation matrix from MNI-space to sphere in POL space
%--------------------------------------------------------------------------
Msphere = eye(4);
Msphere(1:3,4) = -o(1:3);
dipfit.Mpol2sphere = Msphere; % in pol-space: translation of origin to zero
dipfit.Mmni2polsphere = Msphere*inv(dipfit.Mpol2mni); % from MNI-space to pol space (plus translation to origin)

% projecting channels to outer sphere
%--------------------------------------------------------------------------
dist         = sqrt(sum((elc).^2,2));
dipfit.elc   = dipfit.vol.r(4) * elc ./[dist dist dist];

return


% used for fitting sphere to sensors
%==========================================================================
function x = d_elc_sphere(x, elc, r)

% returns sum of squares of distances between sensors elc and sphere.
Nelc = size(elc, 1);
d    = sqrt(sum((elc - repmat(x, Nelc, 1)).^2, 2));
x    = sum((d - r).^2);


=======
function DCM = spm_dcm_erp_dipfit(DCM)
% prepares structures for ECD forward model (both EEG and MEG)
% FORMAT DCM = spm_dcm_erp_dipfit(DCM)
% requires:
% 
%       DCM.xY.Dfile
%       DCM.xY.Ic
%       DCM.M.dipfit.sensorfile - 'ECD'
%       DCM.M.dipfit.type       - 'ECD (EEG)','ECD (MEG)','Imaging','LFP'
%       DCM.M.dipfit.L.pos or DCM.Lpos
%__________________________________________________________________________
% Stefan Kiebel, Karl Friston
% $Id: spm_dcm_erp_dipfit.m 1040 2007-12-21 20:28:30Z karl $

% Get data filename and good channels
%--------------------------------------------------------------------------
try
    DCM.xY.Dfile;
    DCM.M.dipfit.Ic = DCM.xY.Ic;
catch
    errordlg('Please specify data');
    error('')
end

% Get source locations
%--------------------------------------------------------------------------
try
    DCM.M.dipfit.L.pos = DCM.Lpos;
catch
    try
        DCM.Lpos = DCM.M.dipfit.L.pos;
    catch
        errordlg({'Please specify source locations',...
                  'in DCM.Lpos'})
    end
end


% Polhemus file skip
%==========================================================================
try
    Pol_skip = DCM.options.Pol_skip;
catch
    Pol_skip = 2;
    DCM.options.Pol_skip = Pol_skip;
end
%==========================================================================


% D - SPM data structure
%==========================================================================
D       = spm_eeg_ldata(DCM.xY.Dfile);

switch DCM.options.type

    % EEG - ECD
    %----------------------------------------------------------------------
    case 1
    if strcmp(D.modality,'EEG')

        % try to get sen_reg and fid_reg from D.inv
        %------------------------------------------------------------------
        try
            sen_reg = D.inv{1}.datareg.sens_coreg;
            fid_reg = D.inv{1}.datareg.fid_coreg;
        catch

            % Get sensor file
            %--------------------------------------------------------------
            try
                sensorfile = DCM.M.dipfit.sensorfile;
            catch
                [f,p]      = uigetfile({'*.mat;*.pol'},'select sensor file');
                sensorfile = fullfile(p,f);
                DCM.M.dipfit.sensorfile = sensorfile;
            end

            % Polyhmus file
            %--------------------------------------------------------------
            if strcmpi('pol', spm_str_manip(sensorfile, 'e'))
                
                % returns coordinates of fiducials and sensors (filenames)
                % coordinates are in CTF coordinate system
                %----------------------------------------------------------
                [fid, sens] = spm_eeg_inv_ReadPolhemus(sensorfile, Pol_skip);
                
            else
 
                % mat-file must contain 2 matrices with fid and sensor coordinates
                %----------------------------------------------------------
                tmp  = load(sensorfile, '-mat');
                fid  = tmp.fid;
                sens = tmp.sens;
            end

            % Use approximated MRI fiducials (MNI space)... for now
            %--------------------------------------------------------------
            mni_fid = [[0   86 -42];...
                      [-84 -18 -55];...
                      [ 84 -18 -55]];


            % EEG coordinates in MNI space
            % coregistration between EEG and MRI fiducials
            %--------------------------------------------------------------
            [eeg2mri,sen_reg,fid_reg] = spm_eeg_inv_datareg(sens, fid, mni_fid);

            try
                Ic = DCM.M.dipfit.Ic;
            catch
                D  = spm_eeg_ldata(DCM.xY.Dfile);
                Ic = setdiff(D.channels.eeg, D.channels.Bad);
                DCM.M.dipfit.Ic = Ic;
            end
            sen_reg  = sen_reg(Ic, :);

        end

        % evaluate dipfit and save in M
        %------------------------------------------------------------------
        DCM.M.dipfit       = spm_dipfit(DCM.M.dipfit,sen_reg,fid_reg);
        DCM.M.dipfit.type  = 'ECD (EEG)';
    


    % MEG - ECD
    %----------------------------------------------------------------------
    else

        % not much to do because the sensors location/orientations were
        % already read at the time of conversion.
        %------------------------------------------------------------------
        DCM.M.dipfit.vol.r = [85];
        DCM.M.dipfit.vol.o = [0 0 0];
        try
            % done by coreg in source reconstruction GUI
            DCM.M.grad         = D.inv{1}.datareg.grad_coreg;            
        catch
            DCM.M.grad         = D.channels.grad;
        end
        DCM.M.dipfit.type  = 'ECD (MEG)';
        
    end
        
    % Imaging (distributed source reconstruction)
    %----------------------------------------------------------------------
    case 2

        % Load Gain or Lead field matrix
        %------------------------------------------------------------------
        try
            DCM.val = D.val;
            gainmat = D.inv{D.val}.forward.gainmat;
            try
                L       = load(gainmat);
            catch
                [p,f]   = fileparts(gainmat);
                L       = load(f);
                gainmat = fullfile(pwd,f);
            end
            name    = fieldnames(L);
            L       = sparse(getfield(L, name{1}));
        catch
            errordlg({'Please create and save a foward model',...
                      'using spm_eeg_inv_imag_api'})
            error('')
        end

        % centers
        %------------------------------------------------------------------
        xyz = DCM.M.dipfit.L.pos;
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
end

return


function dipfit = spm_dipfit(dipfit,sen_reg,fid_reg)
%==========================================================================
% computes dipfit model for ECD, DCM

% fit origin of outer sphere, with fixed radius, to electrodes
%--------------------------------------------------------------------------
dipfit.vol.r = [71 72 79 85];
dipfit.vol.c = [0.3300 1 0.0042 0.3300];
dipfit.vol.o = [0 0 0]; % origin is always zero (forward model can't handle non-zero origin)

% origin of best fitting sphere (in MNI space)
%--------------------------------------------------------------------------
dipfit.vol.o_sphere = fminsearch(@(x) d_elc_sphere(x, sen_reg, dipfit.vol.r(4)),[0 0 0]);

% compute POL->MNI coordinate transformation, following definition of
% CTF origin, and axes
%--------------------------------------------------------------------------

% origin in MNI space
%--------------------------------------------------------------------------
POLorigin = mean(fid_reg(2:3, :));
POLvx = (fid_reg(1, :) - POLorigin)/norm((fid_reg(1, :) - POLorigin));
POLvz = cross(POLvx, fid_reg(2, :) - fid_reg(3, :));
POLvz = POLvz/norm(POLvz);  % Vector normal to the NLR plane, pointing upwards
POLvy = cross(POLvz,POLvx); POLvy = POLvy/norm(POLvy);

% transformation POL->MNI. Must be the same as computed by registration
% (file RT).
%--------------------------------------------------------------------------
dipfit.Mpol2mni = [[POLvx' POLvy' POLvz' POLorigin']; [0 0 0 1]];
dipfit.MNIelc   = sen_reg;

% original coordinates
%--------------------------------------------------------------------------
elc = inv(dipfit.Mpol2mni)*[dipfit.MNIelc'; ones(1, size(dipfit.MNIelc, 1))];
elc = elc(1:3,:)';

% centre of sphere in POL space
%--------------------------------------------------------------------------
o   = inv(dipfit.Mpol2mni)*[dipfit.vol.o_sphere'; 1];
elc = elc - kron(ones(size(elc, 1), 1), o(1:3)');

% transformation matrix from MNI-space to sphere in POL space
%--------------------------------------------------------------------------
Msphere = eye(4);
Msphere(1:3,4) = -o(1:3);
dipfit.Mpol2sphere = Msphere; % in pol-space: translation of origin to zero
dipfit.Mmni2polsphere = Msphere*inv(dipfit.Mpol2mni); % from MNI-space to pol space (plus translation to origin)

% projecting channels to outer sphere
%--------------------------------------------------------------------------
dist         = sqrt(sum((elc).^2,2));
dipfit.elc   = dipfit.vol.r(4) * elc ./[dist dist dist];

return


% used for fitting sphere to sensors
%==========================================================================
function x = d_elc_sphere(x, elc, r)

% returns sum of squares of distances between sensors elc and sphere.
Nelc = size(elc, 1);
d    = sqrt(sum((elc - repmat(x, Nelc, 1)).^2, 2));
x    = sum((d - r).^2);


>>>>>>> .r1038

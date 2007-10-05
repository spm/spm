function [varargout] = spm_eeg_inv_BSTfwdsol(varargin)

%=======================================================================
% FORMAT D = spm_eeg_inv_BSTparameters(D,val)
%
% This function defines the required options for running the BrainStorm
% forward solution (bst_headmodeler.m)
%
% D                - input struct
% (optional) fields of S:
% D                - filename of EEG/MEG mat-file
%
% Output:
% D                - EEG/MEG struct with filenames of Gain matrices)
%
% See also help lines in bst_headmodeler.m
%==========================================================================
% Copyright (C) 2005 Wellcome Department of Imaging Neuroscience

% Jeremie Mattout & Christophe Phillips
% $Id: spm_eeg_inv_BSTfwdsol.m 934 2007-10-05 12:36:29Z karl $

% Modified by Rik Henson to handle gradiometers (with two positions/orientations
% for component coils) and to allow sphere to be fit to other surfaces, eg
% inner skull rather than scalp				4/6/07

% initialise
%--------------------------------------------------------------------------
[D,val] = spm_eeg_inv_check(varargin{:});
OPTIONS = spm_bst_headmodeler;

% Forward solution
% NB: overlapping_spheres.m does not exist; so it has been removed as an option
%--------------------------------------------------------------------------
try
    OPTIONS.Method = D.inv{val}.forward.method;
catch
    if D.modality == 'MEG'
        str        = {'spherical model'};
        switch spm_input('Forward model','+1','m',str,[1]);
            case 1
                OPTIONS.Method = 'meg_sphere';
            case 2
                OPTIONS.Method = 'meg_os';
        end

    elseif D.modality == 'EEG'

        str        = {'3 {Berg}','or 1 sphere',};
        switch spm_input('Forward model','+1','b',str,[2 1]);
            case 1
                OPTIONS.Method = 'eeg_sphere';
            case 2
                OPTIONS.Method = 'eeg_3sphereBerg';
            case 3
                OPTIONS.Method = 'meg_os';
        end
    end
end

% Added by Rik to allow different options for sphere fitting
try
    sphere2fit = D.inv{val}.forward.sphere2fit;
catch
    if isfield(D.inv{val}.mesh,'tess_iskull')
        sphere2fit = 2;		% Default to inner skull
    else
        sphere2fit = 1;		% Cortex
    end
end
% COULD ADD OPTION TO FIT SPHERE TO POLHEMUS HEAD SHAPE!

% I/O functions
%--------------------------------------------------------------------------
Nvert                      = length(D.inv{val}.mesh.tess_ctx.vert);
[path,nam,ext]             = fileparts(D.inv{val}.mesh.sMRI);
OPTIONS.ImageGridBlockSize = Nvert + 1;
OPTIONS.Verbose            = 1;

% Head Geometry (create tesselation file)
%--------------------------------------------------------------------------
vert = D.inv{val}.mesh.tess_ctx.vert;
face = D.inv{val}.mesh.tess_ctx.face;

% normals
%--------------------------------------------------------------------------
normal = spm_eeg_inv_normals(vert,face);

% convert positions into m
%--------------------------------------------------------------------------
if max(max(vert)) > 1 % non m
    if max(max(vert)) < 20 % cm
        vert = vert/100;
    elseif max(max(vert)) < 200 % mm
        vert = vert/1000;
    end
end
Comment{1}      = 'Cortex Mesh';
Faces{1}        = face;
Vertices{1}     = vert';
Curvature{1}    = [];
VertConn{1}     = cell(Nvert,1);
for i = 1:Nvert
    [ii,jj]     = find(face == i);
    voi = setdiff(unique(face(ii,:)),i);
    voi = reshape(voi,1,length(voi));
    VertConn{1}{i} = voi;
end

if sphere2fit == 1
    % compute the best fitting sphere to CORTEX
    %----------------------------------------------------------------------
    [Center,Radius]    = spm_eeg_inv_BestFitSph(vert);
    OPTIONS.HeadCenter = Center;
    OPTIONS.Radii      = [.88 .93 1]*Radius;
end

if isfield(D.inv{val}.mesh,'tess_iskull')

    % convert positions into m
    %----------------------------------------------------------------------
    vert         = D.inv{val}.mesh.tess_iskull.vert;
    face         = D.inv{val}.mesh.tess_iskull.face;
    if max(max(vert)) > 1 % non m
        if max(max(vert)) < 20 % cm
            vert = vert/100;
        elseif max(max(vert)) < 200 % mm
            vert = vert/1000;
        end
    end
    Comment{2}   = 'Inner Skull Mesh';
    Faces{2}     = face;
    Vertices{2}  = vert';
    Curvature{2} = [];
    VertConn{2}  = [];
    indx         = 3;

    if sphere2fit == 2
        % compute the best fitting sphere to INNER SKULL
        %----------------------------------------------------------------------
        [Center,Radius]    = spm_eeg_inv_BestFitSph(vert);
        OPTIONS.HeadCenter = Center;
        OPTIONS.Radii      = [.88 .93 1]*Radius;
    end
else
    indx         = 2;
end

if isfield(D.inv{val}.mesh,'tess_scalp')

    % convert positions into m
    %----------------------------------------------------------------------
    vert         = D.inv{val}.mesh.tess_scalp.vert;
    face         = D.inv{val}.mesh.tess_scalp.face;
    if max(max(vert)) > 1 % non m
        if max(max(vert)) < 20 % cm
            vert = vert/100;
        elseif max(max(vert)) < 200 % mm
            vert = vert/1000;
        end
    end
    Comment{indx}   = 'Scalp Mesh';
    Faces{indx}     = face;
    Vertices{indx}  = vert';
    Curvature{indx} = [];
    VertConn{indx}  = [];

    OPTIONS.Scalp.iGrid    = indx;

    if sphere2fit == 3
        % compute the best fitting sphere to SCALP
        %----------------------------------------------------------------------
        [Center,Radius]    = spm_eeg_inv_BestFitSph(vert);
        OPTIONS.HeadCenter = Center;
        OPTIONS.Radii      = [.88 .93 1]*Radius;
    end
end

OPTIONS.Cortex.ImageGrid.Comment   = Comment;
OPTIONS.Cortex.ImageGrid.Curvature = Curvature;
OPTIONS.Cortex.ImageGrid.Faces     = Faces;
OPTIONS.Cortex.ImageGrid.VertConn  = VertConn;
OPTIONS.Cortex.ImageGrid.Vertices  = Vertices;

% Sensor Information
%--------------------------------------------------------------------------
OPTIONS.ChannelType = D.modality;
sens     = D.inv{val}.datareg.sens_coreg';
%if length(sens) ~= size(sens,2) & length(sens) > 3
%    sens = sens';
%end

% convert positions into m (!!!)
%--------------------------------------------------------------------------
if max(max(sens)) > 1 % non m
    if max(max(sens)) < 20 % cm
        sens = sens/100;
    elseif max(max(sens)) < 200 % mm
        sens = sens/1000;
    end
end

% sensor orientations
%--------------------------------------------------------------------------
if (D.modality == 'MEG')
    orientation  = D.inv{val}.datareg.sens_orient_coreg';
end

ncoil = size(sens,1)/3;	% = 2 at least one gradiometer present (see spm_bst_headmodeler)

if ncoil > 1
    OPTIONS.HeadCenter = repmat(OPTIONS.HeadCenter,1,ncoil);
end

for i = 1:length(sens)
    Channel(i) = struct('Loc',[],'Orient',[],'Comment','','Weight',[],'Type','','Name','');

    Channel(i).Loc = reshape(sens(:,i),3,ncoil);

    if exist('orientation') == 1
        Channel(i).Orient = reshape(orientation(:,i),3,ncoil);
    else
        Channel(i).Orient = [];
    end

    if isfield(D.channels,'Weight')
        Channel(i).Weight  = D.channels.Weight(i,:);
    elseif ncoil == 1
        Channel(i).Weight  = 1;
    % (currently only handles first-order gradiometers, ie two coils)
    %----------------------------------------------------------------------
    else
        % this is a gradiometer
        if all(isfinite(Channel(i).Loc(:,2)))	       
            Channel(i).Weight  = [1 -1];
        % magnetometer
        else
            Channel(i).Weight  = [1 0];
        end
    end

    if ncoil > 1
        if ~all(isfinite(Channel(i).Loc(:,2)))		% replace mag NaNs
            Channel(i).Loc(:,2) = Channel(i).Loc(:,1);
            Channel(i).Orient(:,2) = Channel(i).Orient(:,1);
        end
    end
    Channel(i).Comment = num2str(i);
    Channel(i).Type    = D.modality;
    Channel(i).Name    = [D.modality ' ' num2str(i)];
end
OPTIONS.Channel = Channel;


% Set OPTIONS: Source Model
%--------------------------------------------------------------------------
OPTIONS.SourceModel      = -1; % current dipole model
OPTIONS.Cortex.iGrid     =  1;
OPTIONS.GridLoc          = Vertices{1};
OPTIONS.GridOrient       = normal';
OPTIONS.ApplyGridOrient  = 1;
OPTIONS.VolumeSourceGrid = 0;
OPTIONS.SourceLoc        = OPTIONS.GridLoc;
OPTIONS.SourceOrient     = OPTIONS.GridOrient;

% Forward computation
%--------------------------------------------------------------------------
[G,Gxyz] = spm_bst_headmodeler(OPTIONS);

% Save
%--------------------------------------------------------------------------
D.inv{val}.forward.gainmat = fullfile(D.path,[nam '_SPMgainmatrix_' num2str(val) '.mat']);
D.inv{val}.forward.gainxyz = fullfile(D.path,[nam '_SPMgainmatxyz_' num2str(val) '.mat']);

if spm_matlab_version_chk('7.1') >=0
    save(D.inv{val}.forward.gainmat,'-V6','G');
    save(D.inv{val}.forward.gainxyz,'-V6','Gxyz');
else
    save(D.inv{val}.forward.gainmat,'G');
    save(D.inv{val}.forward.gainxyz,'Gxyz');
end

varargout{1} = D;
if nargout > 1
    varargout{2} = OPTIONS;
end

return

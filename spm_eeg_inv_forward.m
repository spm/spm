function D = spm_eeg_inv_forward(varargin)
% Compute M/EEG leadfield
% FORMAT D = spm_eeg_inv_forward(D,val)
%
% D                - input struct
% (optional) fields of S:
% D                - filename of EEG/MEG mat-file
%
% Output:
% D                - EEG/MEG struct with filenames of Gain matrices)
%
%__________________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% Jeremie Mattout & Christophe Phillips
% $Id: spm_eeg_inv_forward.m 3965 2010-07-01 17:16:17Z vladimir $

% initialise
%--------------------------------------------------------------------------
[D,val] = spm_eeg_inv_check(varargin{:});

nvol = numel(D.inv{val}.forward);

if numel(D.inv{val}.datareg) ~= nvol
    error('Separate coregistration is required for every modality');
end

Fgraph  = spm_figure('GetWin','Graphics'); figure(Fgraph); clf

spm('Pointer', 'Watch');drawnow;

for i = 1:nvol
    mesh = spm_eeg_inv_transform_mesh(D.inv{val}.datareg(i).fromMNI*D.inv{val}.mesh.Affine, D.inv{val}.mesh);

    switch D.inv{val}.forward(i).voltype
        case '3-Shell Sphere (experimental)'
            cfg                        = [];
            cfg.feedback               = 'yes';
            cfg.headshape(1) = export(gifti(mesh.tess_scalp),  'ft');
            cfg.headshape(2) = export(gifti(mesh.tess_oskull), 'ft');
            cfg.headshape(3) = export(gifti(mesh.tess_iskull), 'ft');

            % determine the convex hull of the brain, to determine the support points
            pnt  = mesh.tess_ctx.vert;
            tric = convhulln(pnt);
            sel  = unique(tric(:));

            % create a triangulation for only the support points
            cfg.headshape(4).pnt = pnt(sel, :);
            cfg.headshape(4).tri = convhulln(pnt(sel, :));

            vol  = ft_prepare_concentricspheres(cfg);
            vert = spm_eeg_inv_mesh_spherify(mesh.tess_ctx.vert, mesh.tess_ctx.face, 'shift', 'no');
            mesh.tess_ctx.vert = vol.r(1)*vert + repmat(vol.o, size(vert, 1), 1);
            modality = 'EEG';
        case 'EEG BEM'      
            [p, f] = fileparts(mesh.sMRI);
            volfile = fullfile(p, [f '_EEG_BEM.mat']);
            if ~exist(volfile, 'file')
                % Just to show the user that something's happening
                spm_progress_bar('Init', 1, 'preparing EEG BEM model'); drawnow;
                spm('Pointer', 'Watch');drawnow;
                vol = [];
                vol.cond   = [0.3300 0.0041 0.3300];
                vol.source = 1; % index of source compartment
                vol.skin   = 3; % index of skin surface
                % brain
                vol.bnd(1) = export(gifti(D.inv{val}.mesh.tess_iskull), 'ft');
                % skull
                vol.bnd(2) = export(gifti(D.inv{val}.mesh.tess_oskull), 'ft');
                % skin
                vol.bnd(3) = export(gifti(D.inv{val}.mesh.tess_scalp),  'ft');

                % create the BEM system matrix
                cfg = [];
                cfg.method = 'bemcp';
                vol = ft_prepare_bemmodel(cfg, vol);

                spm_progress_bar('Set', 1); drawnow;

                save(volfile, 'vol');

                spm_progress_bar('Clear');
                spm('Pointer', 'Arrow');drawnow;
            end
            vol = volfile;
            modality = 'EEG';
        case 'OpenMEEG BEM'           
            vol = [];
            vol.cond   = [0.3300 0.0041 0.3300];
            vol.source = 1; % index of source compartment
            vol.skin   = 3; % index of skin surface
            % brain
            vol.bnd(1) = export(gifti(D.inv{val}.mesh.tess_iskull), 'ft');
            % skull
            vol.bnd(2) = export(gifti(D.inv{val}.mesh.tess_oskull), 'ft');
            % skin
            vol.bnd(3) = export(gifti(D.inv{val}.mesh.tess_scalp),  'ft');

            cfg = [];
            cfg.method = 'openmeeg';
            vol = ft_prepare_bemmodel(cfg, vol);
            modality = 'EEG';
        case 'Single Sphere'
            cfg                        = [];
            cfg.feedback               = 'yes';
            cfg.grad                   = D.inv{val}.datareg(i).sensors;
            cfg.headshape              = mesh.tess_iskull.vert;
            cfg.singlesphere = 'yes';
            vol  = ft_prepare_localspheres(cfg);
            modality = 'MEG';
        case 'MEG Local Spheres'
            cfg                        = [];
            cfg.feedback               = 'yes';
            cfg.grad                   = D.inv{val}.datareg(i).sensors;
            cfg.headshape              = mesh.tess_scalp.vert;
            cfg.radius                 = 85;
            cfg.maxradius              = 200;
            vol  = ft_prepare_localspheres(cfg);
            modality = 'MEG';
        case  'Single Shell'
            vol = [];
            vol.bnd = export(gifti(mesh.tess_iskull), 'ft');
            vol.type = 'nolte';
            vol = ft_convert_units(vol, 'mm');
            modality = 'MEG';
        otherwise
            error('Unsupported volume model type');
    end

    D.inv{val}.forward(i).vol       = vol;
    D.inv{val}.forward(i).mesh      = mesh.tess_ctx;
    D.inv{val}.forward(i).modality  = modality;

    Fgraph  = spm_figure('GetWin','Graphics'); figure(Fgraph); clf
end

% This is to force recomputing the lead fields
try, D.inv{val} = rmfield(D.inv{val}, 'gainmat'); end

spm('Pointer', 'Arrow');drawnow;
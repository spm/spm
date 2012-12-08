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
%__________________________________________________________________________
% Copyright (C) 2008-2012 Wellcome Trust Centre for Neuroimaging

% Jeremie Mattout & Christophe Phillips
% $Id: spm_eeg_inv_forward.m 5103 2012-12-08 06:44:58Z vladimir $


%-Initialisation
%--------------------------------------------------------------------------
[D, val] = spm_eeg_inv_check(varargin{:});

if numel(D.inv{val}.datareg) ~= numel(D.inv{val}.forward)
    error('Separate coregistration is required for every modality');
end

Fgraph = spm_figure('GetWin','Graphics');
spm_figure('Clear',Fgraph);
spm('Pointer', 'Watch');
if isempty(Fgraph) || spm('CmdLine'), graph = 'no'; else graph = 'yes'; end

for i = 1:numel(D.inv{val}.forward)
    mesh = spm_eeg_inv_transform_mesh(D.inv{val}.datareg(i).fromMNI*D.inv{val}.mesh.Affine, D.inv{val}.mesh);
    
    switch D.inv{val}.forward(i).voltype
        case '3-Shell Sphere (experimental)'
            cfg              = [];
            cfg.feedback     = graph;
            cfg.showcallinfo           = 'no';
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
            
            cfg.method = 'concentricspheres';
            
            vol  = ft_prepare_headmodel(cfg);
            vert = spm_eeg_inv_mesh_spherify(mesh.tess_ctx.vert, mesh.tess_ctx.face, 'shift', 'no');
            mesh.tess_ctx.vert = vol.r(1)*vert + repmat(vol.o, size(vert, 1), 1);
            modality = 'EEG';
            
        case 'EEG BEM'
            volfile = spm_file(mesh.sMRI, 'suffix','_EEG_BEM', 'ext','mat');
            if ~exist(volfile, 'file')
                spm_progress_bar('Init', 1, 'preparing EEG BEM model');
                
                vol        = [];
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
                cfg        = [];
                cfg.method = 'bemcp';
                cfg.showcallinfo = 'no';
                vol = ft_prepare_headmodel(cfg, vol);

                spm_progress_bar('Set', 1); drawnow;

                save(volfile, 'vol');
                
                spm_progress_bar('Clear');
                spm('Pointer', 'Arrow');
            end
            vol = volfile;
            modality = 'EEG';
            
        case 'OpenMEEG BEM'
            vol        = [];
            vol.cond   = [0.3300 0.0041 0.3300];
            vol.source = 1; % index of source compartment
            vol.skin   = 3; % index of skin surface
            % brain
            vol.bnd(1) = export(gifti(D.inv{val}.mesh.tess_iskull), 'ft');
            % skull
            vol.bnd(2) = export(gifti(D.inv{val}.mesh.tess_oskull), 'ft');
            % skin
            vol.bnd(3) = export(gifti(D.inv{val}.mesh.tess_scalp),  'ft');
            
            cfg        = [];
            cfg.method = 'openmeeg';
            cfg.showcallinfo = 'no';
            vol = ft_prepare_headmodel(cfg, vol);
            modality = 'EEG';
        case 'Single Sphere'
            cfg                        = [];
            cfg.feedback               = 'yes';
            cfg.showcallinfo           = 'no';
            cfg.grad                   = D.inv{val}.datareg(i).sensors;
            cfg.headshape              = export(gifti(mesh.tess_scalp), 'ft');
            cfg.method                 = 'singlesphere';
            vol                        = ft_prepare_headmodel(cfg);
            modality                   = 'MEG';
        case 'MEG Local Spheres'
            cfg                        = [];
            cfg.feedback               = 'yes';
            cfg.showcallinfo           = 'no';
            cfg.grad                   = D.inv{val}.datareg(i).sensors;
            cfg.headshape              = export(gifti(mesh.tess_scalp), 'ft');
            cfg.radius                 = 85;
            cfg.maxradius              = 200;
            cfg.method                 = 'localspheres';
            vol  = ft_prepare_headmodel(cfg);
            modality = 'MEG';
        case  'Single Shell'
            vol              = [];
            vol.bnd          = export(gifti(mesh.tess_iskull), 'ft');
            vol.type         = 'singleshell';
            vol              = ft_convert_units(vol, 'mm');
            modality         = 'MEG';
            
        otherwise
            error('Unsupported volume model type');
    end
    
    D.inv{val}.forward(i).vol      = vol;
    D.inv{val}.forward(i).mesh     = mesh.tess_ctx;
    D.inv{val}.forward(i).modality = modality;
    
    spm_figure('Clear',Fgraph);
end

% This is to force recomputing the lead fields
try, D.inv{val} = rmfield(D.inv{val}, 'gainmat'); end

spm('Pointer', 'Arrow');

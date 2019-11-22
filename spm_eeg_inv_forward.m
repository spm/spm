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
% Copyright (C) 2008-2018 Wellcome Trust Centre for Neuroimaging

% Jeremie Mattout & Christophe Phillips
% $Id: spm_eeg_inv_forward.m 7702 2019-11-22 11:32:26Z guillaume $


SVNrev = '$Rev: 7702 $';

%-Startup
%--------------------------------------------------------------------------
spm('FnBanner', mfilename, SVNrev);

%-Initialisation
%--------------------------------------------------------------------------
[D, val] = spm_eeg_inv_check(varargin{:});

if numel(D.inv{val}.datareg) ~= numel(D.inv{val}.forward)
    error('Separate coregistration is required for every modality.');
end

Fgraph = spm_figure('FindWin','Graphics');
spm_figure('Clear',Fgraph);
spm('Pointer', 'Watch');
if isempty(Fgraph) || spm('CmdLine'), graph = 'no'; else, graph = 'yes'; end

for i = 1:numel(D.inv{val}.forward)
    M    = D.inv{val}.datareg(i).fromMNI*D.inv{val}.mesh.Affine;
    
    M    = diag([1e-3 1e-3 1e-3 1])*M; % convert to m
    
    mesh = spm_eeg_inv_transform_mesh(M, D.inv{val}.mesh);
    
    mesh_correction = [];
    
    sens = D.inv{val}.datareg(i).sensors;
    
    if isequal(D.inv{val}.datareg(i).modality, 'MEG')
        sens = ft_datatype_sens(sens, 'amplitude', 'T', 'distance', 'm');
    else
        sens = ft_datatype_sens(sens, 'amplitude', 'V', 'distance', 'm');
    end
        
    switch D.inv{val}.forward(i).voltype
        case 'EEG interpolated'
            vol = D.inv{val}.forward(i).vol;
            modality = 'EEG';
        case '3-Shell Sphere'
            cfg              = [];
            cfg.feedback     = graph;
            cfg.siunits      = 'yes';
            cfg.showcallinfo = 'no';
          
            headshape(1) = export(gifti(mesh.tess_scalp),  'ft');
            headshape(2) = export(gifti(mesh.tess_oskull), 'ft');
            headshape(3) = export(gifti(mesh.tess_iskull), 'ft');
            
            % determine the convex hull of the brain, to determine the support points
            pnt  = mesh.tess_ctx.vert;
            tric = convhulln(pnt);
            sel  = unique(tric(:));
            
            % create a triangulation for only the support points
            headshape(4).pnt = pnt(sel, :);
            headshape(4).tri = convhulln(pnt(sel, :));
            
            cfg.method = 'concentricspheres';
            
            vol  = ft_prepare_headmodel(cfg, headshape);
            
            cfg = [];
            cfg.headmodel = vol;
            cfg.grid.pos  = mesh.tess_ctx.vert;
            cfg.spherify  = 'yes';
            gridsphere    = ft_prepare_sourcemodel(cfg);
            
            mesh_correction    = rmfield(cfg, {'headmodel', 'grid'});
            
            mesh.tess_ctx.vert = gridsphere.pos;
            modality = 'EEG';
            
        case 'EEG BEM'                        
            volfile = spm_file(mesh.sMRI, 'suffix','_EEG_BEM', 'ext','mat');
            vol = [];
            
            if exist(volfile, 'file')
                vol = ft_read_headmodel(volfile);
                if ~isfield(vol, 'unit') || ~isequal(vol.unit, 'm')
                    vol = [];
                end
            end
                
            if isempty(vol)
                
                vol.cond   = [0.3300 0.0041 0.3300];
                vol.source = 1; % index of source compartment
                vol.skin   = 3; % index of skin surface
                % brain
                vol.bnd(1) = export(gifti(mesh.tess_iskull), 'ft');
                % skull
                vol.bnd(2) = export(gifti(mesh.tess_oskull), 'ft');
                % skin
                vol.bnd(3) = export(gifti(mesh.tess_scalp),  'ft');
                
                % create the BEM system matrix
                cfg        = [];
                cfg.method = 'bemcp';
                cfg.showcallinfo = 'no';
                cfg.siunits      = 'yes';
                vol = ft_prepare_headmodel(cfg, vol);
                
                spm_progress_bar('Set', 1);
                
                save(volfile, 'vol', spm_get_defaults('mat.format'));
                
                spm_progress_bar('Clear');
                spm('Pointer', 'Arrow');
            end
            
            cfg = [];
            cfg.headmodel = vol;
            cfg.grid.pos = mesh.tess_ctx.vert;         
            cfg.moveinward = 6e-3; %move to empirically determined BEM safe zone
            gridcorrect = ft_prepare_sourcemodel(cfg);
            
            mesh_correction    = rmfield(cfg, {'headmodel', 'grid'});
            
            mesh.tess_ctx.vert = gridcorrect.pos;
            
            vol = volfile;
            modality = 'EEG';
            
        case 'OpenMEEG BEM'
            vol        = [];
            vol.cond   = [0.3300 0.0041 0.3300];
            vol.source = 1; % index of source compartment
            vol.skin   = 3; % index of skin surface
            % brain
            vol.bnd(1) = export(gifti(mesh.tess_iskull), 'ft');
            % skull
            vol.bnd(2) = export(gifti(mesh.tess_oskull), 'ft');
            % skin
            vol.bnd(3) = export(gifti(mesh.tess_scalp),  'ft');
            
            cfg                         = [];
            cfg.method                 = 'openmeeg';
            cfg.siunits                = 'yes';
            cfg.showcallinfo           = 'no';         
            vol                        = ft_prepare_headmodel(cfg, vol);
            
            cfg                        = [];
            cfg.vol                    = vol;
            cfg.grid.pos               = mesh.tess_ctx.vert;          
            cfg.moveinward             = 6e-3; % smaller shift might suffice for OpenMEEG
            gridcorrect                = ft_prepare_sourcemodel(cfg);
            
            mesh_correction            = rmfield(cfg, {'vol', 'grid'});
            
            mesh.tess_ctx.vert         = gridcorrect.pos;            
            
            modality = 'EEG';
        case 'Single Sphere'
            cfg                        = [];
            cfg.feedback               = 'yes';
            cfg.showcallinfo           = 'no';
            cfg.grad                   = D.inv{val}.datareg(i).sensors;            
            cfg.method                 = 'singlesphere';
            cfg.siunits                = 'yes';
            
            headshape                  = export(gifti(mesh.tess_scalp), 'ft');
            
            vol                        = ft_prepare_headmodel(cfg, headshape);
            modality                   = 'MEG';
        case 'MEG Local Spheres'
            cfg                        = [];
            cfg.feedback               = 'yes';
            cfg.showcallinfo           = 'no';
            cfg.grad                   = sens;           
            cfg.method                 = 'localspheres';
            cfg.siunits                = 'yes';
            
            headshape                  = export(gifti(mesh.tess_scalp), 'ft');
            vol                        = ft_prepare_headmodel(cfg, headshape);
            modality                   = 'MEG';
        case  'Single Shell'
            cfg                        = [];
            cfg.feedback               = 'yes';
            cfg.showcallinfo           = 'no';
            cfg.grad                   = sens;           
            cfg.method                 = 'singleshell';
            cfg.siunits                = 'yes';
            
            headshape                  = export(gifti(mesh.tess_iskull), 'ft');
            
            vol                        = ft_prepare_headmodel(cfg, headshape);
            modality                   = 'MEG';
            
        otherwise
            error('Unsupported volume model type.');
    end
    
    D.inv{val}.forward(i).vol             = vol;
    D.inv{val}.forward(i).mesh            = mesh.tess_ctx;
    D.inv{val}.forward(i).mesh_correction = mesh_correction;
    D.inv{val}.forward(i).modality        = modality;   
    D.inv{val}.forward(i).siunits         = 1;    
    
    D.inv{val}.forward(i).sensors  = sens;  
        
    D.inv{val}.forward(i).toMNI    = D.inv{val}.datareg(i).toMNI*diag([1e3 1e3 1e3 1]);
    D.inv{val}.forward(i).fromMNI  = diag([1e-3 1e-3 1e-3 1])*D.inv{val}.datareg(i).fromMNI;
    
    spm_figure('Clear',Fgraph);
end

% This is to force recomputing the lead fields
try, D.inv{val} = rmfield(D.inv{val}, 'gainmat'); end

fprintf('%-40s: %30s\n','Completed',spm('time'));                       %-#
spm('Pointer', 'Arrow');

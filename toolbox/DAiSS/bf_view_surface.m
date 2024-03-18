function res = bf_view_surface(BF, S)
% Diplay surface plot of DAISS output results
%__________________________________________________________________________

% George O'Neill
% Copyright (C) 2022-2023 Wellcome Centre for Human Neuroimaging


if nargin == 0
    
    imageno           = cfg_entry;
    imageno.tag       = 'imageno';
    imageno.name      = 'Image Number';
    imageno.strtype   = 'r';
    imageno.num       = [1 1];
    imageno.val       = {[1]};
    imageno.help      = {'Which image from output to show'};
    
    cmap            = cfg_entry;
    cmap.tag        = 'cmap';
    cmap.name       = 'Colormap';
    cmap.strtype    = 's';
    cmap.val        = {'jet'};
    cmap.help       = {'Colormap used for data'};
    
    cbar            = cfg_menu;
    cbar.tag        = 'cbar';
    cbar.name       = 'Colorbar';
    cbar.help       = {['Adds a colorbar to plot']};
    cbar.labels     = {'yes','no'};
    cbar.values     = {true, false};
    cbar.val        = {true};
    
    dock            = cfg_menu;
    dock.tag        = 'dock';
    dock.name       = 'Dock figure into SPM window';
    dock.help       = {['Plots figure in SPM Figure window, set to NO for new figure']};
    dock.labels     = {'yes','no'};
    dock.values     = {true, false};
    dock.val        = {true};
    
    inflate            = cfg_menu;
    inflate.tag        = 'inflate';
    inflate.name       = 'Inflate the cortical surface';
    inflate.help       = {['Inflates the surfaces, handy for viewing obfuscated areas']};
    inflate.labels     = {'yes','no'};
    inflate.values     = {true, false};
    inflate.val        = {false};
    
    threshold           = cfg_entry;
    threshold.tag       = 'threshold';
    threshold.name      = 'Threshold fraction';
    threshold.strtype   = 'r';
    threshold.num       = [1 1];
    threshold.val       = {[0]};
    threshold.help      = {'Threshold to a fraction of maximum (between 0 and 1)'};
    
    surface           = cfg_branch;
    surface.tag       = 'surface';
    surface.name      = 'Cortical Surface';
    surface.val       = {imageno,cmap,cbar,dock,inflate,threshold};
    
    res             = surface;
    
    return
elseif nargin < 2
    error('Two input arguments are required');
end

% Imageno option may have been unspecified
if ~isfield(S,'imageno')
    S.imageno = 1;
end

% if BF is a path rather than stucture, import
if isa(BF,'string')
    BF = bf_load(BF);
end

% Dock option may have been unspecified
if ~isfield(S,'dock')
    S.dock = true;
end

% Map option may have been unspecified
if ~isfield(S,'cmap')
    S.cmap = 'jet';
end

% Bar option may have been unspecified
if ~isfield(S,'cbar')
    S.cbar = false;
end

% Inflate option may have been unspecified
if ~isfield(S,'inflate')
    S.inflate = false;
end

% Inflate option may have been unspecified
if ~isfield(S,'threshold')
    S.threshold = 0;
end

% switch BF.data.space
%     case 'mni'
tmp = BF.sources.mesh.individual;
source.faces = tmp.face;
source.vertices = tmp.vert;

o = BF.output.image(S.imageno).val;
% S.threshold = 0.1;
if S.threshold
    o(abs(o) < S.threshold*max(abs(o))) = 0;
end

if S.dock
    F = spm_figure('GetWin','DAiSS Suface Plot');clf;
    if ismac
        set(F,'renderer','zbuffer');
    else
        set(F,'renderer','OpenGL');
    end
else
    F = figure;clf
end
a = axes;

if S.inflate
    H = spm_mesh_render('Disp', spm_mesh_inflate(source), 'Parent', a);
    H = spm_mesh_render('Overlay', H, o);
else
    H = spm_mesh_render('Disp', source, 'Parent', a);
    H = spm_mesh_render('Overlay', H, o);
end

cmap = feval(S.cmap,256);

spm_mesh_render('Colormap', H, cmap);

if S.cbar
spm_mesh_render('ColorBar',H,'on');
end


set(F,'color','w');

res = F;

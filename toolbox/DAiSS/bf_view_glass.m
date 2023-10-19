function res = bf_view_glass(BF, S)
% Diplays glass brain plot of DAISS output results
%__________________________________________________________________________

% George O'Neill
% Copyright (C) 2020-2023 Wellcome Centre for Human Neuroimaging


if nargin == 0
    
    classic         = cfg_menu;
    classic.tag     = 'classic';
    classic.name    = 'Classic glass brain visualisation';
    classic.help    = {['Use the classic glass brain visualisation, '...
        'set to false if you want to try a new experimental version...']};
    classic.labels  = {'yes','no'};
    classic.values     = {true, false};
    classic.val        = {true};
    
    ndips           = cfg_entry;
    ndips.tag       = 'ndips';
    ndips.name      = 'Number of sources';
    ndips.strtype   = 'i';
    ndips.num       = [1 1];
    ndips.val       = {[512]};
    ndips.help      = {'Number of sources plotted, set to inf for all sources'};
    
    cmap            = cfg_entry;
    cmap.tag        = 'cmap';
    cmap.name       = 'Colormap';
    cmap.strtype    = 's';
    cmap.val        = {'gray'};
    cmap.help       = {'Colormap of the glass brain'};
    
    cbar            = cfg_menu;
    cbar.tag        = 'cbar';
    cbar.name       = 'Colourbar';
    cbar.help       = {['Add colourbar to plot']};
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
    
    threshold           = cfg_entry;
    threshold.tag       = 'threshold';
    threshold.name      = 'Threshold fraction';
    threshold.strtype   = 'r';
    threshold.num       = [1 1];
    threshold.val       = {[0]};
    threshold.help      = {'Threshold to a fraction of maximum (between 0 and 1)'};
    
    glass           = cfg_branch;
    glass.tag       = 'glass';
    glass.name      = 'Glass Brain';
    glass.val       = {classic,ndips,cmap,cbar,dock,threshold};
    
    res             = glass;
    
    return
elseif nargin < 2
    error('Two input arguments are required');
end

% if BF is a path rather than stucture, import
if isa(BF,'string')
    BF = bf_load(BF);
end

% cbar option may have been unspecified
if ~isfield(S,'cbar')
    S.dock = false;
end

% Dock option may have been unspecified
if ~isfield(S,'dock')
    S.dock = true;
end

if ~isfield(BF.sources, 'pos')
    error('Source space snafu, email george!')
end

% cbar option may have been unspecified
if ~isfield(S,'thresohld')
    S.threshsold = 0;
end

S.ndips = min(S.ndips,length(BF.sources.pos));
if iscell(S.cmap), S.cmap = cell2mat(S.cmap); end

pos = ft_warp_apply(BF.data.transforms.toMNI, BF.sources.pos);

X = BF.output.image(end).val;
[~, id] = sort(abs(X),'descend');
if S.threshold > 0
    ndips = sum(abs(X) >= S.threshold*max(abs(X)));
    id = id(1:ndips);
else
    id = id(1:S.ndips);
end

if S.dock
    F = spm_figure('GetWin','MEG');clf;
    if ismac
        set(F,'renderer','zbuffer');
    else
        set(F,'renderer','OpenGL');
    end
else
    F = figure;clf
end

if S.classic
    spm_mip(X(id),pos(id,:),6);
    colormap(S.cmap);
else
    opt = [];
    opt.brush = 2;
    opt.cmap = S.cmap;
    opt.colourbar = S.cbar;
    spm_glass(X(id),pos(id,:),opt);
end
axis image
set(F,'color','w');

res = F;

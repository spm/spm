function res = bf_view_glass(BF, S)
% Diplays glass brain plot of DAISS output results
% Copyright (C) 2020 Wellcome Trust Centre for Neuroimaging

% George O'Neill
% $Id: bf_view_glass.m 8213 2022-01-27 15:33:26Z george $

%--------------------------------------------------------------------------
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
    
    dock            = cfg_menu;
    dock.tag        = 'dock';
    dock.name       = 'Dock figure into SPM window';
    dock.help       = {['Plots figure in SPM Figure window, set to NO for new figure']};
    dock.labels     = {'yes','no'};
    dock.values     = {true, false};
    dock.val        = {true};
    
    glass           = cfg_branch;
    glass.tag       = 'glass';
    glass.name      = 'Glass Brain';
    glass.val       = {classic,ndips,cmap,dock};
    
    res             = glass;
    
    return
elseif nargin < 2
    error('Two input arguments are required');
end

% if BF is a path rather than stucture, import
if isa(BF,'string')
    BF = bf_load(BF);
end

% Dock option may have been unspecified
if ~isfield(S,'dock')
    S.dock = true;
end

if ~isfield(BF.sources, 'pos')
    error('Source space snafu, email george!')
end

S.ndips = min(S.ndips,length(BF.sources.pos));
if iscell(S.cmap); S.cmap = cell2mat(S.cmap); end

pos = ft_warp_apply(BF.data.transforms.toMNI, BF.sources.pos);

X = BF.output.image(end).val;
[~, id] = sort(X,'descend');
id = id(1:S.ndips);

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
    spm_glass(X(id),pos(id,:),opt);
end
axis image
set(F,'color','w');

res = F;

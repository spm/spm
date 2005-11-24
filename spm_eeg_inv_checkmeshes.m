function [varargout] = spm_eeg_inv_checkmeshes(varargin);

%=======================================================================
% Generate the tesselated surfaces of the inner-skull and scalp from binary volumes.
%
% FORMAT D = spm_eeg_inv_checkmeshes(S)
% Input:
% S		    - input data struct (optional)
% Output:
% D			- data struct including the new files and parameters
%
% FORMAT spm_eeg_inv_checkmeshes(Mcortex, Miskull, Mscalp);
% Input :
% Mcortex   - filename of the cortical mesh (.mat file containing faces and vertices variables)
% Miskull   - inner-skull mesh
% Mscalp    - scalp mesh
%
% FORMAT [h_ctx, h_skl, h_slp] = spm_eeg_inv_checkmeshes(Mcortex, Miskull, Mscalp);
% Input :
% Mcortex   - filename of the cortical mesh (.mat file containing faces and vertices variables)
% Miskull   - inner-skull mesh
% Mscalp    - scalp mesh
%
% Output :
% h_ctx     - handle for cortex mesh
% h_skl     - handle for inner-skull mesh
% h_slp     - handle for scalp mesh
%=======================================================================
% Copyright (C) 2005 Wellcome Department of Imaging Neuroscience

% Jeremie Mattout
% $Id: spm_eeg_inv_checkmeshes.m 312 2005-11-24 19:35:42Z jeremie $

spm_defaults

Nmeshes = [1 1];

if nargin <= 1
    
    try
        D = varargin{1};
    catch
        D = spm_select(1, '.mat', 'Select EEG/MEG mat file');
        D = spm_eeg_ldata(D);
    end

    if ~isfield(D,'inv')
        error(sprintf('no inverse structure has been created for this data set\n'));
    end
    
    val = length(D.inv);
    
    if isempty(D.inv{val}.mesh.tess_ctx)
        Mcortex = spm_select(1,'.mat','Select cortex mesh');
        D.inv{val}.mesh.tess_ctx = Mcortex;
    else
        Mcortex = D.inv{val}.mesh.tess_ctx;
    end
    
    if isempty(D.inv{val}.mesh.tess_iskull)
        try
            Miskull = spm_select(1,'.mat','Select skull mesh');
            D.inv{val}.mesh.tess_iskull = Miskull;
        catch
            Nmeshes(1) = 0;
        end
    else
        Miskull = D.inv{val}.mesh.tess_iskull;
    end
    
    if isempty(D.inv{val}.mesh.tess_scalp)
        try
            Mscalp = spm_select(1,'.mat','Select scalp mesh');
            D.inv{val}.mesh.tess_scalp = Mscalp;
        catch
            Nmeshes(2) = 0;
        end
    else
        Mscalp = D.inv{val}.mesh.tess_scalp;
    end
    
elseif nargin == 3
    
    Mcortex = varargin{1};
    Miskull = varargin{2};
    Mscalp  = varargin{3};
    
else
    
    error(sprintf('Wrong input arguments\n'));
    
end


F = findobj('Tag', 'Graphics');

if isempty(F)
    F = spm_figure;
end

figure(F);
clf

% Cortex Mesh Display
variabl = load(Mcortex);
face    = getfield(variabl,'face');
vert    = getfield(variabl,'vert');
h_ctx   = patch('vertices',vert,'faces',face,'EdgeColor','b','FaceColor','b');
hold on

% Inner-skull Mesh Display
if Nmeshes(1)
    variabl = load(Miskull);
    face    = getfield(variabl,'face');
    vert    = getfield(variabl,'vert');
    h_skl   = patch('vertices',vert,'faces',face,'EdgeColor','r','FaceColor','none');
end

% Scalp Mesh Display
if Nmeshes(2)
    variabl = load(Mscalp);
    face    = getfield(variabl,'face');
    vert    = getfield(variabl,'vert');
    h_slp   = patch('vertices',vert,'faces',face,'EdgeColor','k','FaceColor','none');
end


axis equal;
axis off;
set(gcf,'color','white');
view(-135,45);
zoom(1.5)
rotate3d
drawnow


if nargout == 1
    varargout{1} = D;
elseif nargout == 3
    varargout{1} = h_ctx;
    varargout{2} = h_skl;
    varargout{3} = h_slp;
else
    return
end    

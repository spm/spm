function [varargout] = spm_eeg_inv_checkmeshes(varargin);

%=======================================================================
% Generate the tesselated surfaces of the inner-skull and scalp from binary volumes.
%
% FORMAT D = spm_eeg_inv_getmeshes(S)
% Input:
% S		    - input data struct (optional)
% Output:
% D			- data struct including the new files and parameters
%
% FORMAT spm_eeg_inv_getmeshes(Mcortex, Miskull, Mscalp);
% Input :
% Mcortex   - filename of the cortical mesh (.mat file containing faces and vertices variables)
% Miskull   - inner-skull mesh
% Mscalp    - scalp mesh
%
% FORMAT [h_ctx, h_skl, h_slp] = spm_eeg_inv_getmeshes(Mcortex, Miskull, Mscalp);
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
% $Id$

spm_defaults


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
        Miskull = spm_select(1,'.mat','Select skull mesh');
        D.inv{val}.mesh.tess_iskull = Miskull;
    else
        Miskull = D.inv{val}.mesh.tess_iskull;
    end
    
    if isempty(D.inv{val}.mesh.tess_scalp)
        Mscalp = spm_select(1,'.mat','Select scalp mesh');
        D.inv{val}.mesh.tess_scalp = Mscalp;
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
name    = fieldnames(variabl);
face    = getfield(variabl,name{2});
vert    = getfield(variabl,name{3});
h_ctx   = patch('vertices',vert,'faces',face,'EdgeColor','b','FaceColor','b');
hold on

% Inner-skull Mesh Display
variabl = load(Miskull);
name    = fieldnames(variabl);
face    = getfield(variabl,name{2});
vert    = getfield(variabl,name{3});
h_skl   = patch('vertices',vert,'faces',face,'EdgeColor','r','FaceColor','none');

% Scalp Mesh Display
variabl = load(Mscalp);
name    = fieldnames(variabl);
face    = getfield(variabl,name{2});
vert    = getfield(variabl,name{3});
h_slp   = patch('vertices',vert,'faces',face,'EdgeColor','k','FaceColor','none');


axis equal;
axis off;
set(gcf,'color','white');
view(-135,45);
% cameramenu('noreset');
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
    
%     
% %=======================================================================
% function ph = spm_eeg_inv_DisplayMesh(vertices,faces,Ecolor,Fcolor)
% 
% Position(1,:) = [100*min(vertices(:,1)),100*mean(vertices(:,2)),mean(vertices(:,3))];
% Position(2,:) = [-100*min(vertices(:,1)),100*mean(vertices(:,2)),mean(vertices(:,3))];
% Position(3,:) = [100*mean(vertices(:,1)),100*min(vertices(:,2)),mean(vertices(:,3))];
% Position(4,:) = [100*mean(vertices(:,1)),-100*min(vertices(:,2)),mean(vertices(:,3))];
% 
% for i = 1:size(Position,1)
%     lh(i) = light('Position',Position(i,:),'Color',[1 1 1],'Style','infinite');
% end
% set(ph,'DiffuseStrength',.6,'SpecularStrength',0,'AmbientStrength',.4,'SpecularExponent',5);
% axis equal;
% axis off;
% set(gcf,'color','white');
% cameramenu
% view(-90,0);
% material dull
% 
% return
% %=======================================================================

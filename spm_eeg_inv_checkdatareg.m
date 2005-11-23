function spm_eeg_inv_checkdatareg(varargin);

%=======================================================================
% Generate the tesselated surfaces of the inner-skull and scalp from binary volumes.
%
% FORMAT [handles...] = spm_eeg_inv_getmeshes(S)
% Input:
% S		    - input data struct (optional)
%=======================================================================
% Copyright (C) 2005 Wellcome Department of Imaging Neuroscience

% Jeremie Mattout
% $Id$

spm_defaults


if nargin == 1
    
    try
        D = varargin{1};
    catch
        D = spm_select(1, '.mat', 'Select EEG/MEG mat file');
        D = spm_eeg_ldata(D);
    end
    
elseif nargin == 0

    D = spm_select(1, '.mat', 'Select EEG/MEG mat file');
    D = spm_eeg_ldata(D);

else
        
    error(sprintf('Wrong input arguments\n'));
    
end

if ~isfield(D,'inv')
    error(sprintf('no inverse structure has been created for this data set\n'));
end
    
val = length(D.inv); 

F = findobj('Tag', 'Graphics');

if isempty(F)
    F = spm_figure;
end

figure(F);
clf


% --- DISPLAY ANATOMY ---
if ~isempty(D.inv{val}.mesh.tess_ctx)
    Mcortex = D.inv{val}.mesh.tess_ctx;
else
    Mcortex = spm_select(1, '.mat', 'Select Cortex Mesh file');
end
if ~isempty(D.inv{val}.mesh.tess_iskull)
    Miskull = D.inv{val}.mesh.tess_iskull;
else
    Miskull = '';
end
if ~isempty(D.inv{val}.mesh.tess_scalp)
    Mscalp = D.inv{val}.mesh.tess_scalp;
else
    Mscalp = spm_select(1, '.mat', 'Select Scalp Mesh file');
end

% Cortex Mesh
variabl = load(Mcortex);
name    = fieldnames(variabl);
face    = getfield(variabl,name{2});
vert    = getfield(variabl,name{3});
h_ctx   = patch('vertices',vert,'faces',face,'EdgeColor','b','FaceColor','b');
axis equal
axis off
view(-135,45)
hold on

% Inner-skull Mesh
if ~isempty(Miskull)
    variabl = load(Miskull);
    name    = fieldnames(variabl);
    face    = getfield(variabl,name{2});
    vert    = getfield(variabl,name{3});
    h_skl   = patch('vertices',vert,'faces',face,'EdgeColor','r','FaceColor','none');
end

% Scalp Mesh
variabl = load(Mscalp);
name    = fieldnames(variabl);
face    = getfield(variabl,name{2});
vert    = getfield(variabl,name{3});
h_slp   = patch('vertices',vert,'faces',face,'EdgeColor','k','FaceColor','none');

% --- DISPLAY SETUP ---
if ~isempty(D.inv{val}.datareg.sens_coreg)
    Lsens = D.inv{val}.datareg.sens_coreg;
else
    Lsens = spm_select(1, '.mat', 'Select sensor file');
end
if ~isempty(D.inv{val}.datareg.fid_coreg)
    Lfid = D.inv{val}.datareg.fid_coreg;
else
    Lfid = '';
end
if ~isempty(D.inv{val}.datareg.fidmri)
    Lfidmri = D.inv{val}.datareg.fidmri;
else
    Lfidmri = '';
end

% Sensors (coreg.)
variabl = load(Lsens);
name    = fieldnames(variabl);
sensors = getfield(variabl,name{1});
h_sens  = plot3(sensors(:,1),sensors(:,2),sensors(:,3),'og');
set(h_sens,'MarkerFaceColor','g','MarkerSize',12,'MarkerEdgeColor','k');

% EEG fiducials or MEG coils (coreg.)
if ~isempty(Lfid)
    variabl = load(Lfid);
    name    = fieldnames(variabl);
    fid     = getfield(variabl,name{1});
    h_fid   = plot3(fid(:,1),fid(:,2),fid(:,3),'oc');
    set(h_fid,'MarkerFaceColor','c','MarkerSize',12,'MarkerEdgeColor','k');
end

% MRI fiducials
if ~isempty(Lfidmri)
    variabl = load(Lfidmri);
    name    = fieldnames(variabl);
    fidmri  = getfield(variabl,name{1});
    h_fidmr = plot3(fidmri(:,1),fidmri(:,2),fidmri(:,3),'dm');
    set(h_fidmr,'MarkerFaceColor','m','MarkerSize',12,'MarkerEdgeColor','k');
end

% cameramenu
zoom(1.5)
rotate3d
drawnow

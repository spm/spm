function spm_eeg_inv_checkdatareg(varargin)
% Display of the coregistred meshes and sensor locations in MRI space for
% quality check by eye.
% Fiducials which were used for rigid registration are also displayed
%
% FORMAT spm_eeg_inv_checkdatareg(D, val, ind)
%__________________________________________________________________________
% Copyright (C) 2005-2014 Wellcome Trust Centre for Neuroimaging

% Jeremie Mattout
% $Id: spm_eeg_inv_checkdatareg.m 7010 2017-02-07 18:15:30Z vladimir $


% Inputs
%--------------------------------------------------------------------------
[D,val] = spm_eeg_inv_check(varargin{:});

datareg = D.inv{val}.datareg;

if nargin < 3
    str = sprintf('%s|', datareg(:).modality);
    str = str(1:(end-1));    
    ind = spm_input('What to display?','+1', 'b',  str, 1:numel(D.inv{val}.datareg), 1);
else
    ind = varargin{3};
end

modality = datareg(ind).modality;
meegfid  = datareg(ind).fid_eeg;
mrifid   = datareg(ind).fid_mri;
mesh     = spm_eeg_inv_transform_mesh(datareg(ind).fromMNI*D.inv{val}.mesh.Affine, D.inv{val}.mesh);
sensors  = datareg(ind).sensors;

% SPM Graphics figure
%--------------------------------------------------------------------------
Fgraph   = spm_figure('GetWin','Graphics'); figure(Fgraph); clf
subplot(2,1,1)

%-Display anatomy
%==========================================================================
Mcortex = mesh.tess_ctx;
Miskull = mesh.tess_iskull;
Mscalp  = mesh.tess_scalp;

% Cortical Mesh
%--------------------------------------------------------------------------
face    = Mcortex.face;
vert    = Mcortex.vert;
h_ctx   = patch('vertices',vert,'faces',face,'EdgeColor','b','FaceColor','b');
hold on

% Inner-skull Mesh
%--------------------------------------------------------------------------
face    = Miskull.face;
vert    = Miskull.vert;
h_skl   = patch('vertices',vert,'faces',face,'EdgeColor','r','FaceColor','none');

% Scalp Mesh
%--------------------------------------------------------------------------
face    = Mscalp.face;
vert    = Mscalp.vert;
h_slp   = patch('vertices',vert,'faces',face,'EdgeColor',[1 .7 .55],'FaceColor','none');

%-Display setup
%==========================================================================
try   
    Lhsp    = meegfid.pnt;
    Lfidmri = mrifid.fid.pnt;
    Lfid    = meegfid.fid.pnt(1:size(Lfidmri, 1), :);
catch
    spm('alert!','Please coregister these data',mfilename);
    return
end

% Headshape locations
%--------------------------------------------------------------------------
if ~isempty(Lhsp)
    h_hsp   = plot3(Lhsp(:,1),Lhsp(:,2),Lhsp(:,3),'dm');
    set(h_hsp,'MarkerFaceColor','r','MarkerSize',4,'MarkerEdgeColor','r');
end

% Sensors (coreg.)
%--------------------------------------------------------------------------
try
    ft_plot_sens(sensors, 'chantype', unique(lower(D.chantype(D.indchantype(modality)))), 'edgecolor', [0 1 0]);
catch
    ft_plot_sens(sensors,'edgecolor', [0 1 0], 'coilshape', 'point', 'coil', true);
end
axis normal;
% 

% EEG fiducials or MEG coils (coreg.)
%--------------------------------------------------------------------------
h_fid   = plot3(Lfid(:,1),Lfid(:,2),Lfid(:,3),'oc');
set(h_fid,'MarkerFaceColor','c','MarkerSize',12,'MarkerEdgeColor','k');

% MRI fiducials
%--------------------------------------------------------------------------
h_fidmr = plot3(Lfidmri(:,1),Lfidmri(:,2),Lfidmri(:,3),'dm');
set(h_fidmr,'MarkerFaceColor','m','MarkerSize',12,'MarkerEdgeColor','k');

% Camera view
%--------------------------------------------------------------------------
axis image off
view(-135,45)
% cameratoolbar('setmode','orbit')
hold off

%-Display channels
%==========================================================================
subplot(2,1,2)

try
    [xy, label] = spm_eeg_project3D(sensors, modality);
    FontSize = 8;
catch
    xy = [0.25; 0.5];
    label = {'2D layout not available'};
    FontSize = 20;
end

% Channel names
%--------------------------------------------------------------------------
text(xy(1, :), xy(2, :), label,...
     'FontSize',FontSize,...
     'Color','r',...
     'FontWeight','bold')
  
axis equal off
cameratoolbar('setmode','orbit')

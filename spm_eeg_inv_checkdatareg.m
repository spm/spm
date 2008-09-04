function spm_eeg_inv_checkdatareg(varargin)
% Display of the coregistred meshes and sensor locations in MRI space for
% quality check by eye.
% Fiducials which were used for rigid registration are also displayed
%
% FORMAT spm_eeg_inv_checkdatareg(mesh, sensors)
%__________________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% Jeremie Mattout
% $Id: spm_eeg_inv_checkdatareg.m 2041 2008-09-04 13:39:40Z jean $

% SPM graphics figure
%--------------------------------------------------------------------------

%%

[D,val] = spm_eeg_inv_check(varargin{:});

switch D.inv{val}.modality
    case 'EEG'
        sens = D.inv{val}.datareg.sensors;
        sensorig = sens;
    case 'MEG'    
        cfg = [];
        cfg.style = '3d';
        cfg.rotate = 0;
        cfg.grad = D.inv{val}.datareg.sensors;
   
        lay = ft_prepare_layout(cfg);
                
        sens = [];
        sens.label = lay.label(:, 1);
        sens.pnt = lay.pos;
        sensorig = cfg.grad;
end
        
modality = D.inv{val}.modality;
meegfid = D.inv{val}.datareg.fid_eeg;
vol = D.inv{val}.forward.vol;
mrifid = D.inv{val}.datareg.fid_mri;
mesh = D.inv{val}.mesh;


Fgraph  = spm_figure('GetWin','Graphics'); figure(Fgraph); clf
subplot(2,1,1)

% Cortical Mesh
%--------------------------------------------------------------------------
face    = mesh.tess_ctx.face;
vert    = mesh.tess_ctx.vert;
h_ctx   = patch('vertices',vert,'faces',face,'EdgeColor','b','FaceColor','b');
hold on

% Scalp Mesh
%--------------------------------------------------------------------------
face    = mrifid.tri;
vert    = mrifid.pnt;
h_slp   = patch('vertices',vert,'faces',face,'EdgeColor',[0 0 0],'FaceColor','none');


% Inner volume Mesh
%--------------------------------------------------------------------------
face    = vol.bnd(end).tri;
vert    = vol.bnd(end).pnt;
h_vol   = patch('vertices',vert,'faces',face,'EdgeColor',[1 .7 .55],'FaceColor','none');


% --- DISPLAY SETUP ---
%==========================================================================
try
    Lsens   = sens.pnt;
    Lhsp    = meegfid.pnt;
    Lfidmri = mrifid.fid.pnt;
    Lfid    = meegfid.fid.pnt(1:size(Lfidmri, 1), :);
    Llabel = sens.label;
catch
    warndlg('please coregister these data')
    return
end


% EEG fiducials or MEG coils (coreg.)
%--------------------------------------------------------------------------
h_fid   = plot3(Lfid(:,1),Lfid(:,2),Lfid(:,3),'oc');
set(h_fid,'MarkerFaceColor','c','MarkerSize',12,'MarkerEdgeColor','k');

% MRI fiducials
%--------------------------------------------------------------------------
h_fidmr = plot3(Lfidmri(:,1),Lfidmri(:,2),Lfidmri(:,3),'dm');
set(h_fidmr,'MarkerFaceColor','m','MarkerSize',12,'MarkerEdgeColor','k');

% headshape locations
%--------------------------------------------------------------------------
if ~isempty(Lhsp)
    h_hsp   = plot3(Lhsp(:,1),Lhsp(:,2),Lhsp(:,3),'dm');
    set(h_hsp,'MarkerFaceColor','r','MarkerSize',4,'MarkerEdgeColor','r');
end

% Sensors (coreg.)
%--------------------------------------------------------------------------
h_sens  = plot3(Lsens(:,1),Lsens(:,2),Lsens(:,3),'og');
set(h_sens,'MarkerFaceColor','g','MarkerSize', 12,'MarkerEdgeColor','k');
% 
% camera view
%--------------------------------------------------------------------------
axis image off
view(-135,45)
% cameratoolbar('setmode','orbit')
hold off
zoom(5/3)

% DISPLAY CHANNELS
%==========================================================================
subplot(2,1,2)

[xy, label] = spm_eeg_project3D(sensorig, modality);

% Channel names
%--------------------------------------------------------------------------
text(xy(1, :), xy(2, :), label,...
     'FontSize',8,...
     'Color','r',...
     'FontWeight','bold')
  
axis equal off
cameratoolbar('setmode','orbit')

function spm_eeg_inv_checkdatareg(S)
% Display of the coregistred meshes and sensor locations in MRI space for
% quality check by eye.
% Fiducials which were used for rigid registration are also displayed
%
% FORMAT spm_eeg_inv_checkdatareg(mesh, sensors)
%__________________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% Jeremie Mattout
% $Id: spm_eeg_inv_checkdatareg.m 1650 2008-05-15 10:22:31Z vladimir $

% SPM graphics figure
%--------------------------------------------------------------------------
Fgraph  = spm_figure('GetWin','Graphics'); figure(Fgraph); clf
subplot(2,1,1)

% Cortical Mesh
%--------------------------------------------------------------------------
face    = S.mesh.tess_ctx.face;
vert    = S.mesh.tess_ctx.vert;
h_ctx   = patch('vertices',vert,'faces',face,'EdgeColor','b','FaceColor','b');
hold on

% Scalp Mesh
%--------------------------------------------------------------------------
face    = S.vol.bnd(1).tri;
vert    = S.vol.bnd(1).pnt;
h_slp   = patch('vertices',vert,'faces',face,'EdgeColor',[1 .7 .55],'FaceColor','none');


% --- DISPLAY SETUP ---
%==========================================================================
try
    Lsens   = S.sens.pnt;
    Lhsp    = S.meegfid.pnt;
    Lfidmri = S.mrifid.fid.pnt;
    Lfid    = S.meegfid.fid.pnt(1:size(Lfidmri, 1), :);
    Llabel = S.sens.label;
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
rotate3d on
hold off
zoom(5/3)

% DISPLAY CHANNELS
%==========================================================================
subplot(2,1,2)

% MRI fiducials
%--------------------------------------------------------------------------
plot3(Lfidmri(:,1),Lfidmri(:,2),Lfidmri(:,3),'dm',...
      'MarkerFaceColor','m',...
      'MarkerSize',12,...
      'MarkerEdgeColor','k');


% Channel names
%--------------------------------------------------------------------------
text(Lsens(:, 1),Lsens(:,2),Lsens(:,3), Llabel,...
     'FontSize',8,...
     'Color','r',...
     'FontWeight','bold')
axis equal off
axis([min(Lsens(:,1)), max(Lsens(:,1)), min(Lsens(:,2)),...
      max(Lsens(:,2)), min(Lsens(:,3)), max(Lsens(:,3))])
view(-140,70)
rotate3d on


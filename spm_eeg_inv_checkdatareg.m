function spm_eeg_inv_checkdatareg(varargin);

%=========================================================================
% Display of the coregistred meshes and sensor locations in MRI space for
% quality check by eye.
% Fiducials which were used for rigid registration are also displayed
%
% FORMAT spm_eeg_inv_checkdatareg(D,[val])
% Input:
% D		    - input data struct (optional)
%==========================================================================
% Copyright (C) 2005 Wellcome Department of Imaging Neuroscience

% Jeremie Mattout
% $Id: spm_eeg_inv_checkdatareg.m 850 2007-07-10 15:43:27Z rik $

% Minor change by Rik to handle sensors consisting of two gradiometer coils 5/6/07

% initialise
%--------------------------------------------------------------------------
[D,val] = spm_eeg_inv_check(varargin{:});

% SPM graphics figure
%--------------------------------------------------------------------------
Fgraph  = spm_figure('GetWin','Graphics'); figure(Fgraph); clf
subplot(2,1,1)

% --- DISPLAY ANATOMY ---
%==========================================================================
Mcortex = D.inv{val}.mesh.tess_ctx;
Miskull = D.inv{val}.mesh.tess_iskull;
Mscalp  = D.inv{val}.mesh.tess_scalp;

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


% --- DISPLAY SETUP ---
%==========================================================================
try
    Lsens   = D.inv{val}.datareg.sens_coreg;
    Lfid    = D.inv{val}.datareg.fid_coreg;
    Lhsp    = D.inv{val}.datareg.hsp_coreg;
    Lfidmri = D.inv{val}.datareg.fid_mri;
catch
    warndlg('please coregister these data')
    return
end

if size(Lsens,2)==6	% gradiometers, two coils
    for ch = 1:size(Lsens,1)
        if all(isfinite(Lsens(ch,4:6)))
          Lsens(ch,1:3) = (Lsens(ch,1:3)+Lsens(ch,4:6))/2;
	end
    end
    Lsens = Lsens(:,1:3);
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
h_hsp   = plot3(Lhsp(:,1),Lhsp(:,2),Lhsp(:,3),'dm');
set(h_hsp,'MarkerFaceColor','r','MarkerSize',4,'MarkerEdgeColor','r');

% Sensors (coreg.)
%--------------------------------------------------------------------------
h_sens  = plot3(Lsens(:,1),Lsens(:,2),Lsens(:,3),'og');
set(h_sens,'MarkerFaceColor','g','MarkerSize',12,'MarkerEdgeColor','k');

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
j       = setdiff(D.channels.eeg, D.channels.Bad);
text(Lsens(:, 1),Lsens(:,2),Lsens(:,3),D.channels.name(j),...
     'FontSize',8,...
     'Color','r',...
     'FontWeight','bold')
axis equal off
axis([min(Lsens(:,1)), max(Lsens(:,1)), min(Lsens(:,2)),...
      max(Lsens(:,2)), min(Lsens(:,3)), max(Lsens(:,3))])
view(-140,70)
rotate3d on

% display field
%--------------------------------------------------------------------------
disp(D.inv{val}.datareg)

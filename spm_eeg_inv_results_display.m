function spm_eeg_inv_results_display(D)
% Displays contrast of evoked responses and power
% FORMAT spm_eeg_inv_results_display((D)
%__________________________________________________________________________

%==========================================================================
Ndip  = 512; % Number of dipoles to display
%==========================================================================

% SPM data structure
%==========================================================================
model = D.inv{D.val};
try
    disp(model.contrast);
catch
    warndlg('please specify a [time-frequency] contrast')
    return
end

% inversion parameters
%--------------------------------------------------------------------------
Is   = model.inverse.Is;                          % Indices of ARD vertices
pst  = model.inverse.pst;                         % preistimulus tim (ms)
Nd   = model.inverse.Nd;                          % number of mesh dipoles

W    = model.contrast.W;
JW   = model.contrast.JW;
GW   = model.contrast.GW;

% sqrt(energy) (G) = abs(JW) for single trials
%--------------------------------------------------------------------------
G    = sqrt(sparse(Is,1,GW,Nd,1));

% display
%==========================================================================
Fgraph = spm_figure('GetWin','Graphics');
clf(Fgraph)
figure(Fgraph)
vert   = model.mesh.tess_mni.vert;

% display
%--------------------------------------------------------------------------
subplot(2,1,1)
[i j]  = sort(-G);
try
    j  = j(1:Ndip);
end
spm_mip(G(j),vert(j,:)',6);
axis image
title({'root mean square (energy)',sprintf('%i voxels',length(j))})

% contrast
%--------------------------------------------------------------------------
subplot(2,1,2)
plot(pst,W)
axis square
xlabel('PST {ms}')
ylabel('contrast')





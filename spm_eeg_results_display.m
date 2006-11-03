function spm_eeg_results_display(D)
% Displays contrast of evoked responses and power
% FORMAT spm_eeg_results_display((D)
%__________________________________________________________________________

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
IW   = model.contrast.IW;

% sqrt(energy) (G) = abs(JW) for single trials
%--------------------------------------------------------------------------
G    = sqrt(sparse(Is,1,GW,Nd,1));

% display
%==========================================================================
Fgraph   = spm_figure('GetWin','Graphics');
clf(Fgraph)
figure(Fgraph)
tess_ctx = model.mesh.tess_ctx;
load(tess_ctx)

% display
%--------------------------------------------------------------------------
subplot(2,3,1)
spm_eeg_inv_render(G,tess_ctx)
view([180 -90])
title('energy')

i     = find(G > max(G)/8);
subplot(2,1,2)
spm_mip(G(i),vert(i,:)',6);
axis image
title('root mean square (energy) ( > max/8)')

% contrast
%--------------------------------------------------------------------------
subplot(2,2,2)
plot(pst,W)
axis square
xlabel('PST {ms}')
ylabel('contrast')



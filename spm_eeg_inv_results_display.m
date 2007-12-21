function spm_eeg_inv_results_display(D)
% Displays contrast of evoked responses and power
% FORMAT spm_eeg_inv_results_display((D)
%__________________________________________________________________________

%==========================================================================
Ndip  = 256; % Number of dipoles to display
%==========================================================================

% SPM data structure
%==========================================================================
try, val = D.val; catch, val = 1; end
try, con = D.con; catch, con = 1; end
model = D.inv{D.val};
con   = min(con,length(model.inverse.J));
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
Ndip = min(Ndip,length(Is));

W    = model.contrast.W;
JW   = model.contrast.JW{con};
GW   = model.contrast.GW{con};

% sqrt(energy) (G) = abs(JW) for single trials
%--------------------------------------------------------------------------
G    = sqrt(sparse(Is,1,GW,Nd,1));

% display
%==========================================================================
Fgraph = spm_figure('GetWin','Graphics');
clf(Fgraph)
figure(Fgraph)

% get vertices (even if not normalised)
%--------------------------------------------------------------------------
try
	vert   = model.mesh.tess_mni.vert;
catch
	warndlg('Displaying on subject mesh - may not be in MNI space');
	vert   = model.mesh.tess_ctx.vert;
end

% display
%--------------------------------------------------------------------------
subplot(2,1,1)
[i j]  = sort(-G);
j      = j(1:Ndip);
spm_mip(G(j),vert(j,:)',6);
axis image

try
    str = sprintf('Energy (%s)',model.contrast.type);
catch
    str = 'Energy';
end
    
title({sprintf('Condition %d',con), str, sprintf('%i voxels',length(j))})

% contrast
%--------------------------------------------------------------------------
subplot(2,1,2)
plot(pst,W)
axis square
xlabel('PST {ms}')
ylabel('contrast')
drawnow





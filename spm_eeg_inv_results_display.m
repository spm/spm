function spm_eeg_inv_results_display(D)
% Displays contrast of evoked responses and power
% FORMAT spm_eeg_inv_results_display((D)
%__________________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% Karl Friston
% $Id: spm_eeg_inv_results_display.m 3976 2010-07-08 14:12:31Z karl $

%==========================================================================
Ndip  = 256; % Number of dipoles to display
%==========================================================================

% SPM data structure
%==========================================================================
try, val = D.val; catch, val = 1; end
try, con = D.con; catch, con = 1; end

if con == 0
    con = 1;
end

model = D.inv{D.val};
try
    con   = min(con,length(model.contrast.GW));
catch
    warndlg('please specify a [time-frequency] contrast')
    return
end

% inversion parameters
%--------------------------------------------------------------------------
Is    = model.inverse.Is;                         % Indices of ARD vertices
pst   = model.inverse.pst;                        % preistimulus tim (ms)
Nd    = model.inverse.Nd;                         % number of mesh dipoles
Ndip  = min(Ndip,length(Is));

try
    W = model.contrast.W{con};
catch
    W = model.contrast.W;
end
JW    = model.contrast.JW{con};
GW    = model.contrast.GW{con};

% just display the first trial (for trial-specific contrasts)
%--------------------------------------------------------------------------
if iscell(GW)
    GW = GW{1};
end

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
vert   = model.mesh.tess_mni.vert;

% display
%--------------------------------------------------------------------------
subplot(2,1,1)
[i j]  = sort(-G);
j      = j(1:Ndip);
spm_mip(G(j),vert(j,:)',6);
axis image

try
    if strcmp(model.contrast.type, 'trials')
        str = sprintf('Energy (%s)', 'first trial');
    else
        str = sprintf('Energy (%s)', model.contrast.type);
    end
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





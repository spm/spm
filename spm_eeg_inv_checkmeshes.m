function [h_ctx,h_skl,h_slp] = spm_eeg_inv_checkmeshes(varargin);

%=======================================================================
% Generate the tesselated surfaces of the inner-skull and scalp from binary volumes.
%
% FORMAT [h_ctx,h_skl,h_slp] = spm_eeg_inv_checkmeshes(S)
% Input:
% S		    - input data struct (optional)
% Output:
% h_ctx     - handle to cortex patch
% h_skl     - handle to skull patch
% h_slp     - handle to scalp patch
%
%==========================================================================
% Copyright (C) 2005 Wellcome Department of Imaging Neuroscience

% Jeremie Mattout
% $Id: spm_eeg_inv_checkmeshes.m 850 2007-07-10 15:43:27Z rik $


% initialise
%--------------------------------------------------------------------------
[D,val] = spm_eeg_inv_check(varargin{:});
try
    disp(D.inv{val}.mesh);
    Mctx  = D.inv{val}.mesh.tess_ctx;
    Mskl  = D.inv{val}.mesh.tess_iskull;
    Mslp  = D.inv{val}.mesh.tess_scalp;
catch
    warndlg('please create meshes')
    return
end

% SPM graphics figure
%--------------------------------------------------------------------------
Fgraph  = spm_figure('GetWin','Graphics'); spm_figure('Clear',Fgraph)

% Cortex Mesh Display
%--------------------------------------------------------------------------
h_ctx   = patch('vertices',Mctx.vert,'faces',Mctx.face,'EdgeColor','b','FaceColor','b');
hold on

% Inner-skull Mesh Display
%--------------------------------------------------------------------------
h_skl   = patch('vertices',Mskl.vert,'faces',Mskl.face,'EdgeColor','r','FaceColor','none');

% Scalp Mesh Display
%--------------------------------------------------------------------------
h_slp   = patch('vertices',Mslp.vert,'faces',Mslp.face,'EdgeColor',[1 .7 .55],'FaceColor','none');

axis image off;
view(-135,45);
rotate3d on
drawnow
hold off


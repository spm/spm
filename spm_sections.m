function spm_sections(xSPM,hReg,img)
% Rendering of regional effects [SPM{.}] on orthogonal sections
% FORMAT spm_sections(xSPM,hReg,img)
%
% xSPM  - structure containing details of excursion set (see spm_getSPM)
% hReg  - handle of MIP register
% img   - filename of background image
%__________________________________________________________________________
%
% spm_sections is called by spm_results_ui and uses variable img to
% create three orthogonal sections through a background image.
% Regional foci from the selected xSPM are rendered on this image.
%__________________________________________________________________________
% Copyright (C) 1994-2011 Wellcome Trust Centre for Neuroimaging

% John Ashburner
% $Id: spm_sections.m 4199 2011-02-10 20:07:17Z guillaume $


if ~nargin, [SPM,xSPM] = spm_getSPM; end
if nargin < 2, hReg = []; end
if nargin < 3 || isempty(img)
    [img, sts] = spm_select(1,'image','Select image for rendering on');
    if ~sts, return; end
end

Fgraph = spm_figure('GetWin','Graphics');
spm_results_ui('Clear',Fgraph);
spm_orthviews('Reset');

global st prevsect
st.Space = spm_matrix([0 0 0  0 0 -pi/2]) * st.Space;
prevsect = img;

h = spm_orthviews('Image', img, [0.05 0.05 0.9 0.45]);
spm_orthviews('AddContext', h); 
spm_orthviews('MaxBB');
if ~isempty(hReg), spm_orthviews('Register', hReg); end
spm_orthviews('AddBlobs', h, xSPM.XYZ, xSPM.Z, xSPM.M);
spm_orthviews('Redraw');

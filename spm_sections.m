function spm_sections(SPM,hReg,spms)
% Rendering of regional effects [SPM{Z}] on orthogonal sections
% FORMAT spm_sections(SPM,hReg)
%
% SPM  - xSPM structure containing details of excursion set
% hReg - handle of MIP register
%
% see spm_getSPM for details
%__________________________________________________________________________
%
% spm_sections is called by spm_results and uses variables in SPM and
% VOL to create three orthogonal sections though a background image.
% Regional foci from the selected SPM are rendered on this image.
%
%__________________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% John Ashburner
% $Id: spm_sections.m 3081 2009-04-22 20:15:38Z guillaume $


if nargin < 3 || isempty(spms)
    spms   = spm_select(1,'image','select image for rendering on');
end

Fgraph = spm_figure('FindWin','Graphics');
spm_results_ui('Clear',Fgraph);
spm_orthviews('Reset');
global st
st.Space = spm_matrix([0 0 0  0 0 -pi/2])*st.Space;
h = spm_orthviews('Image',spms,[0.05 0.05 0.9 0.45]);
spm_orthviews('AddContext',h); 
spm_orthviews MaxBB;
spm_orthviews('register',hReg);
spm_orthviews('addblobs',1,SPM.XYZ,SPM.Z,SPM.M);
spm_orthviews('Redraw');

global prevsect
prevsect = spms;


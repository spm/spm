function spm_image_funks(varargin)
% Perform algebraic functions on images.
% FORMAT spm_image_funks(P,Q,func) OR spm_wi
% P    - matrix of input image filenames.
% Q    - name of output image.
% func - the expression to be evaluated.
%
%_______________________________________________________________________
%
% The images specified in P, are referred to as i1, i2, i3...
% in the expression to be evaluated.
%
% With images of different sizes and orientations, the size and
% orientation of the first is used for the output image.
%
% The image Q is written to current working directory unless a valid
% full pathname is given.
%
%
%                           ----------------
%
% spm_image_funks is grandfathered - use spm_imcalc_ui instead
%_______________________________________________________________________
% %W% John Ashburner %E%

%-Print warning of obsolescence
%-----------------------------------------------------------------------
warning([mfilename,' is grandfathered, use spm_imcalc_ui instead'])


%-Pass on arguments to spm_imcalc_ui
%-----------------------------------------------------------------------
spm_imcalc_ui(varargin{:});

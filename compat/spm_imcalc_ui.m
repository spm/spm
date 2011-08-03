function [Q,Vo] = spm_imcalc_ui(P,Q,f,flags,varargin)
% Perform algebraic functions on images
% FORMAT Q = spm_imcalc_ui(P,Q,f,flags)
% P             - matrix of input image filenames
%                 [user prompted to select files if arg missing or empty]
% Q             - name of output image
%                 [user prompted to enter filename if arg missing or empty]
% f             - expression to be evaluated
%                 [user prompted to enter expression if arg missing or empty]
% flags         - cell vector of flags: {dmtx,mask,type,hold}
% dmtx          - Read images into data matrix?
%                 [defaults (missing or empty) to 0 - no]
% mask          - implicit zero mask?
%                 [defaults (missing or empty) to 0]
% type          - data type for output image (see spm_type)
%                 [defaults (missing or empty) to 4 - 16 bit signed shorts]
% hold          - interpolation hold (see spm_slice_vol)
%                 [defaults (missing or empty) to 0 - nearest neighbour]
%
% Q (output)    - full pathname of image written
% Vo            - structure containing information on output image (see spm_vol)
%__________________________________________________________________________
%
% This function is now deprecated, use spm_imcalc instead.
%__________________________________________________________________________
% Copyright (C) 1998-2011 Wellcome Trust Centre for Neuroimaging

% John Ashburner & Andrew Holmes
% $Id: spm_imcalc_ui.m 4418 2011-08-03 12:00:13Z guillaume $

persistent runonce
if isempty(runonce)
    warning('spm_imcalc_ui is deprecated. Use spm_imcalc instead.');
    runonce = 1;
end

if nargin < 3, spm_imcalc; return; end
if nargin < 4, flags = {}; end

%-GUI setup
%--------------------------------------------------------------------------
[Finter,Fgraph,CmdLine] = spm('FnUIsetup','ImCalc',0);

spm('FigName','ImCalc: working',Finter,CmdLine);
spm('Pointer','Watch');

flags = {flags{:} [] [] [] []};
flags([3 4]) = flags([4 3]);
Vo    = spm_imcalc(P, Q, f, flags(1:4), varargin{:});
Q     = Vo.fname;

%-End
%--------------------------------------------------------------------------
spm('Pointer');
spm('FigName','ImCalc: done',Finter,CmdLine);

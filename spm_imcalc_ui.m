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
% Q (output)    - full pathname of image written
% Vo            - structure containing information on output image (see spm_vol)
%
%_______________________________________________________________________
%
% spm_imcalc_ui uses spm_imcalc as an engine to perform user-specified
% algebraic manipulations on a set of images, with the result being
% written out as an image. The user is prompted to supply images to
% work on, a filename for the output image, and the expression to
% evaluate. The expression should be a standard matlab expression,
% within which the images should be referred to as i1, i2, i3,... etc.
%
% If the dmtx flag is set, then images are read into a data matrix X
% (rather than into seperate variables i1, i2, i3,...). The data matrix
% should be referred to as X, and contains images in rows.
%
% Computation is plane by plane, so in data-matrix mode, X is a NxK
% matrix, where N is the number of input images [prod(size(Vi))], and K
% is the number of voxels per plane [prod(Vi(1).dim(1:2))].
%
% For data types without a representation of NaN, implicit zero masking
% assummes that all zero voxels are to be treated as missing, and
% treats them as NaN. NaN's are written as zero (by spm_write_plane),
% for data types without a representation of NaN.
%
% With images of different sizes and orientations, the size and
% orientation of the first is used for the output image. A warning is
% given in this situation. Images are sampled into this orientation
% using the interpolation specified by the hold parameter.  [default -
% nearest neighbour]
%
% The image Q is written to current working directory unless a valid
% full pathname is given.
%
% Example expressions (f):
%
%        i) Mean of six images (select six images)
%           f = '(i1+i2+i3+i4+i5+i6)/6'
%       ii) Make a binary mask image at threshold of 100
%           f = 'i1>100'
%      iii) Make a mask from one image and apply to another
%           f = 'i2.*(i1>100)'
%                 - here the first image is used to make the mask, which is
%                   applied to the second image
%       iv) Sum of n images
%           f = 'i1 + i2 + i3 + i4 + i5 + ...'
%        v) Sum of n images (when reading data into a data-matrix - use dmtx arg)
%           f = 'sum(X)'
% 
% Parameters can be passed as arguments to override internal defaults
% (for hold, mask & type), or to pre-specify images (P), output
% filename (Q), or expression (f). Pass empty matrices for arguments
% not to be set.
% E.g.  Q = spm_imcalc_ui({},'test','',{[],[],[],1})
%       ... pre-specifies the output filename as 'test.img' in the current
% working directory, and sets the interpolation hold to tri-linear.
%
% Further, if calling spm_imcalc directly, additional variables for use in
% the computation can be passed at the end of the argument list. These
% should be referred to by the names of the arguments passed in the
% expression to be evaluated. E.g. if c is a 1xn vector of weights, then
% for n images, using the (dmtx) data-matrix version, the weighted sum can
% be computed using:
%       Vi= spm_vol(spm_select(inf,'image'));
%       Vo= Vi(1);
%       Vo.fname = 'output.img';
%       Vo.pinfo(1:2) = Inf;
%       Q = spm_imcalc(Vi,Vo,'c*X',{1},c)
% Here we've pre-specified the expression and passed the vector c as an
% additional variable (you'll be prompted to select the n images).
%__________________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% John Ashburner & Andrew Holmes
% $Id: spm_imcalc_ui.m 3691 2010-01-20 17:08:30Z guillaume $

%-GUI setup
%--------------------------------------------------------------------------
SVNid = '$Rev: 3691 $';
[Finter,Fgraph,CmdLine] = spm('FnUIsetup','ImCalc',0);
spm('FnBanner',mfilename,SVNid);

%-Condition arguments
%--------------------------------------------------------------------------
if nargin<4, flags={}; end
if nargin<3, f=''; end
if nargin<2, Q=''; end
if nargin<1, P={}; end

if isempty(P), P = spm_select(Inf,'image','Select images to work on'); end
if isempty(P), error('no input images specified'), end
if isempty(Q), Q = spm_input('Output filename',1,'s'); end
if isempty(f), f = spm_input('Evaluated Function',2,'s'); end

if length(flags)<4, hold=[]; else hold=flags{4}; end
if isempty(hold), hold=0; end
if length(flags)<3, type=[]; else type=flags{3}; end
if isempty(type), type=4; end, if ischar(type), type=spm_type(type); end
if length(flags)<2, mask=[]; else mask=flags{2}; end
if isempty(mask), mask=0; end
if length(flags)<1, dmtx=[]; else dmtx=flags{1}; end
if isempty(dmtx), dmtx=0; end

spm('FigName','ImCalc: working',Finter,CmdLine);
spm('Pointer','Watch')

%-Map input files
%--------------------------------------------------------------------------
Vi = spm_vol(char(P));
if isempty(Vi), error('no input images specified'), end

%-Check for consistency of image dimensions and orientation / voxel size
%--------------------------------------------------------------------------
if length(Vi)>1 && any(any(diff(cat(1,Vi.dim),1,1),1))
    warning(['images don''t all have same dimensions',...
        ' - using those of 1st image']);
end
if any(any(any(diff(cat(3,Vi.mat),1,3),3)))
    warning(['images don''t all have same orientation & voxel size',...
        ' - using 1st image']);
end

%-Work out filename for output image
%--------------------------------------------------------------------------
[p n e] = spm_fileparts(Q);
if isempty(p), p = pwd; end
if ~exist(p,'dir')
    warning('Invalid directory: writing to current directory')
    p = pwd;
end

Vo = struct('fname',   fullfile(p, [n, e]),...
            'dim',     Vi(1).dim(1:3),...
            'dt',      [type spm_platform('bigend')],...
            'mat',     Vi(1).mat,...
            'descrip', 'spm - algebra');

%-Call spm_imcalc to handle computations
%--------------------------------------------------------------------------
args = {dmtx,mask,hold};
Vo   = spm_imcalc(Vi,Vo,f,args);

%-End
%--------------------------------------------------------------------------
spm('Pointer');
spm('FigName','ImCalc: done',Finter,CmdLine);

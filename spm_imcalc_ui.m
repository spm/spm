function Q = spm_imcalc_ui(P,Q,f,hold,mask,type)
% Perform algebraic functions on images
% FORMAT spm_imcalc_ui(P,Q,f,hold,mask,type)
% P          - matrix of input image filenames
%              [user prompted to select files if arg missing or empty]
% Q          - name of output image
%              [user prompted to enter filename if arg missing or empty]
% f          - expression to be evaluated
%              [user prompted to enter expression if arg missing or empty]
% hold       - interpolation hold (see spm_slice_vol)
%              [defaults (missing or empty) to 0 - nearest neighbour]
% mask       - implicit zero mask?
%              [defaults (missing or empty) to 0]
% type       - data type for output image (see spm_type)
%              [defaults (missing or empty) to 4 - 16 bit signed shorts]
%              * currently only supports type 4! *
% Q (output) - full pathname of image written
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
% 
% Parameters can be passed as arguments to override internal defaults
% (for hold, mask & type), or to pre-specify images (P), output
% filename (Q), or expression (f). Pass empty matrices for arguments
% not to be set.
% E.g.	Q = spm_imcalc_ui({},'test','',1)
%       ... pre-specifies the output filename as 'test.img' in the current
% working directory, and sets the interpolation hold to tri-linear.
%
%_______________________________________________________________________
% %W% John Ashburner, Andrew Holmes %E%

%-GUI setup
%-----------------------------------------------------------------------
SCCSid = '%I%';
[Finter,Fgraph,CmdLine] = spm('FnUIsetup','ImCalc',0);
SPMid = spm('FnBanner',mfilename,SCCSid);
spm_help('!ContextHelp',[mfilename,'.m'])

%-Condition arguments
%-----------------------------------------------------------------------
if nargin<6, type=[]; end, if isempty(type), type=4; end
if ischar(type), type=spm_type(type); end
if ~any(type==[4]), error('invalid type'), end
if nargin<5, mask=[]; end, if isempty(mask), mask=0; end
if nargin<4, hold=[]; end, if isempty(hold), hold=0; end
if nargin<3, f=''; end
if nargin<2, Q=''; end
if nargin<1, P={}; end

if isempty(P), P = spm_get(Inf,'.img',{'Select images to work on'}); end
if isempty(P), error('no input images specified'), end
if isempty(Q), Q = spm_input('Output filename',1,'s'); end
if isempty(f), f = spm_input('Evaluated Function',2,'s'); end

spm('FigName','ImCalc: working',Finter,CmdLine);
spm('Pointer','Watch')


%-Map input files
%-----------------------------------------------------------------------
Vi = spm_vol(char(P));
if isempty(Vi), error('no input images specified'), end

%-Check for consistency of image dimensions and orientation / voxel size
%-----------------------------------------------------------------------
if length(Vi)>1 & any(any(diff(cat(1,Vi.dim),1,1),1)&[1,1,1,0])
	warning(['images don''t all have same dimensions',...
		' - using those of 1st image']), end
if any(any(any(diff(cat(3,Vi.mat),1,3),3)))
	warning(['images don''t all have same orientation & voxel size',...
		' - using 1st image']), end


%-Work out filename for output image
%------------------------------------------------------------------
Qdir = spm_str_manip(Q,'Hv');
Qfil = [spm_str_manip(Q,'stv'),'.img'];
if ~exist(Qdir,'dir')
	warning('Invalid directory: writing to current directory')
	Qdir = '.';
end
Q = spm_get('CPath',Qfil,Qdir);

Vo = struct(	'fname',	Q,...
		'dim',		[Vi(1).dim(1:3),type],...
		'mat',		Vi(1).mat,...
		'pinfo'	,	[1,0,0]',...
		'descrip',	'spm - algebra');


%-Call spm_imcalc to handle computations
%------------------------------------------------------------------
spm_create_image(Vo);
Vo.pinfo(1)  = spm_imcalc(Vi,Vo,f,hold,mask);
spm_create_image(Vo);

%-End
%------------------------------------------------------------------
spm('Pointer');
spm('FigName','ImCalc: done',Finter,CmdLine);

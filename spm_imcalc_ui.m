function Q = spm_imcalc_ui(P,Q,f)
% Perform algebraic functions on images.
% FORMAT spm_imcalc_ui(P,Q,f)
% P          - matrix of input image filenames
% Q          - name of output image
% f          - expression to be evaluated
% Q (output) - full pathname of image written
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
%_______________________________________________________________________
% %W% John Ashburner, Andrew Holmes %E%

%-Condition arguments
%-----------------------------------------------------------------------
SCCSid = '%I%';
if nargin<3
	[Finter,Fgraph,CmdLine] = spm('FnUIsetup','ImCalc',0);
	SPMid = spm('FnBanner',mfilename,SCCSid);
	spm_help('!ContextHelp',[mfilename,'.m'])
	if nargin<1, P = spm_get(Inf,'.img',{'Select images to work on'}); end
	if isempty(P), error('no input images specified'), end
	if nargin<2, Q = spm_input('Output filename',1,'s'); end
	if nargin<3, f = spm_input('Evaluated Function',2,'s'); end
	spm('FigName','ImCalc: working',Finter,CmdLine);
	spm('Pointer','Watch')
	Q = spm_imcalc_ui(P,Q,f);
	spm('Pointer');
	spm('FigName','ImCalc: done',Finter,CmdLine);
	return
end

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
		'dim',		[Vi(1).dim(1:3),4],...
		'mat',		Vi(1).mat,...
		'pinfo'	,	[1,0,0]',...
		'descrip',	'spm - algebra');


%-Call spm_imcalc to handle computations
%------------------------------------------------------------------
spm_create_image(Vo);
Vo.pinfo(1)  = spm_imcalc(Vi,Vo,f);
spm_create_image(Vo);

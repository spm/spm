
% user interface routine for spm_sn
% FORMAT spm_sn_ui
%____________________________________________________________________________
%
% If the data to be normalized are MRI data spm_sn_ui calls spm_mri_sn_ui,
% otherwise the data are assumed to be [S]PET like and:-
%
% spm_sn_ui prompts for a list of filenames of images that are passed
% to spm_sn - the spatial normalization routine.  The files specified
% for each subject will be averaged prior to determining the spatial
% transformation and should therefore be realigned if necessary.  Clearly
% this only applies to time-series (e.g. activation studies).  In cross-
% sectional studies there will be one image per subject.
%
% Other images that are to be subject to the same spatial transformation
% (but not used in the averaging) can be specified.  This allows co-registered
% structural MRI scans to be spatially transformed in parallel with the
% functional (e.g. PET) data, or (in radioligand studies) later frames
% (sensitive to Vd) to be normalized on the basis of flow-sensitive frames.
%
% spm_sn_ui calls spm_sn with 7 arguments: A row matrix (P) of all the
% images to be spatially normalized and a vector (q) indicating which
% images are to be used in determining the transformation (i.e. those to be
% entered into the averaging). Q is the template to which the images are
% matched.  bb and Vox define the bounding box and voxel size of the
% normalized image.  I is a 'best guess' pitch of the image[s] from the
% ac-pc line (usually about 8 degrees for cantho-meatal positioning.
% The final parameter (H) can be used to tailor the affine or the non-linear
% compoments of the normalization.  The latter components are (i) a fully 3-
% dimensional non-linear transformation using quadratic basis functions (or
% second order terms in the mapping of one coordinate system onto the other)
% and (ii) a slice-based (two dimensional) piece-wise non-linear transformation
% using a 'more general' (Fourier-like) set of basis functions.
%
% It is possible to change the last four of these arguments or to fix pitch: - 
%
% After specifying all the requisite filnames a panel of defaults will be
% displayed. These defaults can be changed to deal with 'special' requirements
% e.g:
%
% a) normalizing the cerebellum and ventral parts of the brain (bb)
% b) decreasing voxel size in the normalized images to 'match' an MRI scan
%    (Vox)
% c) changing the template to a SPECT or fMRI reference image (Q)
% d) imposing linearity (12 parameter affine) on the normalization (H)
% e) fixing the pitch for data with small numbers of slices (H)
%
% Linear or affine options:
%
% 6  - parameter: - A rigid body transformation
% 7  - parameter: - A size and position only transformation
% 9  - parameter: - A size and orthogonal scaling transformation
% 12 - parameter: - A full affine transformation
%
% Non-linear options:
%
% One or both non-linear components can be disabled.  It is necessary
% to diable both to effect a completely linear (e.g. rigid body)
% transformation.
%
% In general the better the data the more comprehensive the nonlinear
% transformations can be.   The slice based nonlinear component will
% produce artifacts for planes with partial data but these are not usually
% problematic in further analysis.  A potential problem with the bottom
% and top planes may be encountered if the slice-based option is used.
% inappropriate distortions in these circustances could confound the
% interpretation of single subject studies.  We recommend that the 2-D
% option be disabled for single subject studies, particularly those relying
% on precise structure-function relationships.
% 
% There are no restrictions on image size or voxel size.
%
% see also spm_sn.m
%
%__________________________________________________________________________
% %W% %E%

%----------------------------------------------------------------------------
global SWD


% cycle through subjects and images
%----------------------------------------------------------------------------
set(2,'Name','Spatial normalization'); figure(2); clf
if spm_input('Modality',1,'[S]PET|MRI',[0 1])
	spm_mri_sn_ui;
	return;
end
n     = spm_input('number of subjects',1);
for i = 1:n
	p = spm_input(['number of scans: subject ' num2str(i)],2);
	U = spm_get(p,'.img');
	q = spm_input(['number of other images'],3);
	V = spm_get(q,'.img');

	if length(V); U = str2mat(U,V); end
	q = [ones(1,p) zeros(1,q)];
	eval(['P' num2str(i) ' = U; q' num2str(i) ' = q;']);
end


% set spatial paramters of standard space in mm (usually the below defaults)
%----------------------------------------------------------------------------
spms   = [SWD '/spm.img'  ];				% model image string
bb     = [ [-64 64]' [-104 68]' [-28 72]' ];		% bb in Talairach space
Vox    = [2 2 4];					% voxel size {x y and z}
I      = 8;						% pitch estimate
H      = [0 1 1 12];					% Default option


% scale positions according to size of figure 2
%---------------------------------------------------------------------------
set(2,'Units','Pixels');
A    = get(2,'Position');
A    = diag([A(3)/400 A(4)/395 A(3)/400 A(4)/395 ]);

% text frames and labels
%----------------------------------------------------------------------------
uicontrol(2,'style','frame','Position',[10 010 380 290]*A);
uicontrol(2,'style','text', 'Position',[20 274 360 20]*A,'string',...
	'Bounding box {x y z} and voxel size',...
	'ForegroundColor',[1 0 0]);
uicontrol(2,'style','text', 'Position',[20 170 360 20]*A,'string',...
	'Template or "model" image',...
	'ForegroundColor',[1 0 0]);
uicontrol(2,'style','text', 'Position',[20 110 160 20]*A,'string',...
	'pitch estimate {degrees}',...
	'ForegroundColor',[1 0 0]);
uicontrol(2,'style','text', 'Position',[240 110 60 20]*A,'string',...
	' - fix ?',...
	'ForegroundColor',[1 0 0]);
uicontrol(2,'style','text', 'Position',[20 080 165 20]*A,'string',...
	'3-D linear - parameters',...
	'ForegroundColor',[1 0 0]);
uicontrol(2,'style','text', 'Position',[20 050 165 20]*A,'string',...
	'3-D nonlinear quadtratic',...
	'ForegroundColor',[1 0 0]);
uicontrol(2,'style','text', 'Position',[20 020 170 20]*A,'string',...
	'2-D nonlinear  -  Fourier',...
	'ForegroundColor',[1 0 0]);

uicontrol(2,'style','text','Position',[020 250 120 20]*A,'string',...
	'minimum {mm}');
uicontrol(2,'style','text','Position',[020 230 120 20]*A,'string',...
	'maximum {mm}');
uicontrol(2,'style','text','Position',[020 200 120 20]*A,'string',...
	'size {mm}   ');

% Bounding box defaults
%----------------------------------------------------------------------------
c   = ['bb(1,1) = eval(get(get(2,''CurrentObject''),''string''));']; 				uicontrol(2,'style','edit','Position',[160 250 40 16]*A,...
'string',sprintf('%3.0i',bb(1,1)),'callback',c);
c   = ['bb(1,2) = eval(get(get(2,''CurrentObject''),''string''));']; 				uicontrol(2,'style','edit','Position',[220 250 40 16]*A,...
'string',sprintf('%3.0i',bb(1,2)),'callback',c);
c   = ['bb(1,3) = eval(get(get(2,''CurrentObject''),''string''));']; 				uicontrol(2,'style','edit','Position',[280 250 40 16]*A,...
'string',sprintf('%3.0i',bb(1,3)),'callback',c);
c   = ['bb(2,1) = eval(get(get(2,''CurrentObject''),''string''));']; 				uicontrol(2,'style','edit','Position',[160 230 40 16]*A,...
'string',sprintf('%3.0i',bb(2,1)),'callback',c);
c   = ['bb(2,2) = eval(get(get(2,''CurrentObject''),''string''));']; 				uicontrol(2,'style','edit','Position',[220 230 40 16]*A,...
'string',sprintf('%3.0i',bb(2,2)),'callback',c);
c   = ['bb(2,3) = eval(get(get(2,''CurrentObject''),''string''));']; 				uicontrol(2,'style','edit','Position',[280 230 40 16]*A,...
'string',sprintf('%3.0i',bb(2,3)),'callback',c);
c   = ['Vox(1) = eval(get(get(2,''CurrentObject''),''string''));']; 				uicontrol(2,'style','edit','Position',[160 200 40 16]*A,...
'string',sprintf('%3.0i',Vox(1)),'callback',c);
c   = ['Vox(2) = eval(get(get(2,''CurrentObject''),''string''));']; 				uicontrol(2,'style','edit','Position',[220 200 40 16]*A,...
'string',sprintf('%3.0i',Vox(2)),'callback',c);
c   = ['Vox(3) = eval(get(get(2,''CurrentObject''),''string''));']; 				uicontrol(2,'style','edit','Position',[280 200 40 16]*A,...
'string',sprintf('%3.0i',Vox(3)),'callback',c);

% template default
%----------------------------------------------------------------------------
c        = ['spms = get(get(2,''CurrentObject''),''string'');']; 				uicontrol(2,'style','edit','Position',[40 140 320 20]*A,...
	    'string',spms,'callback',c);

% initial [best guess] pitch estimate
%----------------------------------------------------------------------------
c        = ['I = eval(get(get(2,''CurrentObject''),''string''));']; 				uicontrol(2,'style','edit','Position',[190 110 40 20]*A,...
	    'string',sprintf('%3.0i',I),'callback',c);
c        = ['H(1) = get(get(2,''CurrentObject''),''Value'');'];
uicontrol(2,'style','Radiobutton','Position',...
	   [300 110 20 20]*A,'callback',c);

% number of affine parameters
%----------------------------------------------------------------------------
c        = ['H(4) = eval(get(get(2,''CurrentObject''),''string''));']; 				uicontrol(2,'style','edit','Position',[190 080 40 20]*A,...
	    'string',sprintf('%3.0i',H(4)),'callback',c);

% disabling nonlinear components
%----------------------------------------------------------------------------
c        = ['H(2) = get(get(2,''CurrentObject''),''Value'');'];
uicontrol(2,'style','Radiobutton','Position',...
	   [200 48 20 20]*A,'callback',c,'Value',1);
c        = ['H(3) = get(get(2,''CurrentObject''),''Value'');'];
uicontrol(2,'style','Radiobutton','Position',...
	   [200 18 20 20]*A,'callback',c,'Value',1);

% implement spatial normalization
%----------------------------------------------------------------------------
c   =  ['set(2,''Name'',''executing'',''Pointer'',''Watch''); drawnow;',...
	'for i = 1:n;',...
		'eval([''P = P'' num2str(i) ''; q = q'' num2str(i) '';'']);',...
		'spm_sn(P,q,spms,bb,Vox,H,I); end; ',...
	'figure(2); clf; set(2,''Name'','''',''Pointer'',''Arrow'')'];
uicontrol(2,'style','PushButton','Position',[300 030 60 30]*A,...
	'string','OK ?','callback',c);



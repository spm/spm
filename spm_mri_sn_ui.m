
% user interface routine for spm_mri_sn
% FORMAT spm_mri_sn_ui
%____________________________________________________________________________
%
% spm_mri_sn_ui prompts for a list of filenames of images that are passed
% to spm_mri_sn - the spatial normalization routine.
%
% Other images that are to be subject to the same spatial transformation
% (but not used in the solutions) can be specified.  This allows co-registered
% functional (e.g. fMRI) scans to be spatially transformed in parallel with the
% structural data,
%
% spm_mri_sn_ui calls spm_mri_sn with 6 arguments: A row matrix (P) of the
% images to be spatially normalized.  The first image is taken to be the MRI
% used in determining the transformation (This image can be s T1, T2 or T2*
% weighted scan). M is a similar matrix specifying the template to which the
% images are matched.  The first should be a white matter template and the
% second a gray matter template.  Both templates should by congruent.
% bb and Vox define the bounding box and voxel size of the
% normalized image.  I is a 'best guess' pitch, roll and yaw relative to the
% ac-pc plane (usually about 8 degrees pitch for cantho-meatal positioning.
% The final parameter (H) can be used to fix the pitch, roll or yaw to these
% specified estimates.
% 
% There are no restrictions on image size or voxel size.
%
% The 'other' images can have different image and voxel dimensions as long
% as the voxel specified by the ORIGIN in the header correponds to the
% equivalent voxel in the MRI image.
%
% see also spm_mri_sn.m
%
%__________________________________________________________________________
% %W% %E%


% cycle through subjects and images
%----------------------------------------------------------------------------
set(2,'Name','Spatial normalization'); figure(2); clf
n     = spm_input('number of subjects',1);
for i = 1:n
	str = ['Subject ' int2str(i) ':'];
	U   = spm_get(1,'.img',[str ' Select T1/2[*]-weighted MRI']);
	q   = spm_input(['number of other images'],2);
	V   = spm_get(q,'.img',[str ' Select other (e.g. fMRI or PET) images']);

	if length(V); U = str2mat(U,V); end
	eval(['P' num2str(i) ' = U;']);
end

% set spatial paramters of standard space in mm (usually the below defaults)
%----------------------------------------------------------------------------
global SWD
mriW   = [SWD '/mriWs.img'  ];				% white template
mriG   = [SWD '/mriGs.img'  ];				% gray template
bb     = [ [-64 64]' [-104 68]' [-28 72]' ];		% bb in Talairach space
Vox    = [2 2 4];					% voxel size {x y and z}
I      = [0 0 0];					% pitch estimate
H      = [0 0 0];					% MRI options


% scale positions according to size of figure 2
%---------------------------------------------------------------------------
set(2,'Units','Pixels');
A    = get(2,'Position');
A    = diag([A(3)/400 A(4)/395 A(3)/400 A(4)/395 ]);

% text frames and labels
%----------------------------------------------------------------------------
uicontrol(2,'style','frame','Position',[10 010 380 290]*A);
uicontrol(2,'style','text', 'Position',[20 274 360 20]*A,'string',...
	'Bounding box {x y z}','ForegroundColor',[1 0 0]);
uicontrol(2,'style','text', 'Position',[20 200 360 20]*A,'string',...
	'New voxel size {x y z}','ForegroundColor',[1 0 0]);
uicontrol(2,'style','text', 'Position',[20 080 120 20]*A,'string',...
	'pitch {degrees}','ForegroundColor',[1 0 0]);
uicontrol(2,'style','text', 'Position',[20 050 120 20]*A,'string',...
	'roll  {degrees}','ForegroundColor',[1 0 0]);
uicontrol(2,'style','text', 'Position',[20 020 120 20]*A,'string',...
	'yaw   {degrees}','ForegroundColor',[1 0 0]);
uicontrol(2,'style','text', 'Position',[220 080 40 20]*A,'string',...
	'fix?','ForegroundColor',[1 0 0]);
uicontrol(2,'style','text', 'Position',[220 050 40 20]*A,'string',...
	'fix?','ForegroundColor',[1 0 0]);
uicontrol(2,'style','text', 'Position',[220 020 40 20]*A,'string',...
	'fix?','ForegroundColor',[1 0 0]);

uicontrol(2,'style','text', 'Position',[20 136 80 20]*A,'string',...
	'white and');
uicontrol(2,'style','text', 'Position',[20 110 80 20]*A,'string',...
	'gray model');

uicontrol(2,'style','text','Position',[020 250 120 20]*A,'string',...
	'minimum {mm}');
uicontrol(2,'style','text','Position',[020 230 120 20]*A,'string',...
	'maximum {mm}');
uicontrol(2,'style','text','Position',[020 170 120 20]*A,'string',...
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
c   = ['Vox(1) = eval(get(get(2,''CurrentObject''),''string''));']; 				uicontrol(2,'style','edit','Position',[160 170 40 16]*A,...
'string',sprintf('%3.2f',Vox(1)),'callback',c);
c   = ['Vox(2) = eval(get(get(2,''CurrentObject''),''string''));']; 				uicontrol(2,'style','edit','Position',[220 170 40 16]*A,...
'string',sprintf('%3.2f',Vox(2)),'callback',c);
c   = ['Vox(3) = eval(get(get(2,''CurrentObject''),''string''));']; 				uicontrol(2,'style','edit','Position',[280 170 40 16]*A,...
'string',sprintf('%3.2f',Vox(3)),'callback',c);

% template default
%----------------------------------------------------------------------------
c        = ['mriW = get(get(2,''CurrentObject''),''string'');']; 				uicontrol(2,'style','edit','Position',[100 136 280 20]*A,...
	    'string',mriW,'callback',c);
c        = ['mriG = get(get(2,''CurrentObject''),''string'');']; 				uicontrol(2,'style','edit','Position',[100 110 280 20]*A,...
	    'string',mriG,'callback',c);

% initial [best guess] pitch estimate
%----------------------------------------------------------------------------
c        = ['I(1) = eval(get(get(2,''CurrentObject''),''string''));']; 				uicontrol(2,'style','edit','Position',[160 80 40 20]*A,...
	    'string',sprintf('%3.1f',I(1)),'callback',c);
c        = ['I(2) = eval(get(get(2,''CurrentObject''),''string''));']; 				uicontrol(2,'style','edit','Position',[160 50 40 20]*A,...
	    'string',sprintf('%3.1f',I(2)),'callback',c);
c        = ['I(3) = eval(get(get(2,''CurrentObject''),''string''));']; 				uicontrol(2,'style','edit','Position',[160 20 40 20]*A,...
	    'string',sprintf('%3.1f',I(3)),'callback',c);
c        = ['H(1) = get(get(2,''CurrentObject''),''Value'');'];
uicontrol(2,'style','Radiobutton','Position',...
	   [260 80 20 20]*A,'callback',c);
c        = ['H(2) = get(get(2,''CurrentObject''),''Value'');'];
uicontrol(2,'style','Radiobutton','Position',...
	   [260 50 20 20]*A,'callback',c);
c        = ['H(3) = get(get(2,''CurrentObject''),''Value'');'];
uicontrol(2,'style','Radiobutton','Position',...
	   [260 20 20 20]*A,'callback',c);

% implement spatial normalization
%----------------------------------------------------------------------------
c   =  ['set(2,''Name'',''executing'',''Pointer'',''Watch''); drawnow;',...
	'for i = 1:n;',...
		'eval([''P = P'' num2str(i) '';'']);',...
		'spm_mri_sn(P,str2mat(mriW,mriG),bb,Vox,H,I); end; ',...
	'figure(2); clf; set(2,''Name'','''',''Pointer'',''Arrow'')'];
uicontrol(2,'style','PushButton','Position',[310 020 60 30]*A,...
	'string','OK ?','callback',c);



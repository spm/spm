
% surface rendering of regional effects
% FORMAT spm_surface_ui
%___________________________________________________________________________
%
% spm_surface_ui renders a specified SPM{Z} (in image format) beneath
% a surface rendered reference image (usually a structural MRI volume)
% spm_surface_ui gets the filenames, thresholds to be applied to the SPM{Z}
% [and MRI] and the 'view' required
%
% The view is specified in terms of pitch roll and yaw transformations.
% The units of rotation are entered in RADIANS.
%
%__________________________________________________________________________
% %W% %E%


%---------------------------------------------------------------------------
global CWD SWD

% get SPM{Z}, reference image and threshold
%---------------------------------------------------------------------------
eval(['cd ' CWD]);
Z    = spm_get(1,'.img','select SPM{Z} for rendering','SPM');
CWD  = pwd;

eval(['cd ' SWD]);
Q    = spm_get(1,'.img','select image for cortical rendering',[]);

% default thresholds and transformation matrix
%---------------------------------------------------------------------------
set(2,'Pointer','watch')
figure(3); spm_clf;
B    = [0 0 0 0 0 0];
U    = 2.33;
V    = 64;

% create control objects
%---------------------------------------------------------------------------
uicontrol(3,'Style','Frame','Position',[120 20 380 140]);
uicontrol(3,'Style','Text', 'Position',[140 120 040 16],'String','pitch');
uicontrol(3,'Style','Text', 'Position',[140 090 040 16],'String','roll');
uicontrol(3,'Style','Text', 'Position',[140 060 040 16],'String','yaw');
uicontrol(3,'Style','Text', 'Position',[260 120 120 16],'String','SPM{Z} threshold');
uicontrol(3,'Style','Text', 'Position',[250 090 120 16],'String','MRI threshold');


c    = ['B(6) =  eval(get(get(3,''CurrentObject''),''String''));'];
uicontrol(3,'Style','Edit', 'Position',[200 120 40 16],'CallBack',c);
c    = ['B(5) =  eval(get(get(3,''CurrentObject''),''String''));'];
uicontrol(3,'Style','Edit', 'Position',[200 090 40 16],'CallBack',c);
c    = ['B(4) =  eval(get(get(3,''CurrentObject''),''String''));'];
uicontrol(3,'Style','Edit', 'Position',[200 060 40 16],'CallBack',c);

c    = ['U = eval(get(get(3,''CurrentObject''),''String''));'];
uicontrol(3,'Style','Edit', 'Position',[430 120 40 16],'CallBack',c,...
	'String',num2str(U));
c    = ['V = get(get(3,''CurrentObject''),''Value'');'];
uicontrol(3,'Style','Slider', 'Position',[370 090 100 16],'CallBack',c,...
	'Min',0,'Max',255,'Value',V);

c    = ['set(3,''Pointer'',''Watch'');',...
	'spm_surface(Z,Q,U,V,spm_matrix(B));',...
	'set(3,''Pointer'',''Arrow'')'];
uicontrol(3,'Style','Pushbutton', 'Position',[320 40 120 30],'CallBack',c,...
	'String','render','ForegroundColor',[1 0 0]);

% default rendering
%---------------------------------------------------------------------------
spm_surface(Z,Q,U,V,spm_matrix(B));
set(2,'Pointer','arrow')


















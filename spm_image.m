
% image display and header edit
% FORMAT spm_image
%____________________________________________________________________________
%
% spm_image is an interactive facility that allows you to edit or create
% image headers, display images, translate, rotate and scale images and
% write the resulting image to the disk.  Transformed images are written
% to the same subdirectory as the original image with the same name but
% prefixed with a 't'
%
% This facility also allows the creation and amendment of header files
% by providing access to spm_fix_header (a standalone module).
%
% The red lines drawn on the images are the planes that contain the ORIGIN as
% specified in the header.  To change this, or any other parameter, save
% the modified header (and display the image again to check any changes)
%
%__________________________________________________________________________
% %W% %E%


S    = get(0,'ScreenSize');
A    = diag([S(3)/1152 S(4)/900 S(3)/1152 S(4)/900]);

% get the image's filename {P} and set output filename {Q}
%----------------------------------------------------------------------------
P      = spm_get(1,'.img','please select image',[]);
P      = P(P ~= ' ');
d      = max([find(P == '/') 0]);
Q      = [P(1:d) 't' P((d + 1):length(P))];
B      = [0 0 0 0 0 0 1 1 1];				% tranformation matrix

[DIM VOX SCALE TYPE OFFSET ORIGIN DESCRIP] = spm_hread(P);
figure(3); spm_clf


% frames and text
%---------------------------------------------------------------------------
uicontrol(3,'Style','Frame','Position',[60 030 200 320]*A);
uicontrol(3,'Style','Frame','Position',[300 30 280 320]*A);
uicontrol(3,'Style','Text', 'Position',[90 260 100 16]*A,'String','right  {mm}');
uicontrol(3,'Style','Text', 'Position',[90 240 100 16]*A,'String','foward  {mm}');
uicontrol(3,'Style','Text', 'Position',[90 220 100 16]*A,'String','up  {mm}');
uicontrol(3,'Style','Text', 'Position',[90 200 100 16]*A,'String','pitch  {rad}');
uicontrol(3,'Style','Text', 'Position',[90 180 100 16]*A,'String','roll  {rad}');
uicontrol(3,'Style','Text', 'Position',[90 160 100 16]*A,'String','yaw  {rad}');
uicontrol(3,'Style','Text', 'Position',[90 140 100 16]*A,'String','resize  {x}');
uicontrol(3,'Style','Text', 'Position',[90 120 100 16]*A,'String','resize  {y}');
uicontrol(3,'Style','Text', 'Position',[90 100 100 16]*A,'String','resize  {z}');

uicontrol(3,'Style','Text','Position',[330 260 100 16]*A,'String','image {x}');
uicontrol(3,'Style','Text','Position',[330 240 100 16]*A,'String','image {y}');
uicontrol(3,'Style','Text','Position',[330 220 100 16]*A,'String','image {z}');
uicontrol(3,'Style','Text','Position',[330 200 100 16]*A,'String','voxel {x}');
uicontrol(3,'Style','Text','Position',[330 180 100 16]*A,'String','voxel {y}');
uicontrol(3,'Style','Text','Position',[330 160 100 16]*A,'String','voxel {z}');
uicontrol(3,'Style','Text','Position',[330 140 100 16]*A,'String','scaling');
uicontrol(3,'Style','Text','Position',[330 120 100 16]*A,'String','data type');
uicontrol(3,'Style','Text','Position',[330 100 100 16]*A,'String','offset');
uicontrol(3,'Style','Text','Position',[330 080 100 16]*A,'String','origin');
uicontrol(3,'Style','Text','Position',[330 060 100 16]*A,'String','description');

% objects with callback
%---------------------------------------------------------------------------
c    = ['set(3,''Pointer'',''watch'');',...
	'for i = 1:4; delete(gca); end;',...
	'spm_display(P,spm_matrix(B));',...
	'set(3,''Pointer'',''arrow'')'];
uicontrol(3,'Style','Pushbutton','String','display','Callback',c,...
         'Position',[110 300 100 20]*A);
c    = ['set(3,''Pointer'',''watch'');',...
	'spm_write(P,Q,spm_matrix(B));',...
	'set(3,''Pointer'',''arrow'');'];
uicontrol(3,'Style','Pushbutton','String','write image','Callback',c,...
         'Position',[110 60 100 20]*A);
c    = ['spm_hwrite(P,DIM,VOX,SCALE,TYPE,OFFSET,ORIGIN,DESCRIP);'];
uicontrol(3,'Style','Pushbutton','String','save header','Callback',c,...
         'Position',[320 300 110 20]*A);
c    = ['spm_fix_header;'];
uicontrol(3,'Style','Pushbutton','String','Fix header[s]','Callback',c,...
         'Position',[450 300 110 20]*A);


c    = ['B(1) = eval(get(get(3,''CurrentObject''),''string''));'];
uicontrol(3,'Style','edit','Callback',c,'Position',[190 260 40 16]*A);
c    = ['B(2) = eval(get(get(3,''CurrentObject''),''string''));'];
uicontrol(3,'Style','edit','Callback',c,'Position',[190 240 40 16]*A);
c    = ['B(3) = eval(get(get(3,''CurrentObject''),''string''));'];
uicontrol(3,'Style','edit','Callback',c,'Position',[190 220 40 16]*A);
c    = ['B(4) = eval(get(get(3,''CurrentObject''),''string''));'];
uicontrol(3,'Style','edit','Callback',c,'Position',[190 200 40 16]*A);
c    = ['B(5) = eval(get(get(3,''CurrentObject''),''string''));'];
uicontrol(3,'Style','edit','Callback',c,'Position',[190 180 40 16]*A);
c    = ['B(6) = eval(get(get(3,''CurrentObject''),''string''));'];
uicontrol(3,'Style','edit','Callback',c,'Position',[190 160 40 16]*A);
c    = ['B(7) = eval(get(get(3,''CurrentObject''),''string''));'];
uicontrol(3,'Style','edit','Callback',c,'Position',[190 140 40 16]*A);
c    = ['B(8) = eval(get(get(3,''CurrentObject''),''string''));'];
uicontrol(3,'Style','edit','Callback',c,'Position',[190 120 40 16]*A);
c    = ['B(9) = eval(get(get(3,''CurrentObject''),''string''));'];
uicontrol(3,'Style','edit','Callback',c,'Position',[190 100 40 16]*A);



c    = ['DIM(1) = eval(get(get(3,''CurrentObject''),''string''));'];
uicontrol(3,'Style','edit','Callback',c,'Position',[440 260 40 16]*A,...
	'String',sprintf('%0.0f',DIM(1)));
c    = ['DIM(2) = eval(get(get(3,''CurrentObject''),''string''));'];
uicontrol(3,'Style','edit','Callback',c,'Position',[440 240 40 16]*A,...
	'String',sprintf('%0.0f',DIM(2)));
c    = ['DIM(3) = eval(get(get(3,''CurrentObject''),''string''));'];
uicontrol(3,'Style','edit','Callback',c,'Position',[440 220 40 16]*A,...
	'String',sprintf('%0.0f',DIM(3)));
c    = ['VOX(1) = eval(get(get(3,''CurrentObject''),''string''));'];
uicontrol(3,'Style','edit','Callback',c,'Position',[440 200 40 16]*A,...
	'String',sprintf('%0.2f',VOX(1)));
c    = ['VOX(2) = eval(get(get(3,''CurrentObject''),''string''));'];
uicontrol(3,'Style','edit','Callback',c,'Position',[440 180 40 16]*A,...
	'String',sprintf('%0.2f',VOX(2)));
c    = ['VOX(3) = eval(get(get(3,''CurrentObject''),''string''));'];
uicontrol(3,'Style','edit','Callback',c,'Position',[440 160 40 16]*A,...
	'String',sprintf('%0.2f',VOX(3)));
c    = ['SCALE = eval(get(get(3,''CurrentObject''),''string''));'];
uicontrol(3,'Style','edit','Callback',c,'Position',[440 140 120 16]*A,...
	'String',sprintf('%0.4g',SCALE));
c    = ['TYPE = eval(get(get(3,''CurrentObject''),''string''));'];
uicontrol(3,'Style','edit','Callback',c,'Position',[440 120 40 16]*A,...
	'String',sprintf('%0.0f',TYPE));
c    = ['OFFSET = eval(get(get(3,''CurrentObject''),''string''));'];
uicontrol(3,'Style','edit','Callback',c,'Position',[440 100 40 16]*A,...
	'String',sprintf('%0.0f',OFFSET));
c    = ['ORIGIN(1) = eval(get(get(3,''CurrentObject''),''string''));'];
uicontrol(3,'Style','edit','Callback',c,'Position',[440 80 40 16]*A,...
	'String',sprintf('%i',ORIGIN(1)));
c    = ['ORIGIN(2) = eval(get(get(3,''CurrentObject''),''string''));'];
uicontrol(3,'Style','edit','Callback',c,'Position',[480 80 40 16]*A,...
	'String',sprintf('%i',ORIGIN(2)));
c    = ['ORIGIN(3) = eval(get(get(3,''CurrentObject''),''string''));'];
uicontrol(3,'Style','edit','Callback',c,'Position',[520 80 40 16]*A,...
	'String',sprintf('%i',ORIGIN(3)));
c    = ['DESCRIP = get(get(3,''CurrentObject''),''string'');'];
uicontrol(3,'Style','edit','Callback',c,'Position',[440 60 120 16]*A,...
	'String',DESCRIP);







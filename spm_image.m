
% image display and header edit
% FORMAT spm_image
%_______________________________________________________________________
%
% spm_image is an interactive facility that allows you to edit or create
% image headers, display images, translate, rotate and scale images and
% write the resulting image to the disk.  Transformed images are written
% to the same subdirectory as the original image with the same name but
% prefixed with a 't'
%
% This facility also allows the creation and amendment of header files
% by providing access to spm_header_edit (a standalone module).
%
% The red lines drawn on the images are the planes that contain the ORIGIN as
% specified in the header.  To change this, or any other parameter, save
% the modified header (and display the image again to check any changes)
%
%_______________________________________________________________________
% %W% Karl Friston %E%


% get the image's filename {P} and set output filename {Q}
%-----------------------------------------------------------------------
P      = spm_get(1,'.img','please select image',[]);
P      = P(P ~= ' ');
d      = max([find(P == '/') 0]);
Q      = [P(1:d) 't' P((d + 1):length(P))];
B      = [0 0 0 0 0 0 1 1 1];			% tranformation matrix

[DIM VOX SCALE TYPE OFFSET ORIGIN DESCRIP] = spm_hread(P);


% locate Graphics window and clear it
%-----------------------------------------------------------------------
F=spm_figure('FindWin','Graphics');
if isempty(F), F=spm_figure('Create','Graphics'); end
spm_clf(F)
figure(F)
WS = spm('GetWinScale');

% frames and text
%-----------------------------------------------------------------------
uicontrol(F,'Style','Frame','Position',[60 030 200 320].*WS);
uicontrol(F,'Style','Frame','Position',[300 30 280 320].*WS);
uicontrol(F,'Style','Text', 'Position',[90 260 100 016].*WS,...
	'String','right  {mm}');
uicontrol(F,'Style','Text', 'Position',[90 240 100 016].*WS,...
	'String','foward  {mm}');
uicontrol(F,'Style','Text', 'Position',[90 220 100 016].*WS,...
	'String','up  {mm}');
uicontrol(F,'Style','Text', 'Position',[90 200 100 016].*WS,...
	'String','pitch  {rad}');
uicontrol(F,'Style','Text', 'Position',[90 180 100 016].*WS,...
	'String','roll  {rad}');
uicontrol(F,'Style','Text', 'Position',[90 160 100 016].*WS,...
	'String','yaw  {rad}');
uicontrol(F,'Style','Text', 'Position',[90 140 100 016].*WS,...
	'String','resize  {x}');
uicontrol(F,'Style','Text', 'Position',[90 120 100 016].*WS,...
	'String','resize  {y}');
uicontrol(F,'Style','Text', 'Position',[90 100 100 016].*WS,...
	'String','resize  {z}');

uicontrol(F,'Style','Text','Position',[330 260 100 016].*WS,...
	'String','image {x}');
uicontrol(F,'Style','Text','Position',[330 240 100 016].*WS,...
	'String','image {y}');
uicontrol(F,'Style','Text','Position',[330 220 100 016].*WS,...
	'String','image {z}');
uicontrol(F,'Style','Text','Position',[330 200 100 016].*WS,...
	'String','voxel {x}');
uicontrol(F,'Style','Text','Position',[330 180 100 016].*WS,...
	'String','voxel {y}');
uicontrol(F,'Style','Text','Position',[330 160 100 016].*WS,...
	'String','voxel {z}');
uicontrol(F,'Style','Text','Position',[330 140 100 016].*WS,...
	'String','scaling');
uicontrol(F,'Style','Text','Position',[330 120 100 016].*WS,...
	'String','data type');
uicontrol(F,'Style','Text','Position',[330 100 100 016].*WS,...
	'String','offset');
uicontrol(F,'Style','Text','Position',[330 080 100 016].*WS,...
	'String','origin');
uicontrol(F,'Style','Text','Position',[330 060 100 016].*WS,...
	'String','description');

% objects with callback
%-----------------------------------------------------------------------
c    = ['set(gcf,''Pointer'',''watch'');',...
	'for i = 1:4; delete(gca); end;',...
	'spm_display(P,spm_matrix(B));',...
	'set(gcf,''Pointer'',''arrow'')'];
uicontrol(F,'Style','Pushbutton','String','display','Callback',c,...
         'Position',[110 300 100 020].*WS);
         
c    = ['set(gcf,''Pointer'',''watch'');',...
	'spm_write(P,Q,spm_matrix(B));',...
	'set(gcf,''Pointer'',''arrow'');'];
uicontrol(F,'Style','Pushbutton','String','write image','Callback',c,...
         'Position',[110 60 100 020].*WS);
         
c    = ['spm_hwrite(P,DIM,VOX,SCALE,TYPE,OFFSET,ORIGIN,DESCRIP);'];
uicontrol(F,'Style','Pushbutton','String','save header','Callback',c,...
         'Position',[320 300 110 020].*WS);
         
c    = ['spm_header_edit;'];
uicontrol(F,'Style','Pushbutton','String','Fix header[s]','Callback',c,...
         'Position',[450 300 110 020].*WS,'Interruptible','yes');


c    = ['B(1) = eval(get(gco,''String''));'];
uicontrol(F,'Style','edit','Callback',c,'Position',[190 260 040 016].*WS);
c    = ['B(2) = eval(get(gco,''String''));'];
uicontrol(F,'Style','edit','Callback',c,'Position',[190 240 040 016].*WS);
c    = ['B(3) = eval(get(gco,''String''));'];
uicontrol(F,'Style','edit','Callback',c,'Position',[190 220 040 016].*WS);
c    = ['B(4) = eval(get(gco,''String''));'];
uicontrol(F,'Style','edit','Callback',c,'Position',[190 200 040 016].*WS);
c    = ['B(5) = eval(get(gco,''String''));'];
uicontrol(F,'Style','edit','Callback',c,'Position',[190 180 040 016].*WS);
c    = ['B(6) = eval(get(gco,''String''));'];
uicontrol(F,'Style','edit','Callback',c,'Position',[190 160 040 016].*WS);
c    = ['B(7) = eval(get(gco,''String''));'];
uicontrol(F,'Style','edit','Callback',c,'Position',[190 140 040 016].*WS);
c    = ['B(8) = eval(get(gco,''String''));'];
uicontrol(F,'Style','edit','Callback',c,'Position',[190 120 040 016].*WS);
c    = ['B(9) = eval(get(gco,''String''));'];
uicontrol(F,'Style','edit','Callback',c,'Position',[190 100 040 016].*WS);



c    = ['DIM(1) = eval(get(gco,''String''));'];
uicontrol(F,'Style','edit','Callback',c,'Position',[440 260 040 016].*WS,...
	'String',sprintf('%0.0f',DIM(1)));
c    = ['DIM(2) = eval(get(gco,''String''));'];
uicontrol(F,'Style','edit','Callback',c,'Position',[440 240 040 016].*WS,...
	'String',sprintf('%0.0f',DIM(2)));
c    = ['DIM(3) = eval(get(gco,''String''));'];
uicontrol(F,'Style','edit','Callback',c,'Position',[440 220 040 016].*WS,...
	'String',sprintf('%0.0f',DIM(3)));
c    = ['VOX(1) = eval(get(gco,''String''));'];
uicontrol(F,'Style','edit','Callback',c,'Position',[440 200 040 016].*WS,...
	'String',sprintf('%0.2f',VOX(1)));
c    = ['VOX(2) = eval(get(gco,''String''));'];
uicontrol(F,'Style','edit','Callback',c,'Position',[440 180 040 016].*WS,...
	'String',sprintf('%0.2f',VOX(2)));
c    = ['VOX(3) = eval(get(gco,''String''));'];
uicontrol(F,'Style','edit','Callback',c,'Position',[440 160 040 016].*WS,...
	'String',sprintf('%0.2f',VOX(3)));
c    = ['SCALE = eval(get(gco,''String''));'];
uicontrol(F,'Style','edit','Callback',c,'Position',[440 140 120 16].*WS,...
	'String',sprintf('%0.4g',SCALE));
c    = ['TYPE = eval(get(gco,''String''));'];
uicontrol(F,'Style','edit','Callback',c,'Position',[440 120 040 016].*WS,...
	'String',sprintf('%0.0f',TYPE));
c    = ['OFFSET = eval(get(gco,''String''));'];
uicontrol(F,'Style','edit','Callback',c,'Position',[440 100 040 016].*WS,...
	'String',sprintf('%0.0f',OFFSET));
c    = ['ORIGIN(1) = eval(get(gco,''String''));'];
uicontrol(F,'Style','edit','Callback',c,'Position',[440 080 040 016].*WS,...
	'String',sprintf('%i',ORIGIN(1)));
c    = ['ORIGIN(2) = eval(get(gco,''String''));'];
uicontrol(F,'Style','edit','Callback',c,'Position',[480 080 040 016].*WS,...
	'String',sprintf('%i',ORIGIN(2)));
c    = ['ORIGIN(3) = eval(get(gco,''String''));'];
uicontrol(F,'Style','edit','Callback',c,'Position',[520 080 040 016].*WS,...
	'String',sprintf('%i',ORIGIN(3)));
c    = ['DESCRIP = get(gco,''String'');'];
uicontrol(F,'Style','edit','Callback',c,'Position',[440 060 120 016].*WS,...
	'String',DESCRIP);



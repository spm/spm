function [s] = spm_hwrite(P,DIM,VOX,SCALE,TYPE,OFFSET,ORIGIN,DESCRIP)
% writes a header
% FORMAT [s] = spm_hwrite(P,DIM,VOX,SCALE,TYPE,OFFSET,ORIGIN,DESCRIP);
%
% P       - filename 	     (e.g 'spm' or 'spm.img')
% DIM     - image size       [i j k [l]] (voxels)
% VOX     - voxel size       [x y z [t]] (mm [sec])
% SCALE   - scale factor
% TYPE    - datatype (integer - see spm_type)
% OFFSET  - offset (bytes)
% ORIGIN  - [i j k] of origin  (default = [0 0 0])
% DESCRIP - description string (default = 'spm compatible')
%
% s       - number of elements successfully written (should be 348)
%___________________________________________________________________________
%
% spm_hwrite writes variables from working memory into a SPM/ANALYZE
% compatible header file.  The 'originator' field of the ANALYZE format has
% been changed to ORIGIN in the SPM version of the header. funused1
% of the ANALYZE format is used for SCALE
%
% see also dbh.h (ANALYZE) spm_hread.m and spm_type.m
%
%__________________________________________________________________________
% %W% %E%


% ensure correct suffix {.hdr} and open header file
%---------------------------------------------------------------------------
P               = P(P ~= ' ');
q    		= length(P);
if P(q - 3) == '.'; P = P(1:(q - 4)); end
P     		= [P '.hdr'];
fid             = fopen(P,'w');

%---------------------------------------------------------------------------
data_type 	= ['dsr      ' 0];

P     		= [P '                  '];
db_name		= [P(1:17) 0];

% set header variables
%---------------------------------------------------------------------------
DIM		= DIM(:)'; if size(DIM,2) < 4; DIM = [DIM 1]; end
VOX		= VOX(:)'; if size(VOX,2) < 4; VOX = [VOX 0]; end
dim		= [4 DIM(1:4) 0 0 0];	
pixdim		= [0 VOX(1:4) 0 0 0];
vox_offset      = OFFSET;
funused1	= SCALE;
glmax		= 1;
glmin		= 0;
bitpix 		= 0;
descrip         = zeros(1,80);
aux_file        = ['none                   ' 0];
origin          = [0 0 0 0 0];

%---------------------------------------------------------------------------
if TYPE == 1;   bitpix = 1;  glmax = 1;        glmin = 0;	end
if TYPE == 2;   bitpix = 8;  glmax = 255;      glmin = 0;	end
if TYPE == 4;   bitpix = 16; glmax = 32767;    glmin = 0;  	end
if TYPE == 8;   bitpix = 32; glmax = (2^31-1); glmin = 0;	end
if TYPE == 16;  bitpix = 32; glmax = 1;        glmin = 0;	end
if TYPE == 64;  bitpix = 64; glmax = 1;        glmin = 0;	end

%---------------------------------------------------------------------------
if nargin >= 7; origin = [ORIGIN(:)' 0 0];  end
if nargin <  8; DESCRIP = 'spm compatible'; end

d          	= 1:min([length(DESCRIP) 79]);
descrip(d) 	= DESCRIP(d);

fseek(fid,0,'bof');

% write (struct) header_key
%---------------------------------------------------------------------------
fwrite(fid,348,		'int32');
fwrite(fid,data_type,	'char' );
fwrite(fid,db_name,	'char' );
fwrite(fid,0,		'int32');
fwrite(fid,0,		'int16');
fwrite(fid,'r',		'char' );
fwrite(fid,'0',		'char' );

% write (struct) image_dimension
%---------------------------------------------------------------------------
fseek(fid,40,'bof');

fwrite(fid,dim,		'int16');
fwrite(fid,'mm',	'char' );
fwrite(fid,0,		'char' );
fwrite(fid,0,		'char' );

fwrite(fid,zeros(1,8),	'char' );
fwrite(fid,0,		'int16');
fwrite(fid,TYPE,	'int16');
fwrite(fid,bitpix,	'int16');
fwrite(fid,0,		'int16');
fwrite(fid,pixdim,	'float');
fwrite(fid,vox_offset,	'float');
fwrite(fid,funused1,	'float');
fwrite(fid,0,		'float');
fwrite(fid,0,		'float');
fwrite(fid,0,		'float');
fwrite(fid,0,		'float');
fwrite(fid,0,		'int32');
fwrite(fid,0,		'int32');
fwrite(fid,glmax,	'int32');
fwrite(fid,glmin,	'int32');

% write (struct) image_dimension
%---------------------------------------------------------------------------
fwrite(fid,descrip,	'char');
fwrite(fid,aux_file,    'char');
fwrite(fid,0,           'char');
fwrite(fid,origin,      'int16');
fwrite(fid,zeros(1,85), 'char');

s   = ftell(fid);
fclose(fid);


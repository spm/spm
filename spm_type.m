function T = spm_type(x)
% translates data type specifiers between SPM & Matlab representations
% FORMAT T = spm_type(x)
% x    - specifier
% T    - type
%_______________________________________________________________________
%
% Format specifiers are based on ANALYZE.  If the input is
% a number then the corresponding matlab string is returned
% If the input is a string then the appropriate TYPE is
% returned
%
% TYPE  bits    bytes   MATLAB  ANALYZE                 UNIX
% 0     0       -       -       DT_NONE                 -
% 1     1       -       uint1   DT_BINARY               -
% 2     8       1       uint8   DT_UNSIGNED_CHAR        (unsigned char)
% 4     16      2       int16   DT_SIGNED_SHORT         (short)
% 8     32      4       int32   DT_SIGNED_INT           (int)
% 16    32      4       float   DT_FLOAT                (float)
% 64    64      8       double  DT_DOUBLE               (double)
%
%_______________________________________________________________________
% %W% John Ashburner, Andrew Holmes %E%

if nargin==0, error('insufficient arguments'), end

if nargin==1
	if ischar(x)
		switch lower(x)
		case 'uint1',  T =  1;
		case 'uint8',  T =  2;
		case 'int16',  T =  4;
		case 'int32',  T =  8;
		case 'float',  T = 16;
		case 'double', T = 64;
		otherwise,     T =  0;
		end
	else
		switch x
		case  1,   T = 'uint1';
		case  2,   T = 'uint8';
		case  4,   T = 'int16';
		case  8,   T = 'int32';
		case 16,   T = 'float';
		case 64,   T = 'double';
		otherwise, T= 'unknown';
		end
	end
end

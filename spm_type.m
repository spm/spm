function [T] = spm_type(x)
% translates data type specifiers within SPM 
% FORMAY [T] = spm_type(x)
% x    - specifier {strings must be at least 5 characters}
% T    - type
%___________________________________________________________________________
%
% Format specifiers are based on ANALYZE.  If the input is
% a number then the corresponding matlab string is returned
% If the input is a string then the appropriate TYPE is
% returned
%
%---------------------------------------------------------------------------
% TYPE  bits    bytes   MATLAB  ANALYZE                 UNIX
% 0     0       -       -       DT_NONE                 -
% 1     1       -       uint1   DT_BINARY               -
% 2     8       1       uint8   DT_UNSIGNED_CHAR        (unsigned char)
% 4     16      2       int16   DT_SIGNED_SHORT         (short)
% 8     32      4       int32   DT_SIGNED_INT           (int)
% 16    32      4       float   DT_FLOAT                (float)
% 64    64      8       double  DT_DOUBLE               (double)
%---------------------------------------------------------------------------
%
%___________________________________________________________________________
% %W% %E%


if isstr(x)
	T = 0;
	x = x(1:5);
	if all(x == 'uint1');  T = 1;  end
	if all(x == 'uint8');  T = 2;  end
	if all(x == 'int16');  T = 4;  end
	if all(x == 'int32');  T = 8;  end
	if all(x == 'float');  T = 16; end
	if all(x == 'doubl');  T = 64; end
else
	T = 'unkmown';
	if (x == 1);  T = 'uint1';  end
	if (x == 2);  T = 'uint8';  end
	if (x == 4);  T = 'int16';  end
	if (x == 8);  T = 'int32';  end
	if (x == 16); T = 'float';  end
	if (x == 64); T = 'double'; end
end

function en = end(a,k,n)
% Overloaded end function for spm_file_array objects.
% _______________________________________________________________________
% @(#)end.m	1.1 John Ashburner 04/06/28

if n>length(a.dim),
	a.dim = [a.dim(1:(k-1)) prod(a.dim(k:end))];
end;
en = a.dim(k);

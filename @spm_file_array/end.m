function en = end(a,k,n)
% Overloaded end function for spm_file_array objects.
% _______________________________________________________________________
% %W% John Ashburner %E%

if n>length(a.dim),
	a.dim = [a.dim(1:(k-1)) prod(a.dim(k:end))];
end;
en = a.dim(k);

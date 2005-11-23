function acc = paccuracy(V,p)
% if ~spm_type(V.dim(4),'intt'),
if ~spm_type(V.dt(1),'intt'),
	acc = 0;
else,
	if size(V.pinfo,2)==1,
		acc = abs(V.pinfo(1,1));
	else,
		acc = abs(V.pinfo(1,p));
	end;
end;

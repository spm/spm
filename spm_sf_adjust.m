function spm_sf_adjust(X)
% adjusts SCALE (scale factors) in a series of headers
% FORMAT spm_sf_adjust(X)
% X  - vector of coeficients that are applied to exisiting SCALE
%___________________________________________________________________________
%
%__________________________________________________________________________
% %W% %E%

X     = X(:);
n     = size(X,1);
P     = spm_get(length(X));
for i = 1:n
	d   = P(i,:);
	[DIM VOX SCALE TYPE OFFSET ORIGIN DESCRIP] = spm_hread(d);
	spm_hwrite(d,DIM,VOX,SCALE*X(i),TYPE,OFFSET,ORIGIN,DESCRIP);
end

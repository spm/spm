function V = spm_write_plane(V,dat,n)
if isfield(V,'n'), n = num2cell([n V.n]); else, n = {n}; end;
S     = struct('type','()','subs',{{':',':',n{:}}});
V.private.dat = subsasgn(V.private.dat,S,dat);

function spm_def2det_ui
%
%_______________________________________________________________________
% John Ashburner $Id$

P    = spm_get(Inf,{'*y_*.img','noexpand'},'Select deformation fields');
n    = size(P,1);
spm_progress_bar('Init',n,'Writing Jacobian Determinants','volumes completed');
for i=1:n,
	Pi = [repmat([deblank(P(i,:)) ','],3,1) num2str([1 2 3]')];
        doit(Pi);
        spm_progress_bar('Set',i);
end;
spm_progress_bar('Clear')
return;
%_______________________________________________________________________

%_______________________________________________________________________
function doit(V)
if ischar(V), V  = spm_vol(V); end;
y1 = spm_load_float(V(1));
y2 = spm_load_float(V(2));
y3 = spm_load_float(V(3));

dt = spm_def2det(y1,y2,y3,V(1).mat);

VO         = V(1);
VO.fname   = prepend(V(1).fname, 'j');
VO.dim(4)  = spm_type('float');
VO.pinfo   = [1 0 0]';
VO.descrip = 'Jacobian determinant';
spm_write_vol(VO,dt);
return;
%_______________________________________________________________________

%_______________________________________________________________________
function out = prepend(in, pre)
[pth,nme,ext,ver] = fileparts(in);
out = fullfile(pth,[pre nme ext ver]);
return;
%_______________________________________________________________________

function spm_def2det_ui
%
%_______________________________________________________________________
% %W% John Ashburner %E%

n       = spm_input('Number of subjects','+0', 'n', '1', 1)';
for i=1:n,
        P{i}    = spm_get(3,'*y?_*.img',['Select deformation field ' num2str(i)]);
end;

spm_progress_bar('Init',n,'Writing Jacobian Determinants','volumes completed');
for i=1:n,
        doit(P{i});
        spm_progress_bar('Set',i);
end;
spm_progress_bar('Clear')
return;
%_______________________________________________________________________

%_______________________________________________________________________
function doit(P)
V  = spm_vol(P);
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

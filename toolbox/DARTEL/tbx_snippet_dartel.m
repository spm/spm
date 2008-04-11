% $Id: tbx_snippet_dartel.m 1381 2008-04-11 19:10:56Z john $
%_______________________________________________________________________
%
% The remainder of this file was appended to an automatically
% generated tbx_cfg_dartel.m file.  Note that the vfiles fields need
% to be modified by hand.
%_______________________________________________________________________

%_______________________________________________________________________

function dep = vout_initial_import(job)
cls = [job.GM, job.WM, job.CSF];
kk = 1;
for k=1:3,
    if cls(k),
        dep(kk)            = cfg_dep;
        dep(kk).sname      = sprintf('Imported Tissue %d', k);
        dep(kk).src_output = substruct('.','cfiles','{}',{':',k});
        dep(kk).tgt_spec   = cfg_findspec({{'filter','nifti'}});
        kk = kk + 1;
    end
end
%_______________________________________________________________________

%_______________________________________________________________________
function chk = check_runjob(job)
n1 = numel(job.images);
n2 = numel(job.images{1});
chk = '';
for i=1:n1,
    if numel(job.images{i}) ~= n2,
        chk = 'Incompatible number of images';
        break;
    end;
end;
%_______________________________________________________________________

%_______________________________________________________________________
function dep = vout_runjob(job)
dep(1)            = cfg_dep;
dep(1).sname      = 'Template';
dep(1).src_output = substruct('.','template','{}',':');
dep(1).tgt_spec   = cfg_findspec({{'filter','nifti'}});

dep(2)            = cfg_dep;
dep(2).sname      = 'Flow Fields';
dep(2).src_output = substruct('.','files','{}',':');
dep(2).tgt_spec   = cfg_findspec({{'filter','nifti'}});
%_______________________________________________________________________

%_______________________________________________________________________
function dep = vout_runjob1(job)
dep            = cfg_dep;
dep.sname      = 'Flow Fields';
dep.src_output = substruct('.','files','{}',':');
dep.tgt_spec   = cfg_findspec({{'filter','nifti'}});
%_______________________________________________________________________

%_______________________________________________________________________
function chk = check_norm(job)
chk = '';
PU = job.flowfields;
PI = job.images;
n1 = numel(PU);
for i=1:numel(PI),
    if numel(PI{i}) ~= n1,
        chk = 'Incompatible number of images';
        break;
    end
end
%_______________________________________________________________________

%_______________________________________________________________________
function dep = vout_norm(job)
if job.jactransf,
    sname = 'Warped Images - Jacobian Transformed';
else
    sname = 'Warped Images';
end
PU    = job.flowfields;
PI    = job.images;
for m=1:numel(PI),
    dep(m)            = cfg_dep;
    dep(m).sname      = sprintf('%s (Image %d)',sname,m);
    dep(m).src_output = substruct('.','files','{}',{':',m});
    dep(m).tgt_spec   = cfg_findspec({{'filter','nifti'}});
end
n = numel(PI);
for m=1:numel(PU),
    dep(m+n)            = cfg_dep;
    dep(m+n).sname      = sprintf('%s (Deformation %d)',sname,m);
    dep(m+n).src_output = substruct('.','files','{}',{m,':'});
    dep(m+n).tgt_spec   = cfg_findspec({{'filter','nifti'}});
end
%_______________________________________________________________________

%_______________________________________________________________________
function dep = vout_invnorm(job)
PU    = job.flowfields;
PI    = job.images;

for m=1:numel(PI),
    dep(m)            = cfg_dep;
    dep(m).sname      = sprintf('Inverse Warped Images (Image %d)',m);
    dep(m).src_output = substruct('.','files','{}',{':',m});
    dep(m).tgt_spec   = cfg_findspec({{'filter','nifti'}});
end
n = numel(PI);
for m=1:numel(PU),
    dep(m+n)            = cfg_dep;
    dep(m+n).sname      = sprintf('Inverse Warped Images (Deformation %d)',m);
    dep(m+n).src_output = substruct('.','files','{}',{m,':'});
    dep(m+n).tgt_spec   = cfg_findspec({{'filter','nifti'}});
end
%_______________________________________________________________________

%_______________________________________________________________________
function dep = vout_jacdet(job)
dep            = cfg_dep;
dep.sname      = 'Jacobian Determinant Fields';
dep.src_output = substruct('.','files','{}',':');
dep.tgt_spec   = cfg_findspec({{'filter','nifti'}});
%_______________________________________________________________________

%_______________________________________________________________________
function chk = check_resids(job)
chk = '';
PU = job.flowfields;
PI = job.images;
n1 = numel(PU);
for i=1:numel(PI),
    if numel(PI{i}) ~= n1,
        chk = 'Incompatible number of images';
        break;
    end
end
%_______________________________________________________________________

%_______________________________________________________________________
function dep = vout_resids(job)
dep = cfg_dep;
dep.sname      = 'Residual Files';
dep.src_output = substruct('.','files','{}',':');
dep.tgt_spec   = cfg_findspec({{'filter','nifti'}});
%_______________________________________________________________________





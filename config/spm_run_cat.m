function spm_run_cat(job)
% SPM job execution function
% takes a harvested job data structure and call SPM functions to perform
% computations on the data.
% Input:
% job    - harvested job data structure (see matlabbatch help)
% Output:
% out    - computation results, usually a struct variable.
%_______________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% John Ashburner
% $Id: spm_run_cat.m 3430 2009-09-29 16:55:38Z guillaume $


V    = spm_vol(strvcat(job.vols{:}));
ind  = cat(1,V.n);
N    = cat(1,V.private);

mx   = -Inf;
mn   = Inf;
for i=1:numel(V),
    dat      = V(i).private.dat(:,:,:,ind(i,1),ind(i,2));
    dat      = dat(isfinite(dat));
    mx       = max(mx,max(dat(:)));
    mn       = min(mn,min(dat(:)));
end;

sf         = max(mx,-mn)/32767;
ni         = nifti;
ni.dat     = file_array('4D.nii',[V(1).dim numel(V)],'INT16-BE',0,sf,0);
ni.mat     = N(1).mat;
ni.mat0    = N(1).mat;
ni.descrip = '4D image';
create(ni);
for i=1:size(ni.dat,4),
    ni.dat(:,:,:,i) = N(i).dat(:,:,:,ind(i,1),ind(i,2));
    spm_get_space([ni.dat.fname ',' num2str(i)], V(i).mat);
end;


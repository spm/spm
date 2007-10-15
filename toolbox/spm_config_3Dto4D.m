function c = spm_config_3Dto4D(varargin)
% Configuration file for concatenation jobs
%_______________________________________________________________________
% Copyright (C) 2005 Wellcome Department of Imaging Neuroscience

% John Ashburner
% $Id: spm_config_3Dto4D.m 946 2007-10-15 16:36:06Z john $

vols.type = 'files';
vols.name = '3D Volumes';
vols.tag  = 'vols';
vols.num  = [1 Inf];
vols.filter = 'image';
vols.help = {'Select the volumes to concatenate'};

c.type = 'branch';
c.name = '3D to 4D';
c.tag  = 'cat';
c.val  = {vols};
c.prog = @spm_3Dto4D;
c.help = {'Concatenate a number of 3D volumes into a single 4D file.',...
'Note that output time series are stored as big-endian int16.'};
%_______________________________________________________________________

%_______________________________________________________________________
function spm_3Dto4D(job)

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


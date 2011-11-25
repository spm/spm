function out = spm_shoot_divergence(job)
% Compute divergences from velocity fields
% FORMAT spm_shoot_divergence(job)
% job.velocities - Filenames of initial velocity fields
%
% Copyright (C) 2011 Wellcome Trust Centre for Neuroimaging

% John Ashburner
% $Id: spm_shoot_divergence.m 4573 2011-11-25 23:01:01Z john $

P = strvcat(job.velocities);
out = cell(size(P,1),1);
for i=1:size(P,1)
    Nii   = nifti(deblank(P(i,:)));
    d     = size(Nii.dat);
    krn   = {[-1;0;1]/2,[-1,0,1]/2,reshape([-1,0,1],[1 1 3])/2};
    dv    = zeros(d(1:3));
    for dm=1:3,
        dv = dv - convn(Nii.dat(:,:,:,1,dm),krn{dm},'same');
    end
    dv(:,:,1)=0; dv(:,:,end)=0;
    dv(:,1,:)=0; dv(:,end,:)=0;
    dv(1,:,:)=0; dv(end,:,:)=0;
    Nio            = Nii;
    [pth,nam,ext]  = fileparts(deblank(P(i,:)));
    Nio.dat.fname  = fullfile(pth,['d' nam(1:end) ext]);
    Nio.dat.dim    = d(1:3);
    Nio.descrip    = 'Divergence Field';
    create(Nio);
    Nio.dat(:,:,:) = dv;
    out{i}         = Nio.dat.fname;
end


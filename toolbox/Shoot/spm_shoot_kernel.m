function out = spm_shoot_kernel(job)
% Generate kernel matrix from initial velocity fields
% FORMAT spm_shoot_kernel(job)
% job.velocities - Initial velocity fields
% job.dotprod    - Part of filename for results
%
% k(x_1,x_2) = <x_1,L x_2> = <L x_1, x_2>
%
% This is very slow, and is not in a form that would be
% suited to weighting according to location in the image.
% For this, the "square root" of L would need to be used
% in order to convert the flow fields into (e.g.) their
% Jacobian tensor fields.  For linear elasticity, this
% field would be decomposed by J = (J+J')/2 + (J-J')/2.
% The elements of the symetric part (along with its trace)
% would then be used to generate the kernel.
%_______________________________________________________________________
% Copyright (C) 2009 Wellcome Trust Centre for Neuroimaging

% John Ashburner
% $Id: spm_shoot_kernel.m 4026 2010-07-29 13:45:50Z john $

defs   = spm_shoot_defaults;
rparam = defs.rparam;
P      = strvcat(job.velocities);
[pth,nam,ext] = fileparts(job.dotprod);
ofname = fullfile(pwd,['dp_' nam '.mat']);
out.fname = {ofname};

N   = nifti(P);
dm  = size(N(1).dat);
vx  = sqrt(sum(N(1).mat(1:3,1:3).^2));
prm = [rparam(1) vx rparam(2:end)];
n   = numel(N);
K   = zeros(n,n);
spm_progress_bar('Init',n*n,'Generating matrix','Elements done');
for i=1:n,
    x1 = single(squeeze(N(i).dat(:,:,:,end,:)));
    x1 = dartel3('vel2mom',x1,prm);
    for j=i:n,
        x2       =squeeze(N(j).dat(:,:,:,end,:,:));
        d        = x1(:)'*x2(:);
        K(i,j)   = d;
        K(j,i)   = d;

        if ~rem(j,8)
            spm_progress_bar('Set',...
                n*n-(n+1-i)*(n+1-i)+(j-i)*2+1);
        end
    end
    input        = job;
    input.rparam = rparam;
    typ          = 'initvel';
    save(ofname,'K','input','typ');
end
spm_progress_bar('Clear');


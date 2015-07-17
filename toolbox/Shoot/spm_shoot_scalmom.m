function out = spm_shoot_scalmom(job)
% Generate ``scalar momenta'' for use as features in pattern recognition
% FORMAT out = spm_shoot_scalmom(job)
%
% See:
% Singh, Nikhil, P. Fletcher, J. Preston, Linh Ha, Richard King,
% J. Marron, Michael Wiener, and Sarang Joshi. "Multivariate
% statistical analysis of deformation momenta relating anatomical
% shape to neuropsychological measures." Medical Image Computing
% and Computer-Assisted Intervention-MICCAI 2010 (2010): 529-537.
%
% Singh, Nikhil, Angela Wang, Preethi Sankaranarayanan, P. Fletcher,
% and Sarang Joshi. "Genetic, Structural and Functional Imaging
% Biomarkers for Early Detection of Conversion from MCI to AD."
% Medical Image Computing and Computer-Assisted Intervention-MICCAI
% 2012 (2012): 132-140.
%__________________________________________________________________________
% Copyright (C) 2013-2015 Wellcome Trust Centre for Neuroimaging

% John Ashburner
% $Id: spm_shoot_scalmom.m 6501 2015-07-17 14:32:09Z spm $

Pt = strvcat(job.template);
Nt = nifti(Pt);
d  = size(Nt.dat);
Py = strvcat(job.deformations);
N  = size(Py,1);
if numel(job.images)~=d(4)-1,
    error('Incompatible number of components between template and images.');
end
Pj = strvcat(job.jacobians);
if size(Pj,1) ~= N,
    error('Incompatible number of deformations and Jacobians.');
end

Pc = cell(1,numel(job.images));
for i=1:numel(job.images),
    Pc{i} = strvcat(job.images{i});
    if size(Pc{i},1) ~= N,
        error('Incompatible number of deformations and images.');
    end
end

fwhm = job.fwhm;

spm_progress_bar('Init',size(Py,1),'Computing Scalar Momenta','Subjects done');
out.scalmom = cell(size(Py,1),1);

R  = null(ones(1,d(4)));     % Weights for linear combination of momentum
A  = zeros([d(1:3),d(4)-1]); % Linear combination of momentum

for j=1:size(Py,1), % Loop over subjects

    Ny = nifti(Py(j,:)); % Header info of deformation and Jacobians
    Nj = nifti(Pj(j,:));

    % Load all the imported tissues
    F  = cell(d(4)-1,1);
    for i=1:(d(4)-1),
        Nc   = nifti(Pc{i}(j,:));
        if i==1,
            if sum(sum((Nc.mat - Nc.mat0).^2)) > 1e-4 && sum(sum((Nc.mat - Ny.mat).^2)) == 0,
                % An "imported image"
                Mat = Nc.mat0;
            else
                % A more typical image
                Mat = Nc.mat;
            end
        end
        F{i} = single(Nc.dat(:,:,:));
    end

    for z=1:d(3) % Loop over slices

        % Load deformation and make it map to voxels instead of mm
        y  = reshape(affind(single(Ny.dat(:,:,z,:,:)),inv(Mat)),[d(1:2),1,3]);

        % Load Jacobian determinants
        jd = squeeze(single(Nj.dat(:,:,z)));

        % Data for all tissues
        x = zeros(d([1 2 4]),'single');
   
        % Background class (from which other classes are subtracted)
        fe = ones(d(1:2),'single');

        % Loop over imported data
        for i=1:d(4)-1,
            f  = spm_diffeo('samp',F{i},y);   % Warp the imported tissue
            f(~isfinite(f)) = 0;              % Assume all values are zero outside FOV
            fe = fe-f;                        % Subtract from background
            t  = single(Nt.dat(:,:,z,i));     % Slice of template
            x(:,:,i) = (f-t).*jd;             % Compute scalar momentum
        end

        % Deal with background class (no imported background)
        t = single(Nt.dat(:,:,z,d(4))); % Background slice of template
        x(:,:,d(4))     = (fe-t).*jd;   % Compute scalar momentum (background)
       %x(~isfinite(x)) = 0;            % Remove NaNs - should not be needed

        % There is redundancy in using all tissues because they sum to 1 at each
        % voxel. Reduce to N-1 tissues by projecting into the null space.
        for j1=1:(d(4)-1),
            A(:,:,z,j1) = 0;
            for j2=1:d(4),
                A(:,:,z,j1) = A(:,:,z,j1) + R(j2,j1)*x(:,:,j2);
            end
        end
    end

    % Write output data
    [pth,nam]  = fileparts(Pc{1}(j,:));
    No         = nifti;
    No.dat     = file_array(fullfile(pth,['a_' nam '.nii']),size(A),'FLOAT32',0,1,0);
    No.mat0    = Nt.mat0;
    No.mat     = Nt.mat;
    No.descrip = sprintf('Scalar Momentum (%g)', fwhm);
    create(No);
    if fwhm>0,
        vx = sqrt(sum(No.mat(1:3,1:3).^2));
        for j1=1:(d(4)-1),
            a = A(:,:,:,j1);
            spm_smooth(a,a,fwhm./vx);
            A(:,:,:,j1) = a;
        end
    end
    No.dat(:,:,:,:,:) = A;

    out.scalmom{j}    = No.dat.fname;
   %fprintf('%d\t%s\n', j, nam);
    spm_progress_bar('Set',j);
end
spm_progress_bar('Clear');


function y1 = affind(y0,M)
y1 = zeros(size(y0),'single');
for d=1:3,
    y1(:,:,:,d) = y0(:,:,:,1)*M(d,1) + y0(:,:,:,2)*M(d,2) + y0(:,:,:,3)*M(d,3) + M(d,4);
end


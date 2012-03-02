Py    = spm_select(Inf,'^y_.*\.nii$', 'Select deformations');
N     = size(Py,1);
Pj    = spm_select(N,'^j_.*\.nii$', 'Select Jacobians');
Pc{1} = spm_select(N,'^rc1.*\.nii$','Select imported GM');
Pc{2} = spm_select(N,'^rc2.*\.nii$','Select imported WM');
Pt    = spm_select(1,'nifti','Select template');

Nt = nifti(Pt);
d  = size(Nt.dat);
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
        y  = reshape(affind(single(Ny.dat(:,:,z,:,:)),inv(Mat)),[d(1:2),1,d(4)]);

        % Load Jacobian determinants
        jd = squeeze(single(Nj.dat(:,:,z)));

        % Data for all tissues
        x = zeros(d([1 2 4]),'single');
   
        % Background class (from which other classes are subtracted)
        fe = ones(d(1:2),'single');

        % Loop over imported data
        for i=1:d(4)-1,
            f  = shoot3('samp',F{i},y);   % Warp the imported tissue
            fe = fe-f;                    % Subtract from background
            t  = single(Nt.dat(:,:,z,i)); % Slice of template
            x(:,:,i) = (f-t).*jd;         % Compute scalar momentum
        end

        % Deal with background class (no imported background)
        t = single(Nt.dat(:,:,z,d(4))); % Background slice of template
        x(:,:,d(4))     = (fe-t).*jd;   % Compute scalar momentum (background)
        x(~isfinite(x)) = 0;            % Remove NaNs

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
    No.dat     = file_array(fullfile('.',['a_' nam '.nii']),size(A),'FLOAT32',0,1,0);
    No.mat0    = Nt.mat0;
    No.mat     = Nt.mat;
    No.descrip = 'Scalar Momentum';
    create(No);
    No.dat(:,:,:,:,:) = A;

    fprintf('%d\t%s\n', j, nam);
end



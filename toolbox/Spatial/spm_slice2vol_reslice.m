function spm_slice2vol_reslice(Nii,Q,fwhm)
% Slice-to-volume alignment reslicing
% FORMAT spm_slice2vol_reslice(Nii,Q,fwhm)
%
% Nii  - NIfTI data structure encoding volumes to align
%        Most all have the same dimensions
% Q    - A 3D array of slicewise motion parameters
% fwhm - Smoothing FWHM (mm)
%__________________________________________________________________________

% John Ashburner
% Copyright (C) 2021-2022 Wellcome Centre for Human Neuroimaging

if nargin<3 || fwhm==0
    fwhm = 0;
    pref = 'r';
    str  = 'slice-to-vol aligned';
else
    pref = 'sr';
    str  = sprintf('smoothed (%g mm), slice-to-vol aligned', fwhm);
end
if isa(Nii,'char'), Nii=nifti(Nii); end
Mat = Nii(1).mat;

offset = 0;
for i=1:numel(Nii)
    Nio             = nifti;
    Nio.dat         = Nii(i).dat;
    [pth,nam,ext]   = fileparts(Nii(i).dat.fname);
    Nio.dat.fname   = fullfile(pth,[pref nam ext]);
    Nio.mat         = Mat;
    Nio.mat0        = Mat;
    Nio.mat_intent  = 2;
    Nio.mat0_intent = 2;
    Nio.descrip     = str;
    create(Nio);
    d4              = size(Nio.dat,4);
    Qi              = Q(:,:,(offset+1):(offset+d4));
    offset          = offset + d4;
    reslice_volumes(Nii(i),Qi,Nio,fwhm);
end


%==========================================================================
function reslice_volumes(Nii,Q,Nio,fwhm)
%==========================================================================
d   = [size(Nii.dat) 1 1];
d   = d(1:4);
Mat = Nii.mat;
B   = bases;
[i1,i2] = ndgrid(single(1:d(1)),single(1:d(2)));
fprintf('writing:');
for n=1:d(4)
    if ~rem(n,10), fprintf('.'); end
    phi = zeros([d(1:3),3],'single');
    for i3=1:d(3)
        q             = Q(:,i3,n);
        Mr            = spm_dexpm(q,B);
        Ms            = Mat\Mr*Mat;
        phi(:,:,i3,1) = Ms(1,1)*i1 + Ms(1,2)*i2 + (Ms(1,3)*i3+Ms(1,4));
        phi(:,:,i3,2) = Ms(2,1)*i1 + Ms(2,2)*i2 + (Ms(2,3)*i3+Ms(2,4));
        phi(:,:,i3,3) = Ms(3,1)*i1 + Ms(3,2)*i2 + (Ms(3,3)*i3+Ms(3,4));
    end
    f      = single(Nii.dat(:,:,:,n));
    msk    = (f==0);
    f(msk) = NaN;
    [G,H]  = spm_diffeo('push',f,phi);

    % Smoothing
    vx  = sqrt(sum(Mat(1:3,1:3).^2));
    krn = max(fwhm./vx,0.25);
    spm_smooth(G,G,krn); % Side effects
    spm_smooth(H,H,krn); % Side effects

    f1     = H.\G;
    Nio.dat(:,:,:,n) = f1;
    md = ceil(size(f,3)/2);
   %imagesc([f(:,:,md)' f1(:,:,md)']); axis image xy off; drawnow
end
fprintf('\n');


%==========================================================================
function B = bases
%==========================================================================
% SE(3) - Special Euclidean
B        = zeros(4,4,6);
B(1,4,1) =  1;
B(2,4,2) =  1;
B(3,4,3) =  1;
B(1,2,4) =  1;
B(2,1,4) = -1;
B(1,3,5) =  1;
B(3,1,5) = -1;
B(2,3,6) =  1;
B(3,2,6) = -1;

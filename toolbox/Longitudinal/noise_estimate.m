function noise = noise_estimate(Scans)

noise = zeros(size(Scans,1),1);
for i=1:size(Scans,1),
    Nii=nifti(Scans(i,:));
    if spm_type(Nii.dat.dtype(1:(end-3)),'intt'),
        f  = Nii.dat(:,:,:);
        f(f==max(f(:)))=0;
        x=0:Nii.dat.scl_slope:max(f(:));
        [h,x]=hist(f(f~=0),x);
    else
        f  = Nii.dat(:,:,:);
        x=(0:1023)*(max(f(:))/1023);
        f(f==max(f(:)))=0;
        [h,x]=hist(f(f~=0 & isfinite(f)),x);
    end
    [mg,nu,sd]=rice_mixture(h(:),x(:),2);
    noise(i) = min(sd);
    fprintf('%d %g\n', i, noise(i));
end


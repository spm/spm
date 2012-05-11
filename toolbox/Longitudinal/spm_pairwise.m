function out = spm_pairwise(job)
% Longitudinal registration of image pairs
% FORMAT out = spm_pairwise(job)
% See tbx_cfg_longitudinal.m for a description of the various fields. 
%_______________________________________________________________________
% Copyright (C) 2012 Wellcome Trust Centre for Neuroimaging

% John Ashburner
% $Id: spm_pairwise.m 4738 2012-05-11 16:41:35Z ged $

N = numel(job.vols1);
if numel(job.vols2) ~= N, error('Incompatible numbers of scans.'); end
if numel(job.tdif) == 1,
    tdif = repmat(abs(job.tdif),N,1);
else
    if numel(job.tdif) ~= N, error('Incompatible numbers of time differences.'); end
    tdif = abs(job.tdif(:));
end
if any(tdif > 50), error('Time differences should be in years.'); end;

if ~isfinite(job.noise),
    % Make an estimate of the scanner noise
    Scans = strvcat(job.vols1{:},job.vols2{:});
    noise = mean(noise_estimate(Scans).^2);
    fprintf('Noise estimate = %g\n', sqrt(noise));
else
    noise = job.noise;
end
prec    = 1/noise;
bparam  = [0 0 job.bparam];
wparam0 = job.wparam;

output = {};
if job.write_avg, output = {output{:}, 'wavg'}; end
if job.write_jac, output = {output{:},  'jac'}; end
if job.write_div, output = {output{:},  'div'}; end
if job.write_def, output = {output{:}, 'wdef'}; end

for i=1:numel(tdif),
    wparam = kron(wparam0,1./(abs(tdif(i)/2)+1/365));
    sparam = round(3*abs(tdif(i)/2)+2);
    Nii    = nifti(strvcat(job.vols1{i},job.vols2{i}));
    [pth,nam1] = fileparts(Nii(1).dat.fname);
    [pth,nam2] = fileparts(Nii(2).dat.fname);
    fprintf('*** %s <=> %s ***\n', nam1, nam2);

    dat    = spm_groupwise_ls(Nii, output, prec, wparam, bparam, sparam);

    if isfield(dat,'jac')
        d           = [size(dat.jac{1}) 1]; d = d(1:3);
        nam         = fullfile(pth,['jd_' nam1 '_' nam2 '.nii']);
        Nio         = nifti;
        Nio.dat     = file_array(nam,d,'float32',0,1,0);
        Nio.mat     = dat.mat;
        Nio.mat0    = Nio.mat;
        Nio.mat_intent  = 'Aligned';
        Nio.mat0_intent = Nio.mat_intent;
        Nio.descrip = 'Jacobian Difference';
        create(Nio);
        Nio.dat(:,:,:) = dat.jac{2} - dat.jac{1};
        out.jac{i}     = nam;
    end

    if isfield(dat,'div')
        d           = [size(dat.div{1}) 1]; d = d(1:3);
        nam         = fullfile(pth,['dv_' nam1 '_' nam2 '.nii']);
        Nio         = nifti;
        Nio.dat     = file_array(nam,d,'float32',0,1,0);
        Nio.mat     = dat.mat;
        Nio.mat0    = Nio.mat;
        Nio.mat_intent  = 'Aligned';
        Nio.mat0_intent = Nio.mat_intent;
        Nio.descrip = 'Div';
        create(Nio);
        Nio.dat(:,:,:) = dat.div{2} - dat.div{1};
        out.div{i}     = nam; 
    end

    if job.write_avg,
        out.avg{i} = dat.avg;
    end;
    if job.write_def,
        out.def1{i} = dat.def{1};
        out.def2{i} = dat.def{2};
    end;

    clear dat
end
return


function noise = noise_estimate(Scans)
noise = zeros(size(Scans,1),1);
for i=1:size(Scans,1),
    Nii = nifti(Scans(i,:));
    if spm_type(Nii.dat.dtype(1:(end-3)),'intt'),
        f      = Nii.dat(:,:,:);
        f(f==max(f(:))) = 0;
        x      = 0:Nii.dat.scl_slope:max(f(:));
        [h,x]  = hist(f(f~=0),x);
    else
        f      = Nii.dat(:,:,:);
        x      = (0:1023)*(max(f(:))/1023);
        f(f==max(f(:))) = 0;
        [h,x]  = hist(f(f~=0 & isfinite(f)),x);
    end
    [mg,nu,sd] = spm_rice_mixture(h(:),x(:),2);
    noise(i)   = min(sd);
    %fprintf('%d %g\n', i, noise(i));
end


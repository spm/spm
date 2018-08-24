function out = spm_pairwise(job)
% Longitudinal registration of image pairs
% FORMAT out = spm_pairwise(job)
% See tbx_cfg_longitudinal.m for a description of the various fields. 
%_______________________________________________________________________
% Copyright (C) 2012 Wellcome Trust Centre for Neuroimaging

% John Ashburner
% $Id: spm_pairwise.m 7408 2018-08-24 14:54:57Z john $

N = numel(job.vols1);
if numel(job.vols2) ~= N, error('Incompatible numbers of scans.'); end
if numel(job.tdif) == 1
    tdif = repmat(abs(job.tdif),N,1);
else
    if numel(job.tdif) ~= N, error('Incompatible numbers of time differences.'); end
    tdif = abs(job.tdif(:));
end
if any(tdif > 50), error('Time differences should be in years.'); end

if numel(job.noise)==1
    noise = repmat(job.noise,[N 2]);
elseif size(job.noise,2) ~= 2
    error('Incompatible numbers of noise estimates for each subject.');
elseif size(job.noise,1) ~= N
    error('Incompatible numbers of noise estimates and subjects.');
else
    noise = job.noise;
end

for i=find(~isfinite(noise(:,1)))'
    % Make an estimate of the scanner noise
    noise(i,1) = spm_noise_estimate(job.vols1{i});
    fprintf('Estimated noise sd for "%s" = %g\n', job.vols1{i}, noise(i,1));
end
for i=find(~isfinite(noise(:,2)))'
    % Make an estimate of the scanner noise
    noise(i,2) = spm_noise_estimate(job.vols2{i});
    fprintf('Estimated noise sd for "%s" = %g\n', job.vols2{i}, noise(i,2));
end

bparam  = [0 0 job.bparam];
wparam0 = job.wparam;

output = {};
if job.write_avg, output = [output, {'wavg'}]; end
if job.write_jac, output = [output, {'jac'} ];  end
if job.write_div, output = [output, {'div'} ];  end
if job.write_def, output = [output, {'wdef'}]; end

for i=1:numel(tdif)
    wparam = kron(wparam0,1./(abs(tdif(i)/2)+1/365));
    sparam = round(3*abs(tdif(i)/2)+2);
    Nii    = nifti(strvcat(job.vols1{i},job.vols2{i}));
    [pth1,nam1] = fileparts(Nii(1).dat.fname);
    [pth2,nam2] = fileparts(Nii(2).dat.fname);
    fprintf('*** %s <=> %s ***\n', nam1, nam2);

    prec   = noise(i,:).^(-2);
    dat    = spm_groupwise_ls(Nii, output, prec, wparam, bparam, sparam);

    if isfield(dat,'jac')
        d           = [size(dat.jac{1}) 1]; d = d(1:3);
        nam         = fullfile(pth1,['jd_' nam1 '_' nam2 '.nii']);
        Nio         = nifti;
        Nio.dat     = file_array(nam,d,'float32',0,1,0);
        Nio.mat     = dat.mat;
        Nio.mat0    = Nio.mat;
        Nio.mat_intent  = 'Aligned';
        Nio.mat0_intent = Nio.mat_intent;
        Nio.descrip = 'Jacobian Difference';
        create(Nio);
        Nio.dat(:,:,:) = (dat.jac{2} - dat.jac{1})/tdif(i);
        out.jac{i}     = nam;
    end

    if isfield(dat,'div')
        d           = [size(dat.div{1}) 1]; d = d(1:3);
        nam         = fullfile(pth1,['dv_' nam1 '_' nam2 '.nii']);
        Nio         = nifti;
        Nio.dat     = file_array(nam,d,'float32',0,1,0);
        Nio.mat     = dat.mat;
        Nio.mat0    = Nio.mat;
        Nio.mat_intent  = 'Aligned';
        Nio.mat0_intent = Nio.mat_intent;
        Nio.descrip = 'Div';
        create(Nio);
        Nio.dat(:,:,:) = (dat.div{2} - dat.div{1})/tdif(i);
        out.div{i}     = nam; 
    end

    if job.write_avg
        out.avg{i} = dat.avg;
    end
    if job.write_def
        out.def1{i} = dat.def{1};
        out.def2{i} = dat.def{2};
    end

    clear dat
end
return


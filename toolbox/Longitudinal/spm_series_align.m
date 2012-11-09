function out = spm_series_align(job)
% Longitudinal registration of image series
% FORMAT out = spm_series_align(job)
%_______________________________________________________________________
% Copyright (C) 2012 Wellcome Trust Centre for Neuroimaging

% John Ashburner
% $Id: spm_series_align.m 5044 2012-11-09 13:40:35Z john $

N = numel(job.vols);
tim = job.times(:);
if numel(tim) ~= N,
    error('Incompatible numbers of times and scans.');
end
if any(abs(diff(tim)) > 50),
    error('Time differences should be in years.');
end;

if numel(job.noise)==1,
    noise = repmat(job.noise,[N,1]);
elseif numel(job.noise) ~= N,
    error('Incompatible numbers of noise estimates and scans.');
else
    noise = job.noise(:);
end
for i=find(~isfinite(noise(:)))',
    % Make an estimate of the scanner noise
    noise(i,1) = spm_noise_estimate(job.vols{i});
    fprintf('Estimated noise sd for "%s" = %g\n', job.vols{i}, noise(i,1));
end
prec   = noise.^(-2);


bparam    = [0 0 job.bparam];
wparam0   = job.wparam;

midtim = median(tim);
tim    = tim - midtim;
wparam = kron(wparam0,1./(abs(tim)+1/365));
sparam = round(3*abs(tim)+2);
Nii    = nifti(strvcat(job.vols));

output = {};
if job.write_avg, output = [output, {'wavg'}]; end
if job.write_jac, output = [output, {'wjac'} ];  end
if job.write_div, output = [output, {'wdiv'} ];  end
if job.write_def, output = [output, {'wdef'}]; end

out    = spm_groupwise_ls(Nii, output, prec, wparam, bparam, sparam);
return


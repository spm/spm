function job = spm_rewrite_job(job)
% Rewrite a job for SPM12
%__________________________________________________________________________
% Copyright (C) 2012 Wellcome Trust Centre for Neuroimaging

% $Id: spm_rewrite_job.m 4904 2012-09-06 15:08:56Z guillaume $


try
    job.spatial.preproc.data;
    fprintf('Conversion Segment -> Old Segment\n');                     %-#
    job = struct('tools', struct('oldseg', job.spatial.preproc));
end

try
    job.spatial.normalise.est.subj(1).source ;
    fprintf('Conversion Normalise:Est -> Old Normalise:Est\n');         %-#
    job = struct('tools', struct('oldnorm', job.spatial.normalise));
end

try
    job.spatial.normalise.write.subj(1).matname;
    fprintf('Conversion Normalise:Write -> Old Normalise:Write\n');     %-#
    job = struct('tools', struct('oldnorm', job.spatial.normalise));
end

try
    job.spatial.normalise.estwrite.subj(1).source;
    fprintf('Conversion Normalise:EstWrite -> Old Normalise:EstWrite\n');%-#
    job = struct('tools', struct('oldnorm', job.spatial.normalise));
end

function job = spm_rewrite_job(job)
% Rewrite a job for SPM12
%__________________________________________________________________________
% Copyright (C) 2012-2013 Wellcome Trust Centre for Neuroimaging

% $Id: spm_rewrite_job.m 5813 2013-12-20 18:54:10Z guillaume $


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

try
    job.tools.preproc8;
    fprintf('Conversion Tools:New Segment -> Spatial:Segment\n');       %-#
    job = struct('spatial',struct('preproc',job.tools.preproc8));
end

try
    for i=1:numel(job.stats.con.consess)
        try
            con = job.stats.con.consess{i}.tcon.convec;
            job.stats.con.consess{i}.tcon = rmfield(job.stats.con.consess{i}.tcon,'convec');
            job.stats.con.consess{i}.tcon.weights = con;
        end
        try
            con = job.stats.con.consess{i}.fcon.convec;
            job.stats.con.consess{i}.fcon = rmfield(job.stats.con.consess{i}.fcon,'convec');
            job.stats.con.consess{i}.fcon.weights = con;
        end
        try
            job.stats.con.consess{i}.fcon.weights{1};
            fprintf('Conversion to new syntax: Contrast Manager:F-contrast\n'); %-#
            try
                con = cat(1,job.stats.con.consess{i}.fcon.weights{:});
            catch
                fprintf('Error concatenating F-contrast vectors.\n');   %-#
            end
            job.stats.con.consess{i}.fcon.weights = con;
        end
    end
end

try
    if isequal(job.stats.results.print, true)
        job.stats.results.print = spm_get_defaults('ui.print');
    end
end

try
    job.tools.sendmail;
    job = struct('util', job.tools);
end

try
    job.util.spm_surf;
    job = struct('util', struct('render', struct('extract',job.util.spm_surf)));
end

try
    job.util.dicom;
    job.util = struct('import', job.util);
end
try
    job.util.minc;
    job.util = struct('import', job.util);
end
try
    job.util.ecat;
    job.util = struct('import', job.util);
end

try
    util = char(fieldnames(job.util));
    if ismember(util, {'cdir','md','deletefiles','movefile'})
        ws = warning('off','backtrace');
        warning(['spm.util.%s was DEPRECATED in SPM8 and has now been REMOVED.\n' ...
            'Please update your batches to use BasicIO instead.'],util);
        warning(ws);
        %job = struct([]);
    end
end

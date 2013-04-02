function results = spm_cfg_eeg_inv_results
% Configuration file for exporting results of source reconstruction
%__________________________________________________________________________
% Copyright (C) 2010-2013 Wellcome Trust Centre for Neuroimaging

% Vladimir Litvak
% $Id: spm_cfg_eeg_inv_results.m 5377 2013-04-02 17:07:57Z vladimir $


D        = cfg_files;
D.tag    = 'D';
D.name   = 'M/EEG datasets';
D.filter = 'mat';
D.num    = [1 Inf];
D.help   = {'Select the M/EEG mat files.'};

val         = cfg_entry;
val.tag     = 'val';
val.name    = 'Inversion index';
val.strtype = 'n';
val.num     = [1 1];
val.help    = {'Index of the cell in D.inv where the inversion results are stored.'};
val.val     = {1};

woi         = cfg_entry;
woi.tag     = 'woi';
woi.name    = 'Time window of interest';
woi.strtype = 'r';
woi.num     = [Inf 2];
woi.val     = {[-Inf Inf]};
woi.help    = {'Time window to average over (ms).'};

foi         = cfg_entry;
foi.tag     = 'foi';
foi.name    = 'Frequency window of interest';
foi.strtype = 'r';
foi.num     = [1 2];
foi.val     = {[0 0]};
foi.help    = {'Frequency window (Hz).'};

ctype        = cfg_menu;
ctype.tag    = 'ctype';
ctype.name   = 'Contrast type';
ctype.help   = {'Contrast type: evoked activity, induced activity or single trials.'};
ctype.labels = {'Evoked', 'Induced', 'Single trials'};
ctype.values = {'evoked', 'induced', 'trials'};
ctype.val    = {'evoked'};

space        = cfg_menu;
space.tag    = 'space';
space.name   = 'Output space';
space.help   = {'Output space to write the results in.'};
space.labels = {'MNI', 'Native'};
space.values = {1, 0};
space.val    = {1};

fmt        = cfg_menu;
fmt.tag    = 'format';
fmt.name   = 'Output file format';
fmt.help   = {'Output file format.'};
fmt.labels = {'Image', 'Mesh'};
fmt.values = {'image', 'mesh'};
fmt.val    = {'image'};

smoothing         = cfg_entry;
smoothing.tag     = 'smoothing';
smoothing.name    = 'Cortical smoothing';
smoothing.help    = {'Number of iterations for cortical smoothing.'};
smoothing.strtype = 'r';
smoothing.num     = [1 1];
smoothing.val     = {8};

results          = cfg_exbranch;
results.tag      = 'results';
results.name     = 'Inversion results';
results.val      = {D, val, woi, foi, ctype, space, fmt, smoothing};
results.help     = {'Export results of imaging source reconstruction.'};
results.prog     = @run_results;
results.vout     = @vout_results;
results.modality = {'EEG'};


%==========================================================================
function  out = run_results(job)

contrast = [];
contrast.fboi = job.foi;
contrast.type = job.ctype;
contrast.space = job.space;
contrast.format = job.format;
contrast.smoothing = job.smoothing;

files = {};

for i = 1:numel(job.D)
    D = spm_eeg_load(job.D{i});
    
    D.val = job.val;      
            
    if ~isfield(D.inv{D.val}, 'inverse') || ~isfield(D.inv{D.val}.inverse, 'J')
        error('Imaging source reconstruction is missing for subject %d.', i);
    end
    
    contrast.woi = fix(sort(job.woi,2));
    
    D.inv{D.val}.contrast = contrast;
    
    D = spm_eeg_inv_results(D);
    D = spm_eeg_inv_Mesh2Voxels(D);
    
    save(D);
    
    fname = D.inv{D.val}.contrast.fname;
    for j = 1:numel(fname)
        if iscell(fname{j})
            for k = 1:numel(fname{j})
                files = [files; fname{j}(k)];
            end
        else
             files = [files; fname(j)];
        end
    end
end

out.files = files;


%==========================================================================
function dep = vout_results(job)
if strcmpi(job.format,'image')
    dep(1)            = cfg_dep;
    dep(1).sname      = 'Exported source reconstruction (images)';
    dep(1).src_output = substruct('.','files');
    dep(1).tgt_spec   = cfg_findspec({{'filter','image','strtype','e'}});
else
    dep(1)            = cfg_dep;
    dep(1).sname      = 'Exported source reconstruction (meshes)';
    dep(1).src_output = substruct('.','files');
    dep(1).tgt_spec   = cfg_findspec({{'filter','mesh','strtype','e'}});
end

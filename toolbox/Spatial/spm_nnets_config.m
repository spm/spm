function nnets = spm_nnets_config
% SPM Configuration file for various neural network tools
%__________________________________________________________________________

% John Ashburner
% Copyright (C) 2026 Wellcome Centre for Human Neuroimaging

    dep      = cfg_exbranch;
    dep.tag  = 'nndepth';
    dep.name = 'Cortical depth';
    dep.val  = @depth_cfg;
    dep.prog = @(job)spm_run_depth('run',job);
    dep.vout = @(job)spm_run_depth('vout',job);
    dep.help = {'Generate approximate maps of cortical depth.','',...
               ['Uses a form of U-net to predict approximate cortical depth from maps of ' ...
                'GM and WM, which may be useful for layer-based fMRI analyses. ' ...
                'Output values are in the range of 0 (pial surface) to 1 (GM/WM interface). ' ...
                'The topology of the thresholded depth map may not be exactly correct, ' ...
                'but the errors may be small enough for layer-based fMRI applications. ' ...
                'This could be done by selecting voxels between various depth thresholds.'], '', ...
               ['The network is applied at 1 mm isotropic resolution, and is intended to work on data in any orientation.'], '', ...
               ['This functionality is implemented in the Julia language and several libraries, which will be ' ...
                'installed in a toolbox/Spatial/Apps folder when the code is first run.']};

    seg      = cfg_exbranch;
    seg.tag  = 'nnseg';
    seg.name = 'Segment approx';
    seg.val  = @nnseg_cfg;
    seg.prog = @(job)spm_run_nnseg('run',job);
    seg.vout = @(job)spm_run_nnseg('vout',job);
    seg.help = {'Generate approximate GM/WM/CSF segmentations.','',...
               ['Uses a form of U-net to predict approximate maps of GM, WM snd CSF. '], '', ...
               ['The network is applied at 1 mm isotropic resolution, and is intended to work ' ...
                'on data in any orientation and with a variety of contrasts. '...
                'The network was trained using a linear combination of T1w, T2w and PD-weighted IXI scans '...
                'in log-space, with the objective of reproducing similar segmentations to SPM.'], '', ...
               ['Note that the extremely heavy augmentation means that the accuracy for ' ...
                'subcortical structures is currently poor. The network currently lacks transformers, ' ...
                'which may (possibly) help deal with images of different contrasts.'], '', ...
               ['This functionality is implemented in the Julia language and several libraries, which will be ' ...
                'installed in a toolbox/Spatial/Apps folder when the code is first run.']};

    upd      = cfg_exbranch;
    upd.tag  = 'update';
    upd.name = 'Julia update';
    upd.prog = @(cfg)spm_julia('update');
    upd.help = {'Update the various installed Julia libraries and trained networks.'};

    nnets        = cfg_choice;
    nnets.tag    = 'nnets';
    nnets.name   = 'NNets';
    nnets.values = {dep seg upd};
    nnets.help   = {'A selection of work-in-progress pre-trained neural networks for various spatial processing tasks.'};

end
%==========================================================================

%==========================================================================
function varargout = depth_cfg

    persistent cfg
    if ~isempty(cfg), varargout = {cfg}; return; end

    GM         = cfg_files;
    GM.tag     = 'GM';
    GM.name    = 'GM maps';
    GM.help    = {
    'Select GM maps. When selecting WM maps, the same order of subjects must be specified for each channel and they must be in register (same position, size, voxel dims etc..).'
                 }';
    GM.filter  = 'image';
    GM.ufilter = '^c1.*';
    GM.num     = [1 Inf];
    GM.preview = @(f) spm_check_registration(char(f));

    WM         = cfg_files;
    WM.tag     = 'WM';
    WM.name    = 'WM maps';
    WM.help    = {['Select GM maps. When selecting WM maps, the same order of subjects ' ...
                   'must be specified for each channel and they must be in register '...
                   '(same position, size, voxel dims etc..).'] };
    WM.filter  = 'image';
    WM.ufilter = '^c2.*';
    WM.num     = [1 Inf];
    WM.preview = @(f) spm_check_registration(char(f));

    [cfg,varargout{1}] = deal({GM WM});
end
%==========================================================================

%==========================================================================
function varargout = nnseg_cfg
    persistent cfg
    if ~isempty(cfg), varargout = {cfg}; return; end

    vols         = cfg_files;
    vols.tag     = 'vols';
    vols.name    = 'Image Volumes';
    vols.help    = {'Select the images to be segmented.'};
    vols.filter  = 'image';
    vols.ufilter = '.*';
    vols.num     = [1 Inf];
    vols.preview = @(f) spm_check_registration(char(f));
    [cfg,varargout{1}] = deal({vols});
end
%==========================================================================

%==========================================================================
function out = spm_run_depth(opt, job)
    switch lower(opt)
    case 'vout'
        out(1)            = cfg_dep;
        out(1).sname      = 'Files';
        out(1).src_output = substruct('.','files');
        out(1).tgt_spec   = cfg_findspec({{'filter','image','strtype','e'}});
    case 'run'
        Ptext = '';
        if length(job.GM) ~= length(job.WM)
            error('Incompatible numbers of images.');
        end
        N = min(length(job.GM), length(job.WM));
        for n=1:N
            p = ['["' deblank(job.GM{n}) '","' deblank(job.WM{n}) '"],'];
            Ptext = [Ptext p];
        end
        spm_julia(['do_depth([' Ptext '])'],'SPMmisc');

        names = job.GM(1:N);
        for i=1:length(names)
            [pth,nam,ext,ind] = spm_fileparts(deblank(names(i)));
            names{i} = fullfile(pth,['depth_' nam(3:end) ext]);
        end
        out = struct('files',{names});

    otherwise
        error('Unrecognised option.');
    end
end
%==========================================================================

%==========================================================================
function out = spm_run_nnseg(opt, job)
    switch lower(opt)
    case 'vout'
        out(1)            = cfg_dep;
        out(1).sname      = 'Files';
        out(1).src_output = substruct('.','files');
        out(1).tgt_spec   = cfg_findspec({{'filter','nifti','strtype','e'}});
        out(2)            = cfg_dep;
        out(2).sname      = 'GM';
        out(2).src_output = substruct('.','GM');
        out(2).tgt_spec   = cfg_findspec({{'filter','image','strtype','e'}});
        out(3)            = cfg_dep;
        out(3).sname      = 'WM';
        out(3).src_output = substruct('.','WM');
        out(3).tgt_spec   = cfg_findspec({{'filter','image','strtype','e'}});
        out(4)            = cfg_dep;
        out(4).sname      = 'CSF';
        out(4).src_output = substruct('.','CSF');
        out(4).tgt_spec   = cfg_findspec({{'filter','image','strtype','e'}});

    case 'run'
        Ptext = '';
        N = length(job.vols);
        for n=1:N
            p = ['["' deblank(job.vols{n}) '",],'];
            Ptext = [Ptext p];
        end
        spm_julia(['apply3Dnet([' Ptext '])'],'SPMmisc');

        gm    = job.vols;
        wm    = job.vols;
        csf   = job.vols;
        names = job.vols;
        for i=1:length(names)
            [pth,nam,ext,ind] = spm_fileparts(deblank(names(i)));
            names{i} = fullfile(pth,['c00' nam ext]);
            gm{i}    = fullfile(pth,['c00' nam ext ',1']);
            wm{i}    = fullfile(pth,['c00' nam ext ',2']);
            csf{i}   = fullfile(pth,['c00' nam ext ',3']);
        end
        out = struct('files',{names},'GM',{gm},'WM',{wm},'CSF',{csf});

    otherwise
        error('Unrecognised option.')
    end
end


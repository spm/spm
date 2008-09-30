function spm_surf_node = spm_cfg_spm_surf
% SPM Configuration file
%_______________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% $Id: spm_cfg_spm_surf.m 2273 2008-09-30 21:40:21Z guillaume $

rev = '$Rev: 2273 $';
% ---------------------------------------------------------------------
% data Grey+white matter image
% ---------------------------------------------------------------------
data         = cfg_files;
data.tag     = 'data';
data.name    = 'Grey+white matter image';
data.help    = {'Images to create rendering/surface from (grey and white matter segments).'};
data.filter  = 'image';
data.ufilter = '.*';
data.num     = [1 Inf];
% ---------------------------------------------------------------------
% mode Output
% ---------------------------------------------------------------------
mode         = cfg_menu;
mode.tag     = 'mode';
mode.name    = 'Output';
mode.help    = {''};
mode.labels  = {'Save Rendering'
                'Save Extracted Surface'
                'Save Rendering and Surface'
                'Save Surface as OBJ format'}';
mode.values  = {1 2 3 4};
mode.val     = {3};
% ---------------------------------------------------------------------
% thresh Surface isovalue(s)
% ---------------------------------------------------------------------
thresh         = cfg_entry;
thresh.tag     = 'thresh';
thresh.name    = 'Surface isovalue(s)';
thresh.help    = {'Enter one or more values at which isosurfaces through the input images will be computed.'};
thresh.strtype = 'e';
thresh.val     = {0.5};
thresh.num     = [1 Inf];
% ---------------------------------------------------------------------
% spm_surf Create Rendering/Surface
% ---------------------------------------------------------------------
spm_surf_node         = cfg_exbranch;
spm_surf_node.tag     = 'spm_surf';
spm_surf_node.name    = 'Create Rendering/Surface';
spm_surf_node.val     = { data mode thresh};
spm_surf_node.help    = {''};
spm_surf_node.prog    = @spm_surf;
spm_surf_node.vout    = @vout_surf;
% ---------------------------------------------------------------------
%
% ---------------------------------------------------------------------
function dep = vout_surf(job)
% fail silently, if job.mode or job.thresh can not be evaluated
try
    cdep = 1;
    if any(job.mode==[1 3]),
        dep(cdep)            = cfg_dep;
        dep(cdep).sname      = 'Render .mat File';
        dep(cdep).src_output = substruct('.','rendfile');
        dep(cdep).tgt_spec   = cfg_findspec({{'filter','mat','strtype','e'}});
        cdep = cdep+1;
    end;

    if any(job.mode==[2 3 4]),
        for k=1:numel(job.thresh)
            if any(job.mode==[2 3]),
                dep(cdep)            = cfg_dep;
                dep(cdep).sname      = sprintf('Surf .mat File (thr=%.02f)', ...
                                               job.thresh(k));
                dep(cdep).src_output = substruct('.','surffile', '()',{k});
                dep(cdep).tgt_spec   = cfg_findspec({{'filter','mat','strtype','e'}});
                cdep = cdep+1;
            end;
            if any(job.mode==4),
                dep(cdep)            = cfg_dep;
                dep(cdep).sname      = sprintf('Surf .obj File (thr=%.02f)', ...
                                               job.thresh(k));
                dep(cdep).src_output = substruct('.','objfile', '()',{k});
                dep(cdep).tgt_spec   = cfg_findspec({{'filter','any','strtype','e'}});
                cdep = cdep+1;
            end;
        end;
    end;
catch
    % something failed, no dependencies
    dep = [];
end;

function rendering = spm_cfg_render
% SPM Configuration file for Render
%__________________________________________________________________________
% Copyright (C) 2013-2016 Wellcome Trust Centre for Neuroimaging

% $Id: spm_cfg_render.m 6925 2016-11-09 17:23:40Z guillaume $

%==========================================================================
% Extract
%==========================================================================

%--------------------------------------------------------------------------
% data Data
%--------------------------------------------------------------------------
data         = cfg_files;
data.tag     = 'data';
data.name    = 'Data';
data.help    = {'Images to create rendering/surface from (usually grey and white matter segmentations).'};
data.filter  = 'image';
data.ufilter = '.*';
data.num     = [1 Inf];

%--------------------------------------------------------------------------
% mode Output
%--------------------------------------------------------------------------
mode        = cfg_menu;
mode.tag    = 'mode';
mode.name   = 'Output';
mode.help   = {'Operation mode.'};
mode.labels = {'Save Rendering'
               'Save Extracted Surface'
               'Save Rendering and Surface'}';
mode.values = {1 2 3};
mode.val    = {3};

%--------------------------------------------------------------------------
% thresh Surface isovalue(s)
%--------------------------------------------------------------------------
thresh         = cfg_entry;
thresh.tag     = 'thresh';
thresh.name    = 'Surface isovalue(s)';
thresh.help    = {...
    'Enter one or more values at which isosurfaces through the input images will be computed.'
    'This is only relevant for extracting surfaces, not rendering.'};
thresh.strtype = 'r';
thresh.val     = {0.5};
thresh.num     = [1 Inf];

%--------------------------------------------------------------------------
% extract Extract Surface
%--------------------------------------------------------------------------
extract      = cfg_exbranch;
extract.tag  = 'extract';
extract.name = 'Extract Surface';
extract.val  = {data mode thresh};
extract.help = {'Surface extraction.'};
extract.prog = @spm_surf;
extract.vout = @vout_extract;

%==========================================================================
% Render
%==========================================================================

%--------------------------------------------------------------------------
% spmmat Select SPM.mat
%--------------------------------------------------------------------------
spmmat         = cfg_files;
spmmat.tag     = 'spmmat';
spmmat.name    = 'Select SPM.mat';
spmmat.help    = {'Select the SPM.mat file that contains the design specification.'};
spmmat.filter  = 'mat';
spmmat.ufilter = '^SPM\.mat$';
spmmat.num     = [1 1];

%--------------------------------------------------------------------------
% contrasts Contrast(s)
%--------------------------------------------------------------------------
contrasts         = cfg_entry;
contrasts.tag     = 'contrasts';
contrasts.name    = 'Contrast(s)';
contrasts.help    = {'Index of contrast(s). If more than one number is entered, analyse a conjunction hypothesis.'};
contrasts.strtype = 'n';
contrasts.num     = [1 Inf];

%--------------------------------------------------------------------------
% threshdesc Threshold type
%--------------------------------------------------------------------------
threshdesc        = cfg_menu;
threshdesc.tag    = 'threshdesc';
threshdesc.name   = 'Threshold type';
threshdesc.help   = {''};
threshdesc.labels = {'FWE' 'none' 'FDR'};
threshdesc.values = {'FWE' 'none' 'FDR'};
threshdesc.val    = {'FWE'};

%--------------------------------------------------------------------------
% thresh Threshold
%--------------------------------------------------------------------------
thresh         = cfg_entry;
thresh.tag     = 'thresh';
thresh.name    = 'Threshold';
thresh.help    = {''};
thresh.strtype = 'r';
thresh.num     = [1 1];
thresh.val     = {0.05};

%--------------------------------------------------------------------------
% extent Extent (voxels)
%--------------------------------------------------------------------------
extent         = cfg_entry;
extent.tag     = 'extent';
extent.name    = 'Extent (voxels)';
extent.help    = {''};
extent.strtype = 'w';
extent.num     = [1 1];
extent.val     = {0};

%--------------------------------------------------------------------------
% contrasts Contrast(s)
%--------------------------------------------------------------------------
contrasts1         = cfg_entry;
contrasts1.tag     = 'contrasts';
contrasts1.name    = 'Contrast(s)';
contrasts1.help    = {'Index of contrast(s) for masking - leave empty for no masking.'};
contrasts1.strtype = 'n';
contrasts1.num     = [1 Inf];

%--------------------------------------------------------------------------
% thresh Mask threshold
%--------------------------------------------------------------------------
thresh1         = cfg_entry;
thresh1.tag     = 'thresh';
thresh1.name    = 'Mask threshold';
thresh1.help    = {''};
thresh1.strtype = 'r';
thresh1.num     = [1 1];
thresh1.val     = {0.05};

%--------------------------------------------------------------------------
% mtype Nature of mask
%--------------------------------------------------------------------------
mtype        = cfg_menu;
mtype.tag    = 'mtype';
mtype.name   = 'Nature of mask';
mtype.help   = {''};
mtype.labels = {'Inclusive' 'Exclusive'};
mtype.values = {0 1};

%--------------------------------------------------------------------------
% mask Mask definition
%--------------------------------------------------------------------------
mask      = cfg_branch;
mask.tag  = 'mask';
mask.name = 'Mask definition';
mask.val  = {contrasts1 thresh1 mtype};
mask.help = {''};

%--------------------------------------------------------------------------
% generic Masking
%--------------------------------------------------------------------------
generic1        = cfg_repeat;
generic1.tag    = 'generic';
generic1.name   = 'Masking';
generic1.help   = {''};
generic1.values = {mask};
generic1.num    = [0 1];

%--------------------------------------------------------------------------
% conspec Contrast query
%--------------------------------------------------------------------------
conspec      = cfg_branch;
conspec.tag  = 'conspec';
conspec.name = 'Contrast query';
conspec.val  = {spmmat contrasts threshdesc thresh extent generic1};
conspec.help = {''};

%--------------------------------------------------------------------------
% generic Contrasts
%--------------------------------------------------------------------------
generic        = cfg_repeat;
generic.tag    = 'generic';
generic.name   = 'Contrasts';
generic.help   = {''};
generic.values = {conspec};
generic.num    = [1 3];

%--------------------------------------------------------------------------
% rendfile Render File
%--------------------------------------------------------------------------
rendfile         = cfg_files;
rendfile.tag     = 'rendfile';
rendfile.name    = 'Render File';
rendfile.help    = {'File containing the images to render on to.'};
rendfile.filter  = {'mat','mesh'};
rendfile.ufilter = '.*';
rendfile.num     = [1 1];

%--------------------------------------------------------------------------
% render Display Surface
%--------------------------------------------------------------------------
render      = cfg_exbranch;
render.tag  = 'display';
render.name = 'Display Surface';
render.val  = {rendfile generic};
render.help = {'Surface rendering.'};
render.prog = @run_render;

%==========================================================================
% rendering Rendering
%==========================================================================
rendering        = cfg_choice;
rendering.tag    = 'render';
rendering.name   = 'Rendering';
rendering.help   = {'Rendering utilities.'};
rendering.values = {extract render};


%==========================================================================
function run_render(job)
for i=1:numel(job.conspec)
    xSPM.swd       = spm_file(job.conspec(i).spmmat{1},'fpath');
    xSPM.Ic        = job.conspec(i).contrasts;
    xSPM.u         = job.conspec(i).thresh;
    xSPM.Im        = [];
    if ~isempty(job.conspec(i).mask)
        xSPM.Im    = job.conspec(i).mask.contrasts;
        xSPM.pm    = job.conspec(i).mask.thresh;
        xSPM.Ex    = job.conspec(i).mask.mtype;
    end
    xSPM.thresDesc = job.conspec(i).threshdesc;
    xSPM.k         = job.conspec(i).extent;
    %xSPM.n        = 1; % conjunction 
    xSPM.units     = {'mm' 'mm' 'mm'};
    [SPM, xSPM]    = spm_getSPM(xSPM);
    dat(i) = struct('XYZ', xSPM.XYZ,...
                    't',   xSPM.Z',...
                    'mat', xSPM.M,...
                    'dim', xSPM.DIM);
end
% Force non-interactive mode...
global prevrend
prevrend = struct('rendfile',job.rendfile{1}, 'brt',1, 'col',eye(3));
spm_render(dat,1,job.rendfile{1});


%==========================================================================
function dep = vout_extract(job)

cdep = 1;
if any(job.mode==[1 3])
    dep(cdep)            = cfg_dep;
    dep(cdep).sname      = 'Render .mat File';
    dep(cdep).src_output = substruct('.','rendfile');
    dep(cdep).tgt_spec   = cfg_findspec({{'filter','mat','strtype','e'}});
    cdep = cdep + 1;
end

if any(job.mode==[2 3])
    for k=1:numel(job.thresh)
        dep(cdep)            = cfg_dep;
        dep(cdep).sname      = sprintf('Surface .gii File (thr=%.02f)', ...
            job.thresh(k));
        dep(cdep).src_output = substruct('.','surffile', '()',{k});
        dep(cdep).tgt_spec   = cfg_findspec({{'filter','mesh','strtype','e'}});
        cdep = cdep + 1;
    end
end

function results = spm_cfg_results
% SPM Configuration file for Results
%__________________________________________________________________________
% Copyright (C) 2005-2012 Wellcome Trust Centre for Neuroimaging

% $Id: spm_cfg_results.m 4771 2012-06-20 14:49:22Z volkmar $


% ---------------------------------------------------------------------
% spmmat Select SPM.mat
% ---------------------------------------------------------------------
spmmat         = cfg_files;
spmmat.tag     = 'spmmat';
spmmat.name    = 'Select SPM.mat';
spmmat.help    = {'Select the SPM.mat file that contains the design specification.'};
spmmat.filter  = 'mat';
spmmat.ufilter = '^SPM\.mat$';
spmmat.num     = [1 1];
% ---------------------------------------------------------------------
% titlestr Results Title
% ---------------------------------------------------------------------
titlestr         = cfg_entry;
titlestr.tag     = 'titlestr';
titlestr.name    = 'Results Title';
titlestr.help    = {'Heading on results page - determined automatically if left empty'};
titlestr.val     = {''};
titlestr.strtype = 's';
titlestr.num     = [0 Inf];
% ---------------------------------------------------------------------
% contrasts Contrast(s)
% ---------------------------------------------------------------------
contrasts         = cfg_entry;
contrasts.tag     = 'contrasts';
contrasts.name    = 'Contrast(s)';
contrasts.help    = {
                     'Index of contrast(s). If more than one number is entered, analyse a conjunction hypothesis.'
                     ''
                     'If only one number is entered, and this number is "Inf", then results are printed for all contrasts found in the SPM.mat file.'
}';
contrasts.strtype = 'e';
contrasts.num     = [1 Inf];
% ---------------------------------------------------------------------
% threshdesc Threshold type
% ---------------------------------------------------------------------
threshdesc         = cfg_menu;
threshdesc.tag     = 'threshdesc';
threshdesc.name    = 'Threshold type';
threshdesc.help    = {''};
threshdesc.labels  = {'FWE' 'none' 'FDR'};
threshdesc.values  = {'FWE' 'none' 'FDR'};
threshdesc.val     = {'FWE'};
% ---------------------------------------------------------------------
% thresh Threshold
% ---------------------------------------------------------------------
thresh         = cfg_entry;
thresh.tag     = 'thresh';
thresh.name    = 'Threshold';
thresh.help    = {''};
thresh.strtype = 'e';
thresh.num     = [1 1];
thresh.val     = {0.05};
% ---------------------------------------------------------------------
% extent Extent (voxels)
% ---------------------------------------------------------------------
extent         = cfg_entry;
extent.tag     = 'extent';
extent.name    = 'Extent (voxels)';
extent.help    = {''};
extent.strtype = 'e';
extent.num     = [1 1];
extent.val     = {0};
% ---------------------------------------------------------------------
% contrasts Contrast(s)
% ---------------------------------------------------------------------
contrasts1         = cfg_entry;
contrasts1.tag     = 'contrasts';
contrasts1.name    = 'Contrast(s)';
contrasts1.help    = {'Index of contrast(s) for masking - leave empty for no masking.'};
contrasts1.strtype = 'e';
contrasts1.num     = [1 Inf];
% ---------------------------------------------------------------------
% thresh Mask threshold
% ---------------------------------------------------------------------
thresh1         = cfg_entry;
thresh1.tag     = 'thresh';
thresh1.name    = 'Mask threshold';
thresh1.help    = {''};
thresh1.strtype = 'e';
thresh1.num     = [1 1];
thresh1.val     = {0.05};
% ---------------------------------------------------------------------
% mtype Nature of mask
% ---------------------------------------------------------------------
mtype         = cfg_menu;
mtype.tag     = 'mtype';
mtype.name    = 'Nature of mask';
mtype.help    = {''};
mtype.labels  = {'Inclusive' 'Exclusive'};
mtype.values  = {0 1};
% ---------------------------------------------------------------------
% mask Mask definition
% ---------------------------------------------------------------------
mask         = cfg_branch;
mask.tag     = 'mask';
mask.name    = 'Mask definition';
mask.val     = {contrasts1 thresh1 mtype};
mask.help    = {''};
% ---------------------------------------------------------------------
% generic Masking
% ---------------------------------------------------------------------
generic1         = cfg_repeat;
generic1.tag     = 'generic';
generic1.name    = 'Masking';
generic1.help    = {''};
generic1.values  = {mask};
generic1.num     = [0 1];
% ---------------------------------------------------------------------
% conspec Contrast query
% ---------------------------------------------------------------------
conspec         = cfg_branch;
conspec.tag     = 'conspec';
conspec.name    = 'Contrast query';
conspec.val     = {titlestr contrasts threshdesc thresh extent generic1};
conspec.help    = {''};
% ---------------------------------------------------------------------
% generic Contrasts
% ---------------------------------------------------------------------
generic         = cfg_repeat;
generic.tag     = 'generic';
generic.name    = 'Contrasts';
generic.help    = {''};
generic.values  = {conspec};
generic.num     = [1 Inf];
% ---------------------------------------------------------------------
% units Units
% ---------------------------------------------------------------------
units         = cfg_menu;
units.tag     = 'units';
units.name    = 'Data type';
units.help    = {['Data type. This option is only meaningful for M/EEG '...
    'data. Keep the default ''Volumetric'' for any other kind of data.']};
units.labels  = {'Volumetric (2D/3D)',...
                 'Scalp-Time',...
                 'Scalp-Frequency',...
                 'Time-Frequency',...
                 'Frequency-Frequency'}';
units.values  = { 1 2 3 4 5 };
units.val     = { 1 };
% ---------------------------------------------------------------------
% basename Basename
% ---------------------------------------------------------------------
basename         = cfg_entry;
basename.tag     = 'basename';
basename.name    = 'Basename';
basename.help    = {'Enter basename of output files ''spm?_????_<basename>.ext''.'};
basename.strtype = 's';
basename.num     = [1 Inf];
% ---------------------------------------------------------------------
% type1 Thresholded SPM
% ---------------------------------------------------------------------
type1         = cfg_branch;
type1.tag     = 'tspm';
type1.name    = 'Thresholded SPM';
type1.val     = {basename};
type1.help    = {'Save filtered SPM{.} as an image.'};
% ---------------------------------------------------------------------
% type2 All clusters (binary)
% ---------------------------------------------------------------------
type2         = cfg_branch;
type2.tag     = 'binary';
type2.name    = 'All clusters (binary)';
type2.val     = {basename};
type2.help    = {'Save filetered SPM{.} as a binary image.'};
% ---------------------------------------------------------------------
% type3 All clusters (n-ary)
% ---------------------------------------------------------------------
type3         = cfg_branch;
type3.tag     = 'nary';
type3.name    = 'All clusters (n-ary)';
type3.val     = {basename};
type3.help    = {'Save filtered SPM{.} as an n-ary image.'};
% ---------------------------------------------------------------------
% none None
% ---------------------------------------------------------------------
none          = cfg_const;
none.tag      = 'none';
none.name     = 'None';
none.val      = { 1 };
none.help     = {'Don''t save filtered images.'};
% ---------------------------------------------------------------------
% write Write filtered images
% ---------------------------------------------------------------------
write         = cfg_choice;
write.tag     = 'write';
write.name    = 'Write filtered images';
write.val     = {none};
write.help    = {''};
write.values  = { none type1 type2 type3 };
% ---------------------------------------------------------------------
% print Print results
% ---------------------------------------------------------------------
print         = cfg_menu;
print.tag     = 'print';
print.name    = 'Print results';
print.help    = {''};
print.labels  = {'Yes' 'No'};
print.values  = {true false};
print.val     = {true};
% ---------------------------------------------------------------------
% results Results Report
% ---------------------------------------------------------------------
results          = cfg_exbranch;
results.tag      = 'results';
results.name     = 'Results Report';
results.val      = {spmmat generic units print write};
results.help     = {''};
results.prog     = @spm_run_results;
results.vout     = @vout_results;
results.modality = {'FMRI' 'PET' 'EEG'};


%==========================================================================
function dep = vout_results(job)

dep(1)            = cfg_dep;
dep(1).sname      = 'xSPM Variable';
dep(1).src_output = substruct('.','xSPMvar');
dep(1).tgt_spec   = cfg_findspec({{'strtype','e'}});

dep(2)            = cfg_dep;
dep(2).sname      = 'TabDat Variable';
dep(2).src_output = substruct('.','TabDatvar');
dep(2).tgt_spec   = cfg_findspec({{'strtype','e'}});

if ~isfield(job.write,'none')
    dep(3)            = cfg_dep;
    dep(3).sname      = 'Filtered image';
    dep(3).src_output = substruct('.','filtered');
    dep(3).tgt_spec   = cfg_findspec({{'filter','image','strtype','e'}});
end

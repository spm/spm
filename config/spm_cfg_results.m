function results = spm_cfg_results
% SPM Configuration file for Results Report
%__________________________________________________________________________
% Copyright (C) 2005-2016 Wellcome Trust Centre for Neuroimaging

% $Id: spm_cfg_results.m 6896 2016-10-03 16:53:31Z guillaume $


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
% titlestr Results Title
%--------------------------------------------------------------------------
titlestr         = cfg_entry;
titlestr.tag     = 'titlestr';
titlestr.name    = 'Results Title';
titlestr.help    = {'Heading on results page - determined automatically if left empty'};
titlestr.val     = {''};
titlestr.strtype = 's';
titlestr.num     = [0 Inf];
titlestr.hidden  = true;

%--------------------------------------------------------------------------
% contrasts Contrast(s)
%--------------------------------------------------------------------------
contrasts         = cfg_entry;
contrasts.tag     = 'contrasts';
contrasts.name    = 'Contrast(s)';
contrasts.help    = {
                     'Index of contrast(s). If more than one number is entered, analyse a conjunction hypothesis.'
                     ''
                     'If only one number is entered, and this number is "Inf", then results are printed for all contrasts found in the SPM.mat file.'
}';
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
% conjunction Conjunction Number
%--------------------------------------------------------------------------
conjunction         = cfg_entry;
conjunction.tag     = 'conjunction';
conjunction.name    = 'Conjunction number';
conjunction.help    = {'Conjunction number. Unused if a simple contrast is entered.'
    'For Conjunction Null, enter 1.'
    'For Global Null, enter the number of selected contrasts.'
    'For Intermediate, enter the number of selected contrasts minus the number of effects under the Null.'}';
conjunction.strtype = 'n';
conjunction.num     = [1 1];
conjunction.val     = {1};

%--------------------------------------------------------------------------
% contrasts Contrast(s)
%--------------------------------------------------------------------------
contrasts1         = cfg_entry;
contrasts1.tag     = 'contrasts';
contrasts1.name    = 'Contrast(s)';
contrasts1.help    = {'Index of contrast(s) for masking.'};
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
% mask Mask using contrast
%--------------------------------------------------------------------------
contrast      = cfg_branch;
contrast.tag  = 'contrast';
contrast.name = 'Contrast';
contrast.val  = {contrasts1 thresh1 mtype};
contrast.help = {'Masking using contrast.'};

%--------------------------------------------------------------------------
% name Mask image
%--------------------------------------------------------------------------
name         = cfg_files;
name.tag     = 'name';
name.name    = 'Mask image(s)';
name.help    = {''};
name.filter  = {'image'};
name.ufilter = '.*';
name.num     = [1 Inf];

%--------------------------------------------------------------------------
% image Mask using image
%--------------------------------------------------------------------------
image      = cfg_branch;
image.tag  = 'image';
image.name = 'Image';
image.val  = {name mtype};
image.help = {'Masking using image(s).'};

%--------------------------------------------------------------------------
% none No Masking
%--------------------------------------------------------------------------
none      = cfg_const;
none.tag  = 'none';
none.name = 'None';
none.val  = { 1 };
none.help = {'No masking.'};

%--------------------------------------------------------------------------
% mask Masking
%--------------------------------------------------------------------------
mask        = cfg_choice;
mask.tag    = 'mask';
mask.name   = 'Masking';
mask.help   = {''};
mask.values = {none contrast image};
mask.val    = {none};

%--------------------------------------------------------------------------
% conspec Contrast query
%--------------------------------------------------------------------------
conspec      = cfg_branch;
conspec.tag  = 'conspec';
conspec.name = 'Contrast query';
conspec.val  = {titlestr contrasts threshdesc thresh extent conjunction mask};
conspec.help = {''};

%--------------------------------------------------------------------------
% generic Contrasts
%--------------------------------------------------------------------------
generic        = cfg_repeat;
generic.tag    = 'generic';
generic.name   = 'Contrasts';
generic.help   = {''};
generic.values = {conspec};
generic.num    = [1 Inf];

%--------------------------------------------------------------------------
% units Units
%--------------------------------------------------------------------------
units        = cfg_menu;
units.tag    = 'units';
units.name   = 'Data type';
units.help   = {['Data type. This option is only meaningful for M/EEG '...
    'data. Keep the default ''Volumetric'' for any other kind of data.']};
units.labels = {'Volumetric (2D/3D)',...
                'Scalp-Time',...
                'Scalp-Frequency',...
                'Time-Frequency',...
                'Frequency-Frequency'}';
units.values = { 1 2 3 4 5 };
units.val    = { 1 };

%--------------------------------------------------------------------------
% basename Basename
%--------------------------------------------------------------------------
basename         = cfg_entry;
basename.tag     = 'basename';
basename.name    = 'Basename';
basename.help    = {'Enter basename of output files ''spm?_????_<basename>.ext''.'};
basename.strtype = 's';
basename.num     = [1 Inf];

%--------------------------------------------------------------------------
% nsubj Number of subjects
%--------------------------------------------------------------------------
nsubj         = cfg_entry;
nsubj.tag     = 'nsubj';
nsubj.name    = 'Number of subjects';
nsubj.help    = {'Number of subjects.'};
nsubj.strtype = 'r';
nsubj.num     = [1 1];

%--------------------------------------------------------------------------
% label Label
%--------------------------------------------------------------------------
grplabel         = cfg_entry;
grplabel.tag     = 'label';
grplabel.name    = 'Label';
grplabel.help    = {'Group label.'};
grplabel.strtype = 's';
grplabel.num     = [0 Inf];

%--------------------------------------------------------------------------
% group 
%--------------------------------------------------------------------------
group      = cfg_branch;
group.tag  = 'group';
group.name = 'Group';
group.val  = {nsubj grplabel};
group.help = {['Number of subjects and labels per group. ', ...
    'For a single subject analysis, enter "1" and "single subject".']};

%--------------------------------------------------------------------------
% groups
%--------------------------------------------------------------------------
groups        = cfg_repeat;
groups.tag    = 'groups';
groups.name   = 'Groups';
groups.help   = {['Number of groups. ', ...
    'For a single subject analysis, specify one group.']};
groups.values = {group};
groups.num    = [1 Inf];

%--------------------------------------------------------------------------
% modality Modality
%--------------------------------------------------------------------------
modality        = cfg_menu;
modality.tag    = 'modality';
modality.name   = 'Modality';
modality.help   = {'Modality.'};
modality.labels = {'Anatomical MRI',...
                   'Functional MRI',...
                   'Diffusion MRI',...
                   'PET',...
                   'SPECT',...
                   'EEG',...
                   'MEG'
}';
modality.values = {'AMRI','FMRI','DMRI','PET','SPECT','EEG','MEG'};

%--------------------------------------------------------------------------
% refspace Reference space
%--------------------------------------------------------------------------
refspace        = cfg_menu;
refspace.tag    = 'refspace';
refspace.name   = 'Reference space';
refspace.help   = {['Reference space. For an experiment completed only ',...
    'within SPM, choose one of the first four options.']};
refspace.labels = {'Subject space (no normalisation)',...
                   'Normalised space (using segment)',...
                   'Normalised space (using old segment)',...
                   'Customised space',...
                   'Other normalised MNI space',...
                   'Other normalised Talairach space',...
}';
refspace.values = {'subject','ixi','icbm','custom','mni','talairach'};

%--------------------------------------------------------------------------
% exports
%--------------------------------------------------------------------------
exports{1}      = cfg_branch;
exports{1}.tag  = 'tspm';
exports{1}.name = 'Thresholded SPM';
exports{1}.val  = { basename };
exports{1}.help = {'Save filtered SPM{.} as an image.'};

exports{end+1}    = cfg_branch;
exports{end}.tag  = 'binary';
exports{end}.name = 'All clusters (binary)';
exports{end}.val  = { basename };
exports{end}.help = {'Save filtered SPM{.} as a binary image.'};

exports{end+1}    = cfg_branch;
exports{end}.tag  = 'nary';
exports{end}.name = 'All clusters (n-ary)';
exports{end}.val  = { basename };
exports{end}.help = {'Save filtered SPM{.} as an n-ary image.'};

pf = spm_print('format');
for i=1:numel(pf)
    exports{end+1}    = cfg_const;
    exports{end}.tag  = pf(i).label{1};
    exports{end}.name = pf(i).name;
    exports{end}.val  = { true };
    exports{end}.help = {pf(i).name};
end
exports{end+1}    = cfg_const;
exports{end}.tag  = 'csv';
exports{end}.name = 'CSV file';
exports{end}.val  = { true };
exports{end}.help = {exports{end}.name};
if ispc
    exports{end+1}    = cfg_const;
    exports{end}.tag  = 'xls';
    exports{end}.name = 'Excel spreadsheet file';
    exports{end}.val  = { true };
    exports{end}.help = {exports{end}.name};
end
exports{end+1}    = cfg_branch;
exports{end}.tag  = 'nidm';
exports{end}.name = 'NIDM (Neuroimaging Data Model)';
exports{end}.val  = {modality refspace groups};
exports{end}.help = {exports{end}.name};

%--------------------------------------------------------------------------
% export Export results
%--------------------------------------------------------------------------
export        = cfg_repeat;
export.tag    = 'export';
export.name   = 'Export results';
export.help   = {['Select the export format you want. PostScript (PS) is '...
               'the only format that allows to append figures to the same ' ...
               'file.']};
export.values = exports;
uiprintdef    = spm_get_defaults('ui.print');
for i=1:numel(exports)
    if strcmp(uiprintdef,exports{i}.tag)
        export.val = exports(i);
        break;
    end
end

%--------------------------------------------------------------------------
% results Results Report
%--------------------------------------------------------------------------
results      = cfg_exbranch;
results.tag  = 'results';
results.name = 'Results Report';
results.val  = {spmmat generic units export};
results.help = {''};
results.prog = @spm_run_results;
results.vout = @vout_results;


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

for i=1:numel(job.export)
    fn = char(fieldnames(job.export{i}));
    if ismember(fn,{'tspm','binary','nary'})
        dep(end+1)          = cfg_dep;
        dep(end).sname      = 'Filtered image';
        dep(end).src_output = substruct('.','filtered');
        dep(end).tgt_spec   = cfg_findspec({{'filter','image','strtype','e'}});
        break; % so far, a single dependency pointing to the last exported image
    end
end

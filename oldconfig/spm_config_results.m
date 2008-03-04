function conf = spm_config_results
% Configuration file for results reporting
%_______________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% Volkmar Glauche
% $Id: spm_config_results.m 1185 2008-03-04 16:31:21Z volkmar $


%-------------------------------------------------------------------------

spm.type = 'files';
spm.name = 'Select SPM.mat';
spm.tag  = 'spmmat';
spm.num  = [1 1];
spm.filter  = 'mat';
spm.ufilter = '^SPM\.mat$';
spm.help   = {'Select the SPM.mat file that contains the design specification.'};

print.type   = 'menu';
print.name   = 'Print results';
print.tag    = 'print';
print.labels = {'Yes','No'};
print.values = {1,0};
print.val    = {1};

%-------------------------------------------------------------------------

threshdesc.type   = 'menu';
threshdesc.name   = 'Threshold type';
threshdesc.tag    = 'threshdesc';
threshdesc.labels = {'FWE','FDR','none'};
threshdesc.values = {'FWE','FDR','none'};

thresh.type    = 'entry';
thresh.name    = 'Threshold';
thresh.tag     = 'thresh';
thresh.strtype = 'e';
thresh.num     = [1 1];
thresh.val     = {.05};

extent.type    = 'entry';
extent.name    = 'Extent (voxels)';
extent.tag     = 'extent';
extent.strtype = 'e';
extent.num     = [1 1];
extent.val     = {0};

titlestr.type    = 'entry';
titlestr.name    = 'Results Title';
titlestr.tag     = 'titlestr';
titlestr.strtype = 's';
titlestr.num     = [1 1];
titlestr.help    = {['Heading on results page - determined automatically if' ...
         ' left empty']};
titlestr.val     = {''};

contrasts.type    = 'entry';
contrasts.name    = 'Contrast(s)';
contrasts.tag     = 'contrasts';
contrasts.strtype = 'e';
contrasts.num     = [1 Inf];
contrasts.help    = {['Index of contrast(s). If more than one number is' ...
            ' entered, analyse a conjunction hypothesis.'], ...
            '',...
            ['If only one number is entered, and this number is' ...
            ' "Inf", then results are printed for all contrasts' ...
            ' found in the SPM.mat file.']};

%-------------------------------------------------------------------------

mthresh = thresh;
mthresh.name = 'Mask threshold';

mcons = contrasts;
mcons.help = {'Index of contrast(s) for masking - leave empty for no masking.'};

mtype.type   = 'menu';
mtype.name   = 'Nature of mask';
mtype.tag    = 'mtype';
mtype.labels = {'Inclusive','Exclusive'};
mtype.values = {0,1};

mask.type = 'branch';
mask.name = 'Mask definition';
mask.tag  = 'mask';
mask.val  = {mcons, mthresh, mtype};

masks.type = 'repeat';
masks.name = 'Masking';
masks.tag  = 'masks';
masks.values = {mask};
masks.num  = [0 1];

%-------------------------------------------------------------------------

conspec.type = 'branch';
conspec.name = 'Contrast query';
conspec.tag  = 'conspec';
conspec.val  = {titlestr, contrasts, threshdesc, thresh, extent, masks};

conspecs.type = 'repeat';
conspecs.name = 'Contrasts';
conspecs.tag  = 'conspecs';
conspecs.values = {conspec};
conspecs.num  = [1 Inf];

%-------------------------------------------------------------------------

conf.type   = 'branch';
conf.name   = 'Results Report';
conf.tag    = 'results';
conf.val    = {spm,conspecs,print};
conf.prog   = @run_results;
conf.modality = {'FMRI','PET'};

return;
%=======================================================================

%=======================================================================
function run_results(job)
cspec = job.conspec;
for k = 1:numel(cspec)
    job.conspec=cspec(k);
    if (numel(cspec(k).contrasts) == 1) && isinf(cspec(k).contrasts)
    tmp=load(job.spmmat{1});
    for l=1:numel(tmp.SPM.xCon)
        cspec1(l) = cspec(k);
        cspec1(l).contrasts = l;
    end;
    job1=job;
    job1.print = 1;
    job1.conspec = cspec1;
    run_results(job1);
    else
    [hReg xSPM SPM]=spm_results_ui('Setup',job);
    if job.print
        spm_list('List',xSPM,hReg);
        spm_figure('Print');
    end;
    assignin('base','hReg',hReg);
    assignin('base','xSPM',xSPM);
    assignin('base','SPM',SPM);
    end;
end;

function review = spm_cfg_model_review
% SPM Configuration file for Model Review
%__________________________________________________________________________
% Copyright (C) 2014-2016 Wellcome Trust Centre for Neuroimaging

% Guillaume Flandin
% $Id: spm_cfg_model_review.m 6952 2016-11-25 16:03:13Z guillaume $


%--------------------------------------------------------------------------
% spmmat Select SPM.mat
%--------------------------------------------------------------------------
spmmat         = cfg_files;
spmmat.tag     = 'spmmat';
spmmat.name    = 'Select SPM.mat';
spmmat.help    = {
                  'Select the SPM.mat file that contains the design specification. '
                  'The directory containing this file is known as the input directory.'
}';
spmmat.filter  = 'mat';
spmmat.ufilter = '^SPM\.mat$';
spmmat.num     = [1 1];

%--------------------------------------------------------------------------
% matrix Design Matrix
%--------------------------------------------------------------------------
matrix          = cfg_const;
matrix.tag      = 'matrix';
matrix.name     = 'Design Matrix';
matrix.val      = { 1 };
matrix.help     = {'Review Design Matrix.'};

%--------------------------------------------------------------------------
% factors Files & Factors
%--------------------------------------------------------------------------
factors      = cfg_const;
factors.tag  = 'factors';
factors.name = 'Files & Factors';
factors.val  = { 1 };
factors.help = {'Review Files & Factors. Only available for second-level models.'};


%--------------------------------------------------------------------------
% orth Design Orthogonality
%--------------------------------------------------------------------------
orth          = cfg_const;
orth.tag      = 'orth';
orth.name     = 'Design Orthogonality';
orth.val      = { 1 };
orth.help     = {'Review Design Orthogonality.'};

%--------------------------------------------------------------------------
% sess Session
%--------------------------------------------------------------------------
sess         = cfg_entry;
sess.tag     = 'sess';
sess.name    = 'Session';
sess.help    = {'Index of session.'};
sess.strtype = 'n';
sess.num     = [1 Inf];

%--------------------------------------------------------------------------
% cond Condition
%--------------------------------------------------------------------------
cond         = cfg_entry;
cond.tag     = 'cond';
cond.name    = 'Condition';
cond.help    = {'Index of condition.'};
cond.strtype = 'n';
cond.num     = [1 Inf];

%--------------------------------------------------------------------------
% condition Condition
%--------------------------------------------------------------------------
condition      = cfg_branch;
condition.tag  = 'condition';
condition.name = 'Condition';
condition.val  = {sess cond};
condition.help = {'Review Condition. Only available for fMRI first-level models.'};

%--------------------------------------------------------------------------
% covariance Covariance Structure
%--------------------------------------------------------------------------
covariance      = cfg_const;
covariance.tag  = 'covariance';
covariance.name = 'Covariance Structure';
covariance.val  = { 1 };
covariance.help = {'Review Covariance Structure.'};

%--------------------------------------------------------------------------
% covariates Covariates
%--------------------------------------------------------------------------
covariates      = cfg_const;
covariates.tag  = 'covariates';
covariates.name = 'Covariates';
covariates.val  = { 1 };
covariates.help = {'Review Covariance Structure. Only available for second-level models.'};


%--------------------------------------------------------------------------
% display Display
%--------------------------------------------------------------------------
display         = cfg_choice;
display.tag     = 'display';
display.name    = 'Display';
display.val     = {matrix};
display.help    = {'Select graphical report.'};
display.values  = {matrix orth factors covariates condition covariance};

%--------------------------------------------------------------------------
% print Print results
%--------------------------------------------------------------------------
print        = cfg_menu;
print.tag    = 'print';
print.name   = 'Print results';
print.help   = {['Select the printing format you want. PostScript (PS) is '...
               'the only format that allows to append figures to the same ' ...
               'file.']};
pf           = spm_print('format');
print.labels = {'No'};
print.values = {false};
for i=1:numel(pf)
    print.labels{end+1} = pf(i).name;
    print.values{end+1} = pf(i).label{1};
end
print.def = @(val)spm_get_defaults('ui.print', val{:});

%--------------------------------------------------------------------------
% review Model Review
%--------------------------------------------------------------------------
review          = cfg_exbranch;
review.tag      = 'review';
review.name     = 'Model review';
review.val      = {spmmat display print};
review.help     = {'Review a General linear Model.'};
review.prog     = @spm_run_model_review;
review.modality = {'FMRI' 'PET' 'EEG'};


%==========================================================================
function spm_run_model_review(job)

load(job.spmmat{1});

action = char(fieldnames(job.display));

switch action
    case {'matrix','factors'}
        try
            filenames = reshape(cellstr(SPM.xY.P),size(SPM.xY.VY));
        catch
            filenames = {};
        end
end

switch action
    case 'matrix'
        spm_DesRep('DesMtx',SPM.xX,filenames,SPM.xsDes);
        
    case 'factors'
        if isfield(SPM.xX,'I')
            spm_DesRep('Files&Factors',...
                filenames,...
                SPM.xX.I,SPM.xC,SPM.xX.sF,SPM.xsDes);
        end
        
    case 'orth'
        spm_DesRep('DesOrth',SPM.xX);
        
    case 'covariance'
        spm_DesRep('xVi', SPM.xVi);
        
    case 'covariates'
        if isfield(SPM,'xC') && ~isempty(SPM.xC)
            spm_DesRep('Covs',SPM.xX,SPM.xC);
        end
        
    case 'condition'
        s = job.display.condition.sess;
        c = job.display.condition.cond;
        if s > length(SPM.Sess)
            error('Session not found.');
        end
        if c > length(SPM.Sess(s).Fc)
            error('Condition not found.');
        end
        spm_DesRep('fMRIDesMtx',SPM,s,c);
end

if ~isequal(job.print, false)
    spm_figure('Print','Graphics','',job.print);
end

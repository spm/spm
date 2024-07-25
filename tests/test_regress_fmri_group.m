function tests = test_regress_fmri_group
% Regression tests for second-level SPM for fMRI

% Manual specification of tests to run
tests = functiontests({@setup;
                       @test_regress_onesample_ttest_classical;
                       @test_regress_twosample_ttest_classical;
                       @test_regress_flex_factorial_classical;
                       @test_regress_flex_factorial_bayesian});

% -------------------------------------------------------------------------
function setup(testCase) 

% Start SPM
spm('defaults','fmri');
spm_jobman('initcfg');
spm_get_defaults('cmdline',true);

% -------------------------------------------------------------------------
function test_regress_onesample_ttest_classical(testCase)

xyz_mm = run_onesample_ttest(false);

% Check the global peak hasn't moved (3mm tolerance)
tol = [3 3 3]';
testCase.assertEqual(xyz_mm,[42 -82 -6]','AbsTol', tol);

% -------------------------------------------------------------------------
function test_regress_twosample_ttest_classical(testCase)

xyz_mm = run_twosample_ttest(false);

% Check the global peak hasn't moved (3mm tolerance)
tol = [3 3 3]';
testCase.assertEqual(xyz_mm,[40 -22 52]','AbsTol', tol);

% -------------------------------------------------------------------------
function test_regress_flex_factorial_classical(testCase)

xyz_mm = run_flex_factorial(false);

% Check the global peak hasn't moved (3mm tolerance)
tol = [3 3 3]';
testCase.assertEqual(xyz_mm,[40 -22 52]','AbsTol', tol);

% -------------------------------------------------------------------------
function test_regress_flex_factorial_bayesian(testCase)

xyz_mm = run_flex_factorial(true);

% Check the global peak hasn't moved (3mm tolerance)
tol = [3 3 3]';
testCase.assertEqual(xyz_mm,[40 -22 52]','AbsTol', tol);

% -------------------------------------------------------------------------
function xyz_mm = run_onesample_ttest(bayesian)
% Run a second level one sample t-test on fMRI data and get the global peak

% Locate input folders
base_dir   = fullfile(spm('Dir'),'tests','data','fMRI','group_analysis_semantic','derivatives');
level1_dir = fullfile(base_dir,'first_level');

% Locate / create output folder
level2_dir = fullfile(spm('Dir'),'tests','output','group_analysis_semantic');
if ~exist(level2_dir,'file')
    mkdir(level2_dir);
end

% Ensure it's empty
if exist(fullfile(level2_dir,'SPM.mat'),'file')
    delete(fullfile(level2_dir,'SPM.mat'));
end

% Locate 1st level contrast images
P = cellstr(spm_select('FPListRec',level1_dir,'con_0009.nii'));

% Make GLM specification batch
clear matlabbatch
matlabbatch{1}.spm.stats.factorial_design.dir = cellstr(level2_dir);
matlabbatch{1}.spm.stats.factorial_design.des.t1.scans = P;
matlabbatch{1}.spm.stats.factorial_design.cov = struct('c', {}, 'cname', {}, 'iCFI', {}, 'iCC', {});
matlabbatch{1}.spm.stats.factorial_design.multi_cov = struct('files', {}, 'iCFI', {}, 'iCC', {});
matlabbatch{1}.spm.stats.factorial_design.masking.tm.tm_none = 1;
matlabbatch{1}.spm.stats.factorial_design.masking.im = 1;
matlabbatch{1}.spm.stats.factorial_design.masking.em = {''};
matlabbatch{1}.spm.stats.factorial_design.globalc.g_omit = 1;
matlabbatch{1}.spm.stats.factorial_design.globalm.gmsca.gmsca_no = 1;
matlabbatch{1}.spm.stats.factorial_design.globalm.glonorm = 1;

% Make GLM estimation batch
matlabbatch{2}.spm.stats.fmri_est.spmmat(1) = cfg_dep('Factorial design specification: SPM.mat File', substruct('.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','spmmat'));
matlabbatch{2}.spm.stats.fmri_est.write_residuals = 0;
matlabbatch{2}.spm.stats.fmri_est.method.Classical = 1;

% Run classical estimation
spm_jobman('run',matlabbatch);

% Identify resulting SPM.mat
spm_mat = cellstr(fullfile(level2_dir,'SPM.mat'));

% Run Bayesian estimation if requested
if bayesian    
    clear matlabbatch
    matlabbatch{1}.spm.stats.fmri_est.spmmat = spm_mat;
    matlabbatch{1}.spm.stats.fmri_est.write_residuals = 0;
    matlabbatch{1}.spm.stats.fmri_est.method.Bayesian2 = 1;
    spm_jobman('run',matlabbatch);
end

% Make contrast
clear matlabbatch
matlabbatch{1}.spm.stats.con.spmmat(1) = spm_mat;
matlabbatch{1}.spm.stats.con.consess{1}.tcon.name = 'Activations';
matlabbatch{1}.spm.stats.con.consess{1}.tcon.weights = 1;
matlabbatch{1}.spm.stats.con.consess{1}.tcon.sessrep = 'none';
matlabbatch{1}.spm.stats.con.delete = 1;
spm_jobman('run',matlabbatch);

xyz_mm = get_global_peak(level2_dir,bayesian);

% -------------------------------------------------------------------------
function xyz_mm = run_twosample_ttest(bayesian)
% Run a second level one sample t-test on fMRI data and get the global peak

% Locate input folders
base_dir   = fullfile(spm('Dir'),'tests','data','fMRI','group_analysis_semantic');
level1_dir = fullfile(base_dir,'derivatives','first_level');

% Locate / create output folder
level2_dir = fullfile(spm('Dir'),'tests','output','group_analysis_semantic');
if ~exist(level2_dir,'file')
    mkdir(level2_dir);
end

% Ensure it's empty
if exist(fullfile(level2_dir,'SPM.mat'),'file')
    delete(fullfile(level2_dir,'SPM.mat'));
end

% Load group membership (0=left hand response, 1=right hand response)
T = load(fullfile(base_dir,'participants.mat'));
T = T.participants;
group = T.response_hand;

% Locate 1st level contrast images
P = cellstr(spm_select('FPListRec',level1_dir,'con_0009.nii'));

% Split into groups
P1 = P(group==0); % left-handed
P2 = P(group==1); % right-handed

% Make GLM specification batch
clear matlabbatch
matlabbatch{1}.spm.stats.factorial_design.dir = cellstr(level2_dir);
matlabbatch{1}.spm.stats.factorial_design.des.t2.scans1 = P1;
matlabbatch{1}.spm.stats.factorial_design.des.t2.scans2 = P2;
matlabbatch{1}.spm.stats.factorial_design.des.t2.dept = 0;
matlabbatch{1}.spm.stats.factorial_design.des.t2.variance = 1;
matlabbatch{1}.spm.stats.factorial_design.des.t2.gmsca = 0;
matlabbatch{1}.spm.stats.factorial_design.des.t2.ancova = 0;
matlabbatch{1}.spm.stats.factorial_design.cov = struct('c', {}, 'cname', {}, 'iCFI', {}, 'iCC', {});
matlabbatch{1}.spm.stats.factorial_design.multi_cov = struct('files', {}, 'iCFI', {}, 'iCC', {});
matlabbatch{1}.spm.stats.factorial_design.masking.tm.tm_none = 1;
matlabbatch{1}.spm.stats.factorial_design.masking.im = 1;
matlabbatch{1}.spm.stats.factorial_design.masking.em = {''};
matlabbatch{1}.spm.stats.factorial_design.globalc.g_omit = 1;
matlabbatch{1}.spm.stats.factorial_design.globalm.gmsca.gmsca_no = 1;
matlabbatch{1}.spm.stats.factorial_design.globalm.glonorm = 1;

% Make GLM estimation batch
matlabbatch{2}.spm.stats.fmri_est.spmmat(1) = cfg_dep('Factorial design specification: SPM.mat File', substruct('.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','spmmat'));
matlabbatch{2}.spm.stats.fmri_est.write_residuals = 0;
matlabbatch{2}.spm.stats.fmri_est.method.Classical = 1;

% Run classical estimation
spm_jobman('run',matlabbatch);

% Identify resulting SPM.mat
spm_mat = cellstr(fullfile(level2_dir,'SPM.mat'));

% Run Bayesian estimation if requested
if bayesian    
    clear matlabbatch
    matlabbatch{1}.spm.stats.fmri_est.spmmat = spm_mat;
    matlabbatch{1}.spm.stats.fmri_est.write_residuals = 0;
    matlabbatch{1}.spm.stats.fmri_est.method.Bayesian2 = 1;
    spm_jobman('run',matlabbatch);
end

% Make contrast
clear matlabbatch
matlabbatch{1}.spm.stats.con.spmmat(1) = spm_mat;
matlabbatch{1}.spm.stats.con.consess{1}.tcon.name = 'Left-Right';
matlabbatch{1}.spm.stats.con.consess{1}.tcon.weights = [1 -1];
matlabbatch{1}.spm.stats.con.consess{1}.tcon.sessrep = 'none';
matlabbatch{1}.spm.stats.con.delete = 1;
spm_jobman('run',matlabbatch);

xyz_mm = get_global_peak(level2_dir,bayesian);

% -------------------------------------------------------------------------
function xyz_mm = run_flex_factorial(bayesian)
% Run a second level one sample t-test on fMRI data and get the global peak

% Locate input folders
base_dir   = fullfile(spm('Dir'),'tests','data','fMRI','group_analysis_semantic');
level1_dir = fullfile(base_dir,'derivatives','first_level');

% Locate / create output folder
level2_dir = fullfile(spm('Dir'),'tests','output','group_analysis_semantic');
if ~exist(level2_dir,'file')
    mkdir(level2_dir);
end

% Ensure it's empty
if exist(fullfile(level2_dir,'SPM.mat'),'file')
    delete(fullfile(level2_dir,'SPM.mat'));
end

% Load group membership (0=left hand response, 1=right hand response)
T = load(fullfile(base_dir,'participants.mat'));
T = T.participants;
handedness = T.handedness;
group      = T.response_hand;

% Factor matrix (subject IDs, handedness, group, unused)
n = size(T,1);
IDs = 1:n;

factor_handedness = zeros(n,1);
factor_handedness(handedness==0) = 1; % left
factor_handedness(handedness==1) = 2; % right

factor_group = zeros(n,1);
factor_group(group==0) = 1; % left
factor_group(group==1) = 2; % right

imatrix = [IDs(:) factor_handedness factor_group ones(n,1)];

% Locate 1st level contrast images
P = cellstr(spm_select('FPListRec',level1_dir,'con_0009.nii'));

% Make GLM specification batch
clear matlabbatch
matlabbatch{1}.spm.stats.factorial_design.dir = cellstr(level2_dir);
matlabbatch{1}.spm.stats.factorial_design.des.fblock.fac(1).name = 'Handedness';
matlabbatch{1}.spm.stats.factorial_design.des.fblock.fac(1).dept = 0;
matlabbatch{1}.spm.stats.factorial_design.des.fblock.fac(1).variance = 1;
matlabbatch{1}.spm.stats.factorial_design.des.fblock.fac(1).gmsca = 0;
matlabbatch{1}.spm.stats.factorial_design.des.fblock.fac(1).ancova = 0;
matlabbatch{1}.spm.stats.factorial_design.des.fblock.fac(2).name = 'Response hand';
matlabbatch{1}.spm.stats.factorial_design.des.fblock.fac(2).dept = 0;
matlabbatch{1}.spm.stats.factorial_design.des.fblock.fac(2).variance = 1;
matlabbatch{1}.spm.stats.factorial_design.des.fblock.fac(2).gmsca = 0;
matlabbatch{1}.spm.stats.factorial_design.des.fblock.fac(2).ancova = 0;
matlabbatch{1}.spm.stats.factorial_design.des.fblock.fsuball.specall.scans = P;
matlabbatch{1}.spm.stats.factorial_design.des.fblock.fsuball.specall.imatrix = imatrix;                                                                          
matlabbatch{1}.spm.stats.factorial_design.des.fblock.maininters{1}.inter.fnums = [1;2];
matlabbatch{1}.spm.stats.factorial_design.cov = struct('c', {}, 'cname', {}, 'iCFI', {}, 'iCC', {});
matlabbatch{1}.spm.stats.factorial_design.multi_cov = struct('files', {}, 'iCFI', {}, 'iCC', {});
matlabbatch{1}.spm.stats.factorial_design.masking.tm.tm_none = 1;
matlabbatch{1}.spm.stats.factorial_design.masking.im = 1;
matlabbatch{1}.spm.stats.factorial_design.masking.em = {''};
matlabbatch{1}.spm.stats.factorial_design.globalc.g_omit = 1;
matlabbatch{1}.spm.stats.factorial_design.globalm.gmsca.gmsca_no = 1;
matlabbatch{1}.spm.stats.factorial_design.globalm.glonorm = 1;

% Make GLM estimation batch
matlabbatch{2}.spm.stats.fmri_est.spmmat(1) = cfg_dep('Factorial design specification: SPM.mat File', substruct('.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','spmmat'));
matlabbatch{2}.spm.stats.fmri_est.write_residuals = 0;
matlabbatch{2}.spm.stats.fmri_est.method.Classical = 1;

% Run classical estimation
spm_jobman('run',matlabbatch);

% Identify resulting SPM.mat
spm_mat = cellstr(fullfile(level2_dir,'SPM.mat'));

% Run Bayesian estimation if requested
if bayesian    
    clear matlabbatch
    matlabbatch{1}.spm.stats.fmri_est.spmmat = spm_mat;
    matlabbatch{1}.spm.stats.fmri_est.write_residuals = 0;
    matlabbatch{1}.spm.stats.fmri_est.method.Bayesian2 = 1;
    spm_jobman('run',matlabbatch);
end

% Make contrast
clear matlabbatch
matlabbatch{1}.spm.stats.con.spmmat(1) = spm_mat;
matlabbatch{1}.spm.stats.con.consess{1}.tcon.name = 'Left-Right Response hand';
matlabbatch{1}.spm.stats.con.consess{1}.tcon.weights = [1 -1 1 -1];
matlabbatch{1}.spm.stats.con.consess{1}.tcon.sessrep = 'none';
matlabbatch{1}.spm.stats.con.delete = 1;
spm_jobman('run',matlabbatch);

xyz_mm = get_global_peak(level2_dir,bayesian);

% -------------------------------------------------------------------------
function xyz_mm = get_global_peak(spm_dir,bayesian)
% Gets the global peak from the SPM (if classical) or from the first 
% log-odds image (if Bayesian)
spm_mat = fullfile(spm_dir,'SPM.mat');

if bayesian
    % The results batch doesn't support Bayesian GLMs, so a manual approach
    
    % Read log-odds image
    [Y,XYZmm] = spm_read_vols(spm_vol(fullfile(spm_dir,'spmP_0001.nii')));
        
    % Locate the global peak
    [~,idx] = max(Y(:));
    xyz_mm = XYZmm(:,idx);   
else    
    % Make results batch
    clear matlabbatch    
    matlabbatch{1}.spm.stats.results.spmmat(1) = cellstr(spm_mat);
    matlabbatch{1}.spm.stats.results.conspec.titlestr = '';
    matlabbatch{1}.spm.stats.results.conspec.contrasts = 1;
    matlabbatch{1}.spm.stats.results.conspec.threshdesc = 'FWE';
    matlabbatch{1}.spm.stats.results.conspec.thresh = 0.05;
    matlabbatch{1}.spm.stats.results.conspec.extent = 0;
    matlabbatch{1}.spm.stats.results.conspec.conjunction = 1;
    matlabbatch{1}.spm.stats.results.conspec.mask.none = 1;
    matlabbatch{1}.spm.stats.results.units = 1;
    matlabbatch{1}.spm.stats.results.export = cell(1, 0);
    
    % Run
    out = spm_jobman('run',matlabbatch);

    % Get the global peak from the SPM results table
    xyz_mm = out{1}.TabDatvar.dat{1,end};
end

% -------------------------------------------------------------------------
% function test_timing_twosample_ttest(testCase)
% % Tests the timing of second level Bayesian modelling
% 
% t = [];
% for i = 1:3
%     tic
%     
%     test_twosample_ttest_bayesian(testCase);
%     
%     t(i) = toc;
% end
% disp(t);
% fprintf('Median time: %2.2f\n seconds',median(t));
%
% % Results:
%
% % Run 1. No parallelisation: 
% % median 96.28 seconds [97.1794   96.2817   96.1885]
% 
% % Run 2. With parallelisation:
% % median 24.29 seconds [24.2643   24.2921   24.3753]
% 
% % Run 3. No parallelisation:
% % median 97.29 seconds [98.1365   97.2127   97.2867]
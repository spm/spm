function factorial_design = spm_cfg_factorial_design
% SPM Configuration file for 2nd-level models
%_______________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% Will Penny
% $Id: spm_cfg_factorial_design.m 3900 2010-05-25 16:17:13Z guillaume $

% ---------------------------------------------------------------------
% dir Directory
% ---------------------------------------------------------------------
dir         = cfg_files;
dir.tag     = 'dir';
dir.name    = 'Directory';
dir.help    = {'Select a directory where the SPM.mat file containing the specified design matrix will be written.'};
dir.filter = 'dir';
dir.ufilter = '.*';
dir.num     = [1 1];

% ---------------------------------------------------------------------
% scans Scans
% ---------------------------------------------------------------------
scans         = cfg_files;
scans.tag     = 'scans';
scans.name    = 'Scans';
scans.help    = {'Select the images.  They must all have the same image dimensions, orientation, voxel size etc.'};
scans.filter = 'image';
scans.ufilter = '.*';
scans.num     = [1 Inf];
% ---------------------------------------------------------------------
% t1 One-sample t-test
% ---------------------------------------------------------------------
t1         = cfg_branch;
t1.tag     = 't1';
t1.name    = 'One-sample t-test';
t1.val     = {scans };
t1.help    = {''};

% ---------------------------------------------------------------------
% scans1 Group 1 scans
% ---------------------------------------------------------------------
scans1         = cfg_files;
scans1.tag     = 'scans1';
scans1.name    = 'Group 1 scans';
scans1.help    = {'Select the images from sample 1.  They must all have the same image dimensions, orientation, voxel size etc.'};
scans1.filter = 'image';
scans1.ufilter = '.*';
scans1.num     = [1 Inf];
% ---------------------------------------------------------------------
% scans2 Group 2 scans
% ---------------------------------------------------------------------
scans2         = cfg_files;
scans2.tag     = 'scans2';
scans2.name    = 'Group 2 scans';
scans2.help    = {'Select the images from sample 2.  They must all have the same image dimensions, orientation, voxel size etc.'};
scans2.filter = 'image';
scans2.ufilter = '.*';
scans2.num     = [1 Inf];
% ---------------------------------------------------------------------
% dept Independence
% ---------------------------------------------------------------------
dept         = cfg_menu;
dept.tag     = 'dept';
dept.name    = 'Independence';
dept.help    = {
                'By default, the measurements are assumed to be independent between levels. '
                ''
                'If you change this option to allow for dependencies, this will violate the assumption of sphericity. It would therefore be an example of non-sphericity. One such example would be where you had repeated measurements from the same subjects - it may then be the case that, over subjects, measure 1 is correlated to measure 2. '
                ''
                'Restricted Maximum Likelihood (REML): The ensuing covariance components will be estimated using ReML in spm_spm (assuming the same for all responsive voxels) and used to adjust the statistics and degrees of freedom during inference. By default spm_spm will use weighted least squares to produce Gauss-Markov or Maximum likelihood estimators using the non-sphericity structure specified at this stage. The components will be found in SPM.xVi and enter the estimation procedure exactly as the serial correlations in fMRI models.'
                ''
}';
dept.labels  = {
               'Yes'
               'No'
}';
dept.values  = {0 1};
dept.val     = {0};
% ---------------------------------------------------------------------
% deptn Independence (default is 'No')
% ---------------------------------------------------------------------
deptn         = cfg_menu;
deptn.tag     = 'dept';
deptn.name    = 'Independence';
deptn.help    = {
                'By default, the measurements are assumed to be dependent between levels. '
                ''
                'If you change this option to allow for dependencies, this will violate the assumption of sphericity. It would therefore be an example of non-sphericity. One such example would be where you had repeated measurements from the same subjects - it may then be the case that, over subjects, measure 1 is correlated to measure 2. '
                ''
                'Restricted Maximum Likelihood (REML): The ensuing covariance components will be estimated using ReML in spm_spm (assuming the same for all responsive voxels) and used to adjust the statistics and degrees of freedom during inference. By default spm_spm will use weighted least squares to produce Gauss-Markov or Maximum likelihood estimators using the non-sphericity structure specified at this stage. The components will be found in SPM.xVi and enter the estimation procedure exactly as the serial correlations in fMRI models.'
                ''
}';
deptn.labels  = {
               'Yes'
               'No'
}';
deptn.values  = {0 1};
deptn.val     = {1};

% ---------------------------------------------------------------------
% variance Variance
% ---------------------------------------------------------------------
variance         = cfg_menu;
variance.tag     = 'variance';
variance.name    = 'Variance';
variance.help    = {
                    'By default, the measurements in each level are assumed to have unequal variance. '
                    ''
                    'This violates the assumption of ''sphericity'' and is therefore an example of ''non-sphericity''.'
                    ''
                    'This can occur, for example, in a 2nd-level analysis of variance, one contrast may be scaled differently from another.  Another example would be the comparison of qualitatively different dependent variables (e.g. normals vs. patients).  Different variances (heteroscedasticy) induce different error covariance components that are estimated using restricted maximum likelihood (see below).'
                    ''
                    'Restricted Maximum Likelihood (REML): The ensuing covariance components will be estimated using ReML in spm_spm (assuming the same for all responsive voxels) and used to adjust the statistics and degrees of freedom during inference. By default spm_spm will use weighted least squares to produce Gauss-Markov or Maximum likelihood estimators using the non-sphericity structure specified at this stage. The components will be found in SPM.xVi and enter the estimation procedure exactly as the serial correlations in fMRI models.'
                    ''
}';
variance.labels = {
                   'Equal'
                   'Unequal'
}';
variance.values = {0 1};
variance.val    = {1};
% ---------------------------------------------------------------------
% gmsca Grand mean scaling
% ---------------------------------------------------------------------
gmsca         = cfg_menu;
gmsca.tag     = 'gmsca';
gmsca.name    = 'Grand mean scaling';
gmsca.help    = {
                 'This option is only used for PET data.'
                 ''
                 'Selecting YES will specify ''grand mean scaling by factor'' which could be eg. ''grand mean scaling by subject'' if the factor is ''subject''. '
                 ''
                 'Since differences between subjects may be due to gain and sensitivity effects, AnCova by subject could be combined with "grand mean scaling by subject" to obtain a combination of between subject proportional scaling and within subject AnCova. '
                 ''
}';
gmsca.labels = {
                'No'
                'Yes'
}';
gmsca.values = {0 1};
gmsca.val    = {0};
% ---------------------------------------------------------------------
% ancova ANCOVA
% ---------------------------------------------------------------------
ancova         = cfg_menu;
ancova.tag     = 'ancova';
ancova.name    = 'ANCOVA';
ancova.help    = {
                  'This option is only used for PET data.'
                  ''
                  'Selecting YES will specify ''ANCOVA-by-factor'' regressors. This includes eg. ''Ancova by subject'' or ''Ancova by effect''. These options allow eg. different subjects to have different relationships between local and global measurements. '
                  ''
}';
ancova.labels = {
                 'No'
                 'Yes'
}';
ancova.values = {0 1};
ancova.val    = {0};
% ---------------------------------------------------------------------
% t2 Two-sample t-test
% ---------------------------------------------------------------------
t2         = cfg_branch;
t2.tag     = 't2';
t2.name    = 'Two-sample t-test';
t2.val     = {scans1 scans2 dept variance gmsca ancova };
t2.help    = {''};

% ---------------------------------------------------------------------
% scans Scans [1,2]
% ---------------------------------------------------------------------
scans         = cfg_files;
scans.tag     = 'scans';
scans.name    = 'Scans [1,2]';
scans.help    = {'Select the pair of images. '};
scans.filter = 'image';
scans.ufilter = '.*';
scans.num     = [2 2];
% ---------------------------------------------------------------------
% pair Pair
% ---------------------------------------------------------------------
pair         = cfg_branch;
pair.tag     = 'pair';
pair.name    = 'Pair';
pair.val     = {scans };
pair.help    = {'Add a new pair of scans to your experimental design'};
% ---------------------------------------------------------------------
% generic Pairs
% ---------------------------------------------------------------------
generic         = cfg_repeat;
generic.tag     = 'generic';
generic.name    = 'Pairs';
generic.help    = {''};
generic.values  = {pair};
generic.num     = [1 Inf];
% ---------------------------------------------------------------------
% pt Paired t-test
% ---------------------------------------------------------------------
pt         = cfg_branch;
pt.tag     = 'pt';
pt.name    = 'Paired t-test';
pt.val     = {generic gmsca ancova};
pt.help    = {''};

% ---------------------------------------------------------------------
% scans Scans
% ---------------------------------------------------------------------
scans         = cfg_files;
scans.tag     = 'scans';
scans.name    = 'Scans';
scans.help    = {'Select the images.  They must all have the same image dimensions, orientation, voxel size etc.'};
scans.filter = 'image';
scans.ufilter = '.*';
scans.num     = [1 Inf];
% ---------------------------------------------------------------------
% c Vector
% ---------------------------------------------------------------------
c         = cfg_entry;
c.tag     = 'c';
c.name    = 'Vector';
c.help    = {'Vector of covariate values'};
c.strtype = 'e';
c.num     = [Inf 1];
% ---------------------------------------------------------------------
% cname Name
% ---------------------------------------------------------------------
cname         = cfg_entry;
cname.tag     = 'cname';
cname.name    = 'Name';
cname.help    = {'Name of covariate'};
cname.strtype = 's';
cname.num     = [1 Inf];
% ---------------------------------------------------------------------
% iCC Centering
% ---------------------------------------------------------------------
iCC         = cfg_menu;
iCC.tag     = 'iCC';
iCC.name    = 'Centering';
iCC.help    = {''};
iCC.labels = {
              'Overall mean'
              'No centering'
}';
iCC.values = {1 5};
iCC.val    = {1};
% ---------------------------------------------------------------------
% mcov Covariate
% ---------------------------------------------------------------------
mcov         = cfg_branch;
mcov.tag     = 'mcov';
mcov.name    = 'Covariate';
mcov.val     = {c cname iCC };
mcov.help    = {'Add a new covariate to your experimental design'};
% ---------------------------------------------------------------------
% generic Covariates
% ---------------------------------------------------------------------
generic         = cfg_repeat;
generic.tag     = 'generic';
generic.name    = 'Covariates';
generic.help    = {'Covariates'};
generic.values  = {mcov };
generic.num     = [0 Inf];
% ---------------------------------------------------------------------
% incint Intercept
% ---------------------------------------------------------------------
incint = cfg_menu;
incint.tag = 'incint';
incint.name = 'Intercept';
incint.help = {['By default, an intercept is always added to the model. If the ',...
    'covariates supplied by the user include a constant effect, the ',...
    'intercept may be omitted.']};
incint.labels = {'Include Intercept','Omit Intercept'};
incint.values = {1,0};
incint.val    = {1};
% ---------------------------------------------------------------------
% mreg Multiple regression
% ---------------------------------------------------------------------
mreg         = cfg_branch;
mreg.tag     = 'mreg';
mreg.name    = 'Multiple regression';
mreg.val     = {scans generic incint};
mreg.help    = {''};
% ---------------------------------------------------------------------
% name Name
% ---------------------------------------------------------------------
name         = cfg_entry;
name.tag     = 'name';
name.name    = 'Name';
name.help    = {'Name of factor, eg. ''Repetition'' '};
name.strtype = 's';
name.num     = [1 Inf];
% ---------------------------------------------------------------------
% levels Levels
% ---------------------------------------------------------------------
levels         = cfg_entry;
levels.tag     = 'levels';
levels.name    = 'Levels';
levels.help    = {'Enter number of levels for this factor, eg. 2'};
levels.strtype = 'e';
levels.num     = [Inf 1];
% ---------------------------------------------------------------------
% fact Factor
% ---------------------------------------------------------------------
fact         = cfg_branch;
fact.tag     = 'fact';
fact.name    = 'Factor';
fact.val     = {name levels dept variance gmsca ancova };
fact.help    = {'Add a new factor to your experimental design'};
% ---------------------------------------------------------------------
% generic Factors
% ---------------------------------------------------------------------
generic         = cfg_repeat;
generic.tag     = 'generic';
generic.name    = 'Factors';
generic.help    = {
                   'Specify your design a factor at a time. '
                   ''
}';
generic.values  = {fact };
generic.num     = [1 Inf];
% ---------------------------------------------------------------------
% levels Levels
% ---------------------------------------------------------------------
levels         = cfg_entry;
levels.tag     = 'levels';
levels.name    = 'Levels';
levels.help    = {
                  'Enter a vector or scalar that specifies which cell in the factorial design these images belong to. The length of this vector should correspond to the number of factors in the design'
                  ''
                  'For example, length 2 vectors should be used for two-factor designs eg. the vector [2 3] specifies the cell corresponding to the 2nd-level of the first factor and the 3rd level of the 2nd factor.'
                  ''
}';
levels.strtype = 'e';
levels.num     = [Inf 1];
% ---------------------------------------------------------------------
% scans Scans
% ---------------------------------------------------------------------
scans         = cfg_files;
scans.tag     = 'scans';
scans.name    = 'Scans';
scans.help    = {'Select the images for this cell.  They must all have the same image dimensions, orientation, voxel size etc.'};
scans.filter = 'image';
scans.ufilter = '.*';
scans.num     = [1 Inf];
% ---------------------------------------------------------------------
% icell Cell
% ---------------------------------------------------------------------
icell         = cfg_branch;
icell.tag     = 'icell';
icell.name    = 'Cell';
icell.val     = {levels scans };
icell.help    = {'Enter data for a cell in your design'};
% ---------------------------------------------------------------------
% scell Cell
% ---------------------------------------------------------------------
scell         = cfg_branch;
scell.tag     = 'icell';
scell.name    = 'Cell';
scell.val     = {scans };
scell.help    = {'Enter data for a cell in your design'};
% ---------------------------------------------------------------------
% generic Specify cells
% ---------------------------------------------------------------------
generic1         = cfg_repeat;
generic1.tag     = 'generic';
generic1.name    = 'Specify cells';
generic1.help    = {
                    'Enter the scans a cell at a time'
                    ''
}';
generic1.values  = {icell };
generic1.num     = [1 Inf];
% ---------------------------------------------------------------------
% generic Specify cells
% ---------------------------------------------------------------------
generic2         = cfg_repeat;
generic2.tag     = 'generic';
generic2.name    = 'Specify cells';
generic2.help    = {
                    'Enter the scans a cell at a time'
                    ''
}';
generic2.values  = {scell };
generic2.num     = [1 Inf];
% ---------------------------------------------------------------------
% anova ANOVA 
% ---------------------------------------------------------------------
anova         = cfg_branch;
anova.tag     = 'anova';
anova.name    = 'One-way ANOVA';
anova.val     = {generic2 dept variance gmsca ancova};
anova.help    = {
              'One-way Analysis of Variance (ANOVA)'
}';
% ---------------------------------------------------------------------
% fd Full factorial
% ---------------------------------------------------------------------
fd         = cfg_branch;
fd.tag     = 'fd';
fd.name    = 'Full factorial';
fd.val     = {generic generic1 };
fd.help    = {
              'This option is best used when you wish to test for all main effects and interactions in one-way, two-way or three-way ANOVAs. Design specification proceeds in 2 stages. Firstly, by creating new factors and specifying the number of levels and name for each. Nonsphericity, ANOVA-by-factor and scaling options can also be specified at this stage. Secondly, scans are assigned separately to each cell. This accomodates unbalanced designs.'
              ''
              'For example, if you wish to test for a main effect in the population from which your subjects are drawn and have modelled that effect at the first level using K basis functions (eg. K=3 informed basis functions) you can use a one-way ANOVA with K-levels. Create a single factor with K levels and then assign the data to each cell eg. canonical, temporal derivative and dispersion derivative cells, where each cell is assigned scans from multiple subjects.'
              ''
              'SPM will also automatically generate the contrasts necessary to test for all main effects and interactions. '
              ''
}';
% ---------------------------------------------------------------------
% name Name
% ---------------------------------------------------------------------
name         = cfg_entry;
name.tag     = 'name';
name.name    = 'Name';
name.help    = {'Name of factor, eg. ''Repetition'' '};
name.strtype = 's';
name.num     = [1 Inf];
% ---------------------------------------------------------------------
% fac Factor
% ---------------------------------------------------------------------
fac         = cfg_branch;
fac.tag     = 'fac';
fac.name    = 'Factor';
fac.val     = {name dept variance gmsca ancova };
fac.help    = {
               'Add a new factor to your design.'
               ''
               'If you are using the ''Subjects'' option to specify your scans and conditions, you may wish to make use of the following facility. There are two reserved words for the names of factors. These are ''subject'' and ''repl'' (standing for replication). If you use these factor names then SPM can automatically create replication and/or subject factors without you having to type in an extra entry in the condition vector.'
               ''
               'For example, if you wish to model Subject and Task effects (two factors), under Subjects->Subject->Conditions you can type in simply [1 2 1 2] to specify eg. just the ''Task'' factor level. You do not need to eg. for the 4th subject enter the matrix [1 4; 2 4; 1 4; 2 4]. '
               ''
}';
% ---------------------------------------------------------------------
% generic Factors
% ---------------------------------------------------------------------
generic         = cfg_repeat;
generic.tag     = 'generic';
generic.name    = 'Factors';
generic.help    = {
                   'Specify your design a factor at a time.'
                   ''
}';
generic.values  = {fac };
generic.num     = [1 Inf];
% ---------------------------------------------------------------------
% scans Scans
% ---------------------------------------------------------------------
scans         = cfg_files;
scans.tag     = 'scans';
scans.name    = 'Scans';
scans.help    = {'Select the images to be analysed.  They must all have the same image dimensions, orientation, voxel size etc.'};
scans.filter = 'image';
scans.ufilter = '.*';
scans.num     = [1 Inf];
% ---------------------------------------------------------------------
% conds Conditions
% ---------------------------------------------------------------------
conds         = cfg_entry;
conds.tag     = 'conds';
conds.name    = 'Conditions';
conds.help    = {''};
conds.strtype = 'e';
conds.num     = [Inf Inf];
% ---------------------------------------------------------------------
% fsubject Subject
% ---------------------------------------------------------------------
fsubject         = cfg_branch;
fsubject.tag     = 'fsubject';
fsubject.name    = 'Subject';
fsubject.val     = {scans conds };
fsubject.help    = {'Enter data and conditions for a new subject'};
% ---------------------------------------------------------------------
% generic Subjects
% ---------------------------------------------------------------------
generic1         = cfg_repeat;
generic1.tag     = 'generic';
generic1.name    = 'Subjects';
generic1.help    = {''};
generic1.values  = {fsubject };
generic1.num     = [1 Inf];
% ---------------------------------------------------------------------
% scans Scans
% ---------------------------------------------------------------------
scans         = cfg_files;
scans.tag     = 'scans';
scans.name    = 'Scans';
scans.help    = {'Select the images to be analysed.  They must all have the same image dimensions, orientation, voxel size etc.'};
scans.filter = 'image';
scans.ufilter = '.*';
scans.num     = [1 Inf];
% ---------------------------------------------------------------------
% imatrix Factor matrix
% ---------------------------------------------------------------------
imatrix         = cfg_entry;
imatrix.tag     = 'imatrix';
imatrix.name    = 'Factor matrix';
imatrix.help    = {'Specify factor/level matrix as a nscan-by-4 matrix. Note that the first column of I is reserved for the internal replication factor and must not be used for experimental factors.'};
imatrix.strtype = 'e';
imatrix.num     = [Inf Inf];
% ---------------------------------------------------------------------
% specall Specify all
% ---------------------------------------------------------------------
specall         = cfg_branch;
specall.tag     = 'specall';
specall.name    = 'Specify all';
specall.val     = {scans imatrix };
specall.help    = {
                   'Specify (i) all scans in one go and (ii) all conditions using a factor matrix, I. This option is for ''power users''. The matrix I must have four columns and as as many rows as scans. It has the same format as SPM''s internal variable SPM.xX.I. '
                   ''
                   'The first column of I denotes the replication number and entries in the other columns denote the levels of each experimental factor.'
                   ''
                   'So, for eg. a two-factor design the first column denotes the replication number and columns two and three have entries like 2 3 denoting the 2nd level of the first factor and 3rd level of the second factor. The 4th column in I would contain all 1s.'
}';
% ---------------------------------------------------------------------
% fsuball Specify Subjects or all Scans & Factors
% ---------------------------------------------------------------------
fsuball         = cfg_choice;
fsuball.tag     = 'fsuball';
fsuball.name    = 'Specify Subjects or all Scans & Factors';
fsuball.val     = {generic1 };
fsuball.help    = {''};
fsuball.values  = {generic1 specall };
% ---------------------------------------------------------------------
% fnum Factor number
% ---------------------------------------------------------------------
fnum         = cfg_entry;
fnum.tag     = 'fnum';
fnum.name    = 'Factor number';
fnum.help    = {'Enter the number of the factor.'};
fnum.strtype = 'e';
fnum.num     = [1 1];
% ---------------------------------------------------------------------
% fmain Main effect
% ---------------------------------------------------------------------
fmain         = cfg_branch;
fmain.tag     = 'fmain';
fmain.name    = 'Main effect';
fmain.val     = {fnum };
fmain.help    = {'Add a main effect to your design matrix'};
% ---------------------------------------------------------------------
% fnums Factor numbers
% ---------------------------------------------------------------------
fnums         = cfg_entry;
fnums.tag     = 'fnums';
fnums.name    = 'Factor numbers';
fnums.help    = {'Enter the numbers of the factors of this (two-way) interaction.'};
fnums.strtype = 'e';
fnums.num     = [2 1];
% ---------------------------------------------------------------------
% inter Interaction
% ---------------------------------------------------------------------
inter         = cfg_branch;
inter.tag     = 'inter';
inter.name    = 'Interaction';
inter.val     = {fnums };
inter.help    = {'Add an interaction to your design matrix'};
% ---------------------------------------------------------------------
% maininters Main effects & Interactions
% ---------------------------------------------------------------------
maininters         = cfg_repeat;
maininters.tag     = 'maininters';
maininters.name    = 'Main effects & Interactions';
maininters.help    = {''};
maininters.values  = {fmain inter };
maininters.num     = [1 Inf];
% ---------------------------------------------------------------------
% anovaw ANOVA within subject
% ---------------------------------------------------------------------
anovaw         = cfg_branch;
anovaw.tag     = 'anovaw';
anovaw.name    = 'One-way ANOVA - within subject';
anovaw.val     = {generic1 deptn variance gmsca ancova};
anovaw.help    = {
              'One-way Analysis of Variance (ANOVA) - within subject'
}';
% ---------------------------------------------------------------------
% fblock Flexible factorial
% ---------------------------------------------------------------------
fblock         = cfg_branch;
fblock.tag     = 'fblock';
fblock.name    = 'Flexible factorial';
fblock.val     = {generic fsuball maininters };
fblock.help    = {
                  'Create a design matrix a block at a time by specifying which main effects and interactions you wish to be included.'
                  ''
                  'This option is best used for one-way, two-way or three-way ANOVAs but where you do not wish to test for all possible main effects and interactions. This is perhaps most useful for PET where there is usually not enough data to test for all possible effects. Or for 3-way ANOVAs where you do not wish to test for all of the two-way interactions. A typical example here would be a group-by-drug-by-task analysis where, perhaps, only (i) group-by-drug or (ii) group-by-task interactions are of interest. In this case it is only necessary to have two-blocks in the design matrix - one for each interaction. The three-way interaction can then be tested for using a contrast that computes the difference between (i) and (ii).'
                  ''
                  'Design specification then proceeds in 3 stages. Firstly, factors are created and names specified for each. Nonsphericity, ANOVA-by-factor and scaling options can also be specified at this stage.'
                  ''
                  'Secondly, a list of scans is produced along with a factor matrix, I. This is an nscan x 4 matrix of factor level indicators (see xX.I below). The first factor must be ''replication'' but the other factors can be anything. Specification of I and the scan list can be achieved in one of two ways (a) the ''Specify All'' option allows I to be typed in at the user interface or (more likely) loaded in from the matlab workspace. All of the scans are then selected in one go. (b) the ''Subjects'' option allows you to enter scans a subject at a time. The corresponding experimental conditions (ie. levels of factors) are entered at the same time. SPM will then create the factor matrix I. This style of interface is similar to that available in SPM2.'
                  ''
                  'Thirdly, the design matrix is built up a block at a time. Each block can be a main effect or a (two-way) interaction. '
                  ''
}';
% ---------------------------------------------------------------------
% des Design
% ---------------------------------------------------------------------
des         = cfg_choice;
des.tag     = 'des';
des.name    = 'Design';
des.val     = {t1 };
des.help    = {''};
des.values  = {t1 t2 pt mreg anova anovaw fd fblock };
% ---------------------------------------------------------------------
% c Vector
% ---------------------------------------------------------------------
c         = cfg_entry;
c.tag     = 'c';
c.name    = 'Vector';
c.help    = {
             'Vector of covariate values.'
             'Enter the covariate values ''''per subject'''' (i.e. all for subject 1, then all for subject 2, etc). Importantly, the ordering of the cells of a factorial design has to be the same for all subjects in order to be consistent with the ordering of the covariate values.'
}';
c.strtype = 'e';
c.num     = [Inf 1];
% ---------------------------------------------------------------------
% cname Name
% ---------------------------------------------------------------------
cname         = cfg_entry;
cname.tag     = 'cname';
cname.name    = 'Name';
cname.help    = {'Name of covariate'};
cname.strtype = 's';
cname.num     = [1 Inf];
% ---------------------------------------------------------------------
% iCFI Interactions
% ---------------------------------------------------------------------
iCFI         = cfg_menu;
iCFI.tag     = 'iCFI';
iCFI.name    = 'Interactions';
iCFI.help    = {
                'For each covariate you have defined, there is an opportunity to create an additional regressor that is the interaction between the covariate and a chosen experimental factor. '
                ''
}';
iCFI.labels = {
               'None'
               'With Factor 1'
               'With Factor 2'
               'With Factor 3'
}';
iCFI.values = {1 2 3 4};
iCFI.val    = {1};
% ---------------------------------------------------------------------
% iCC Centering
% ---------------------------------------------------------------------
iCC         = cfg_menu;
iCC.tag     = 'iCC';
iCC.name    = 'Centering';
iCC.help    = {
               'The appropriate centering option is usually the one that corresponds to the interaction chosen, and ensures that main effects of the interacting factor aren''t affected by the covariate. You are advised to choose this option, unless you have other modelling considerations. '
               ''
}';
iCC.labels = {
              'Overall mean'
              'Factor 1 mean'
              'Factor 2 mean'
              'Factor 3 mean'
              'No centering'
              'User specified value'
              'As implied by ANCOVA'
              'GM'
}';
iCC.values = {1 2 3 4 5 6 7 8};
iCC.val    = {1};
% ---------------------------------------------------------------------
% cov Covariate
% ---------------------------------------------------------------------
cov         = cfg_branch;
cov.tag     = 'cov';
cov.name    = 'Covariate';
cov.val     = {c cname iCFI iCC };
cov.help    = {'Add a new covariate to your experimental design'};
% ---------------------------------------------------------------------
% generic Covariates
% ---------------------------------------------------------------------
generic         = cfg_repeat;
generic.tag     = 'generic';
generic.name    = 'Covariates';
generic.help    = {
                   'This option allows for the specification of covariates and nuisance variables. Unlike SPM94/5/6, where the design was partitioned into effects of interest and nuisance effects for the computation of adjusted data and the F-statistic (which was used to thresh out voxels where there appeared to be no effects of interest), SPM does not partition the design in this way anymore. The only remaining distinction between effects of interest (including covariates) and nuisance effects is their location in the design matrix, which we have retained for continuity.  Pre-specified design matrix partitions can be entered. '
                   ''
}';
generic.values  = {cov };
generic.num     = [0 Inf];
% ---------------------------------------------------------------------
% tm_none None
% ---------------------------------------------------------------------
tm_none         = cfg_const;
tm_none.tag     = 'tm_none';
tm_none.name    = 'None';
tm_none.val     = {1};
tm_none.help    = {'No threshold masking'};
% ---------------------------------------------------------------------
% athresh Threshold
% ---------------------------------------------------------------------
athresh         = cfg_entry;
athresh.tag     = 'athresh';
athresh.name    = 'Threshold';
athresh.help    = {
                   'Enter the absolute value of the threshold.'
                   ''
}';
athresh.strtype = 'e';
athresh.num     = [1 1];
athresh.val     = {100};
% ---------------------------------------------------------------------
% tma Absolute
% ---------------------------------------------------------------------
tma         = cfg_branch;
tma.tag     = 'tma';
tma.name    = 'Absolute';
tma.val     = {athresh };
tma.help    = {
               'Images are thresholded at a given value and only voxels at which all images exceed the threshold are included. '
               ''
               'This option allows you to specify the absolute value of the threshold.'
               ''
}';
% ---------------------------------------------------------------------
% rthresh Threshold
% ---------------------------------------------------------------------
rthresh         = cfg_entry;
rthresh.tag     = 'rthresh';
rthresh.name    = 'Threshold';
rthresh.help    = {
                   'Enter the threshold as a proportion of the global value'
                   ''
}';
rthresh.strtype = 'e';
rthresh.num     = [1 1];
rthresh.val     = {.8};
% ---------------------------------------------------------------------
% tmr Relative
% ---------------------------------------------------------------------
tmr         = cfg_branch;
tmr.tag     = 'tmr';
tmr.name    = 'Relative';
tmr.val     = {rthresh };
tmr.help    = {
               'Images are thresholded at a given value and only voxels at which all images exceed the threshold are included. '
               ''
               'This option allows you to specify the value of the threshold as a proportion of the global value. '
               ''
}';
% ---------------------------------------------------------------------
% tm Threshold masking
% ---------------------------------------------------------------------
tm         = cfg_choice;
tm.tag     = 'tm';
tm.name    = 'Threshold masking';
tm.val     = {tm_none };
tm.help    = {
              'Images are thresholded at a given value and only voxels at which all images exceed the threshold are included. '
              ''
}';
tm.values  = {tm_none tma tmr };
% ---------------------------------------------------------------------
% im Implicit Mask
% ---------------------------------------------------------------------
im         = cfg_menu;
im.tag     = 'im';
im.name    = 'Implicit Mask';
im.help    = {
              'An "implicit mask" is a mask implied by a particular voxel value. Voxels with this mask value are excluded from the analysis. '
              ''
              'For image data-types with a representation of NaN (see spm_type.m), NaN''s is the implicit mask value, (and NaN''s are always masked out). '
              ''
              'For image data-types without a representation of NaN, zero is the mask value, and the user can choose whether zero voxels should be masked out or not.'
              ''
              'By default, an implicit mask is used. '
              ''
}';
im.labels = {
             'Yes'
             'No'
}';
im.values = {1 0};
im.val    = {1};
% ---------------------------------------------------------------------
% em Explicit Mask
% ---------------------------------------------------------------------
em         = cfg_files;
em.tag     = 'em';
em.name    = 'Explicit Mask';
em.val     = {{''}};
em.help    = {
              'Explicit masks are other images containing (implicit) masks that are to be applied to the current analysis.'
              ''
              'All voxels with value NaN (for image data-types with a representation of NaN), or zero (for other data types) are excluded from the analysis. '
              ''
              'Explicit mask images can have any orientation and voxel/image size. Nearest neighbour interpolation of a mask image is used if the voxel centers of the input images do not coincide with that of the mask image.'
              ''
}';
em.filter = 'image';
em.ufilter = '.*';
em.num     = [0 1];
% ---------------------------------------------------------------------
% masking Masking
% ---------------------------------------------------------------------
masking         = cfg_branch;
masking.tag     = 'masking';
masking.name    = 'Masking';
masking.val     = {tm im em };
masking.help    = {
                   'The mask specifies the voxels within the image volume which are to be assessed. SPM supports three methods of masking (1) Threshold, (2) Implicit and (3) Explicit. The volume analysed is the intersection of all masks.'
                   ''
}';
% ---------------------------------------------------------------------
% g_omit Omit
% ---------------------------------------------------------------------
g_omit         = cfg_const;
g_omit.tag     = 'g_omit';
g_omit.name    = 'Omit';
g_omit.val     = {1};
g_omit.help    = {'Omit'};
% ---------------------------------------------------------------------
% global_uval Global values
% ---------------------------------------------------------------------
global_uval         = cfg_entry;
global_uval.tag     = 'global_uval';
global_uval.name    = 'Global values';
global_uval.help    = {
                       'Enter the vector of global values'
                       ''
}';
global_uval.strtype = 'e';
global_uval.num     = [Inf 1];
% ---------------------------------------------------------------------
% g_user User
% ---------------------------------------------------------------------
g_user         = cfg_branch;
g_user.tag     = 'g_user';
g_user.name    = 'User';
g_user.val     = {global_uval };
g_user.help    = {
                  'User defined  global effects (enter your own '
                  'vector of global values)'
}';
% ---------------------------------------------------------------------
% g_mean Mean
% ---------------------------------------------------------------------
g_mean         = cfg_const;
g_mean.tag     = 'g_mean';
g_mean.name    = 'Mean';
g_mean.val     = {1};
g_mean.help    = {
                  'SPM standard mean voxel value'
                  ''
                  'This defines the global mean via a two-step process. Firstly, the overall mean is computed. Voxels with values less than 1/8 of this value are then deemed extra-cranial and get masked out. The mean is then recomputed on the remaining voxels.'
                  ''
}';
% ---------------------------------------------------------------------
% globalc Global calculation
% ---------------------------------------------------------------------
globalc         = cfg_choice;
globalc.tag     = 'globalc';
globalc.name    = 'Global calculation';
globalc.val     = {g_omit };
globalc.help    = {
                   'This option is only used for PET data.'
                   ''
                   'There are three methods for estimating global effects (1) Omit (assumming no other options requiring the global value chosen) (2) User defined (enter your own vector of global values) (3) Mean: SPM standard mean voxel value (within per image fullmean/8 mask) '
                   ''
}';
globalc.values  = {g_omit g_user g_mean };
% ---------------------------------------------------------------------
% gmsca_no No
% ---------------------------------------------------------------------
gmsca_no         = cfg_const;
gmsca_no.tag     = 'gmsca_no';
gmsca_no.name    = 'No';
gmsca_no.val     = {1};
gmsca_no.help    = {'No overall grand mean scaling'};
% ---------------------------------------------------------------------
% gmscv Grand mean scaled value
% ---------------------------------------------------------------------
gmscv         = cfg_entry;
gmscv.tag     = 'gmscv';
gmscv.name    = 'Grand mean scaled value';
gmscv.help    = {
                 'The default value of 50, scales the global flow to a physiologically realistic value of 50ml/dl/min.'
                 ''
}';
gmscv.strtype = 'e';
gmscv.num     = [Inf 1];
gmscv.val     = {50};
% ---------------------------------------------------------------------
% gmsca_yes Yes
% ---------------------------------------------------------------------
gmsca_yes         = cfg_branch;
gmsca_yes.tag     = 'gmsca_yes';
gmsca_yes.name    = 'Yes';
gmsca_yes.val     = {gmscv };
gmsca_yes.help    = {
                     'Scaling of the overall grand mean simply scales all the data by a common factor such that the mean of all the global values is the value specified. For qualitative data, this puts the data into an intuitively accessible scale without altering the statistics. '
                     ''
}';
% ---------------------------------------------------------------------
% gmsca Overall grand mean scaling
% ---------------------------------------------------------------------
gmsca         = cfg_choice;
gmsca.tag     = 'gmsca';
gmsca.name    = 'Overall grand mean scaling';
gmsca.val     = {gmsca_no };
gmsca.help    = {
                 'Scaling of the overall grand mean simply scales all the data by a common factor such that the mean of all the global values is the value specified. For qualitative data, this puts the data into an intuitively accessible scale without altering the statistics. '
                 ''
                 'When proportional scaling global normalisation is used each image is separately scaled such that it''s global value is that specified (in which case the grand mean is also implicitly scaled to that value). So, to proportionally scale each image so that its global value is eg. 20, select <Yes> then type in 20 for the grand mean scaled value.'
                 ''
                 'When using AnCova or no global normalisation, with data from different subjects or sessions, an intermediate situation may be appropriate, and you may be given the option to scale group, session or subject grand means separately. '
                 ''
}';
gmsca.values  = {gmsca_no gmsca_yes };
% ---------------------------------------------------------------------
% glonorm Normalisation
% ---------------------------------------------------------------------
glonorm         = cfg_menu;
glonorm.tag     = 'glonorm';
glonorm.name    = 'Normalisation';
glonorm.help    = {
                   'Global nuisance effects are usually accounted for either by scaling the images so that they all have the same global value (proportional scaling), or by including the global covariate as a nuisance effect in the general linear model (AnCova). Much has been written on which to use, and when. Basically, since proportional scaling also scales the variance term, it is appropriate for situations where the global measurement predominantly reflects gain or sensitivity. Where variance is constant across the range of global values, linear modelling in an AnCova approach has more flexibility, since the model is not restricted to a simple proportional regression. '
                   ''
                   '''Ancova by subject'' or ''Ancova by effect'' options are implemented using the ANCOVA options provided where each experimental factor (eg. subject or effect), is defined. These allow eg. different subjects to have different relationships between local and global measurements. '
                   ''
                   'Since differences between subjects may be due to gain and sensitivity effects, AnCova by subject could be combined with "grand mean scaling by subject" (an option also provided where each experimental factor is originally defined) to obtain a combination of between subject proportional scaling and within subject AnCova. '
                   ''
}';
glonorm.labels = {
                  'None'
                  'Proportional'
                  'ANCOVA'
}';
glonorm.values = {1 2 3};
glonorm.val    = {1};
% ---------------------------------------------------------------------
% globalm Global normalisation
% ---------------------------------------------------------------------
globalm         = cfg_branch;
globalm.tag     = 'globalm';
globalm.name    = 'Global normalisation';
globalm.val     = {gmsca glonorm };
globalm.help    = {
                   'This option is only used for PET data.'
                   ''
                   'Global nuisance effects are usually accounted for either by scaling the images so that they all have the same global value (proportional scaling), or by including the global covariate as a nuisance effect in the general linear model (AnCova). Much has been written on which to use, and when. Basically, since proportional scaling also scales the variance term, it is appropriate for situations where the global measurement predominantly reflects gain or sensitivity. Where variance is constant across the range of global values, linear modelling in an AnCova approach has more flexibility, since the model is not restricted to a simple proportional regression. '
                   ''
                   '''Ancova by subject'' or ''Ancova by effect'' options are implemented using the ANCOVA options provided where each experimental factor (eg. subject or effect), is defined. These allow eg. different subjects to have different relationships between local and global measurements. '
                   ''
                   'Since differences between subjects may be due to gain and sensitivity effects, AnCova by subject could be combined with "grand mean scaling by subject" (an option also provided where each experimental factor is originally defined) to obtain a combination of between subject proportional scaling and within subject AnCova. '
                   ''
}';
% ---------------------------------------------------------------------
% factorial_design Factorial design specification
% ---------------------------------------------------------------------
factorial_design         = cfg_exbranch;
factorial_design.tag     = 'factorial_design';
factorial_design.name    = 'Factorial design specification';
factorial_design.val     = {dir des generic masking globalc globalm };
factorial_design.help    = {
                            'This interface is used for setting up analyses of PET data. It is also used for ''2nd level'' or ''random effects'' analysis which allow one to make a population inference. First level models can be used to produce appropriate summary data, which can then be used as raw data for a second-level analysis. For example, a simple t-test on contrast images from the first-level turns out to be a random-effects analysis with random subject effects, inferring for the population based on a particular sample of subjects.'
                            ''
                            'This interface configures the design matrix, describing the general linear model, data specification, and other parameters necessary for the statistical analysis. These parameters are saved in a configuration file (SPM.mat), which can then be passed on to spm_spm.m which estimates the design. This is achieved by pressing the ''Estimate'' button. Inference on these estimated parameters is then handled by the SPM results section. '
                            ''
                            'A separate interface handles design configuration for fMRI time series.'
                            ''
                            'Various data and parameters need to be supplied to specify the design (1) the image files, (2) indicators of the corresponding condition/subject/group (2) any covariates, nuisance variables, or design matrix partitions (3) the type of global normalisation (if any) (4) grand mean scaling options (5) thresholds and masks defining the image volume to analyse. The interface supports a comprehensive range of options for all these parameters.'
                            ''
}';
factorial_design.prog = @spm_run_factorial_design;
factorial_design.vout = @vout_stats;
%-------------------------------------------------------------------------

%-------------------------------------------------------------------------
function dep = vout_stats(job)
dep(1)            = cfg_dep;
dep(1).sname      = 'SPM.mat File';
dep(1).src_output = substruct('.','spmmat');
dep(1).tgt_spec   = cfg_findspec({{'filter','mat','strtype','e'}});
%dep(2)            = cfg_dep;
%dep(2).sname      = 'SPM Variable';
%dep(2).src_output = substruct('.','spmvar');
%dep(2).tgt_spec   = cfg_findspec({{'strtype','e'}});

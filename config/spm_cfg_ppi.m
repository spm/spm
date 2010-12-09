function ppis = spm_cfg_ppi
% SPM Configuration file for PPIs
%__________________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% Guillaume Flandin
% $Id: spm_cfg_ppi.m 4136 2010-12-09 22:22:28Z guillaume $

% ---------------------------------------------------------------------
% spmmat Select SPM.mat
% ---------------------------------------------------------------------
spmmat         = cfg_files;
spmmat.tag     = 'spmmat';
spmmat.name    = 'Select SPM.mat';
spmmat.help    = {'Select SPM.mat file'};
spmmat.filter = 'mat';
spmmat.ufilter = '^SPM\.mat$';
spmmat.num     = [1 1];

% ---------------------------------------------------------------------
% voi1 Select VOI.mat
% ---------------------------------------------------------------------
voi1         = cfg_files;
voi1.tag     = 'voi';
voi1.name    = 'Select VOI';
voi1.help    = {'physiological variable'};
voi1.filter = 'mat';
voi1.ufilter = '^VOI.*\.mat$';
voi1.num     = [1 1];

% ---------------------------------------------------------------------
% voi2 Select VOI.mat
% ---------------------------------------------------------------------
voi2         = cfg_files;
voi2.tag     = 'voi';
voi2.name    = 'Select VOI';
voi2.help    = {'physiological variables'};
voi2.filter  = 'mat';
voi2.ufilter = '^VOI.*\.mat$';
voi2.num     = [2 2];

% ---------------------------------------------------------------------
% con Matrix of input variables and contrast weights
% ---------------------------------------------------------------------
con         = cfg_entry;
con.tag     = 'u';
con.name    = ' Input variables and contrast weights';
con.help    = {['Matrix of input variables and contrast weights.',... 
'This is an [n x 3] matrix. The first column indexes SPM.Sess.U(i). ',...
'The second column indexes the name of the input or cause, see ',...
'SPM.Sess.U(i).name{j}. The third column is the contrast weight. ',...
' Unless there are parametric effects the second column will generally ',...
'be a 1.']};
con.strtype = 'e';
con.num     = [Inf 3];

% ---------------------------------------------------------------------
% sd Simple deconvolution
% ---------------------------------------------------------------------
sd         = cfg_branch;
sd.tag     = 'sd';
sd.name    = 'Simple deconvolution';
sd.val     = {voi1};
sd.help    = {'Simple deconvolution'};

% ---------------------------------------------------------------------
% phipi Physio-Physiologic Interaction
% ---------------------------------------------------------------------
phipi         = cfg_branch;
phipi.tag     = 'phipi';
phipi.name    = 'Physio-Physiologic Interaction';
phipi.val     = {voi2};
phipi.help    = {'Physio-Physiologic Interaction'};

% ---------------------------------------------------------------------
% ppi Psycho-Physiologic Interaction
% ---------------------------------------------------------------------
ppi         = cfg_branch;
ppi.tag     = 'ppi';
ppi.name    = 'Psycho-Physiologic Interaction';
ppi.val     = {voi1 con};
ppi.help    = {'Psycho-Physiologic Interaction'};

% ---------------------------------------------------------------------
% ppiflag Type of analysis
% ---------------------------------------------------------------------
ppiflag         = cfg_choice;
ppiflag.tag     = 'type';
ppiflag.name    = 'Type of analysis';
ppiflag.help    = {'Type of analysis'};
ppiflag.values  = {sd ppi phipi};

% ---------------------------------------------------------------------
% name Name of PPI
% ---------------------------------------------------------------------
name         = cfg_entry;
name.tag     = 'name';
name.name    = 'Name of PPI';
name.help    = {['Name of the PPI mat file that will be saved in the ',...
'same directory than the SPM.mat file. A ''PPI_'' prefix will be added.']};
name.strtype = 's';
name.num     = [1 Inf];

% ---------------------------------------------------------------------
% showGraphics Display results
% ---------------------------------------------------------------------
showGraphics         = cfg_menu;
showGraphics.tag     = 'disp';
showGraphics.name    = 'Display results';
showGraphics.help    = {'Display results'};
showGraphics.labels  = {'Yes' 'No'};
showGraphics.values  = { 1     0  }; 
showGraphics.val     = { 0 };

% ---------------------------------------------------------------------
% ppis PPI
% ---------------------------------------------------------------------
ppis         = cfg_exbranch;
ppis.tag     = 'ppi';
ppis.name    = 'Physio/Psycho-Physiologic Interaction';
ppis.val     = {spmmat ppiflag name showGraphics};
ppis.help    = {['Bold deconvolution to create physio- or '...
    'psycho-physiologic interactions.']};
ppis.prog    = @run_ppi;
ppis.vout    = @vout_ppi;
%-------------------------------------------------------------------------

%-------------------------------------------------------------------------
function out = run_ppi(job)
switch char(fieldnames(job.type))
    case 'sd'
        PPI = spm_peb_ppi(job.spmmat{1}, 'sd', char(job.type.sd.voi),...
            [], job.name, job.disp);
    case 'ppi'
        PPI = spm_peb_ppi(job.spmmat{1}, 'ppi', char(job.type.ppi.voi),...
            job.type.ppi.u, job.name, job.disp);
    case 'phipi'
        PPI = spm_peb_ppi(job.spmmat{1}, 'phipi', char(job.type.phipi.voi),...
            [], job.name, job.disp);
    otherwise
        error('Unknown type of analysis.');
end
[p n] = fileparts(job.name);
out.ppimat = cellstr(fullfile(fileparts(job.spmmat{1}),['PPI_' n '.mat']));

%-------------------------------------------------------------------------
function dep = vout_ppi(varargin)
dep(1)            = cfg_dep;
dep(1).sname      = ' PPI mat File';
dep(1).src_output = substruct('.','ppimat');
dep(1).tgt_spec   = cfg_findspec({{'filter','mat','strtype','e'}});

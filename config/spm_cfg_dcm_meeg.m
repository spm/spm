function dcm_meeg_spec = spm_cfg_dcm_meeg
% Invert multiple DCMs specified in GUI.
%__________________________________________________________________________
% Copyright (C) 2013-2014 Wellcome Trust Centre for Neuroimaging

% Vladimir Litvak
% $Id: spm_cfg_dcm_meeg.m 6722 2016-02-17 15:04:11Z vladimir $

%--------------------------------------------------------------------------
% D
%--------------------------------------------------------------------------
D        = cfg_files;
D.tag    = 'D';
D.name   = 'M/EEG datasets';
D.filter = 'mat';
D.num    = [1 Inf];
D.help   = {'Select the M/EEG mat files.'};

%--------------------------------------------------------------------------
% DCM
%--------------------------------------------------------------------------
DCM        = cfg_files;
DCM.tag    = 'dcmmat';
DCM.name   = 'DCM files';
DCM.filter = 'mat';
DCM.num    = [1 Inf];
DCM.help   = {'Select pre-specified DCM files.'};

%--------------------------------------------------------------------------
% Prior DCM
%--------------------------------------------------------------------------
pE        = cfg_files;
pE.tag    = 'pE';
pE.name   = 'Priors';
pE.filter = 'mat';
pE.num    = [0 1];
pE.val    = {{''}};
pE.help   = {'Select a DCM file where priors will be taken from (DCM.M.pE, DCM.M.pC)'};

%--------------------------------------------------------------------------
% Initialisation DCM
%--------------------------------------------------------------------------
P        = cfg_files;
P.tag    = 'P';
P.name   = 'Initialisation';
P.filter = 'mat';
P.num    = [0 1];
P.val    = {{''}};
P.help   = {'Select a DCM file for initialising the inversion (at DCM.Ep)'};

%--------------------------------------------------------------------------
% Feedback
%--------------------------------------------------------------------------
feedback        = cfg_menu;
feedback.tag    = 'feedback';
feedback.name   = 'Graphical feedback';
feedback.help   = {'Plot intermediate results during inversions'};
feedback.labels = {'Yes', 'No'};
feedback.values = {1, 0};
feedback.val    = {1};

%--------------------------------------------------------------------------
% Output
%--------------------------------------------------------------------------
output        = cfg_menu;
output.tag    = 'output';
output.name   = 'Output';
output.help   = {'Specify what to output'};
output.labels = {'Single GCM_*.mat file', 'Separate DCM files'};
output.values = {'GCM', 'DCM'};
output.val    = {'GCM'};

%--------------------------------------------------------------------------
% dir Directory
%--------------------------------------------------------------------------
dir         = cfg_files;
dir.tag     = 'dir';
dir.name    = 'Directory';
dir.help    = {'Select the directory where the output will be written.'};
dir.filter  = 'dir';
dir.ufilter = '.*';
dir.num     = [1 1];

%--------------------------------------------------------------------------
% dcm_meeg_spec
%--------------------------------------------------------------------------
dcm_meeg_spec          = cfg_exbranch;
dcm_meeg_spec.tag      = 'meeg';
dcm_meeg_spec.name     = 'DCM for M/EEG';
dcm_meeg_spec.val      = {D, DCM, pE, P, feedback, output, dir};
dcm_meeg_spec.help     = {'Specify DCMs for multiple models and subjects'}';
dcm_meeg_spec.prog     = @meeg_dcm;
dcm_meeg_spec.vout     = @vout_dcm_meeg;
dcm_meeg_spec.modality = {'EEG'};

%==========================================================================
function out = meeg_dcm(job)

if ~isempty(char(job.pE))
    pE = load(char(job.pE));
    pC = pE.DCM.M.pC;
    pE = pE.DCM.M.pE;
else
    pE = [];
    pC = [];
end

if ~isempty(char(job.P))
    P  = load(char(job.P));
    P  = P.DCM.Ep;
else
    P  = [];
end

DCMs     = char(job.dcmmat);
DD       =  char(job.D);
feedback = job.feedback;

GCM = {};
for i = 1:size(DCMs, 1)
    cDCM = getfield(load(deblank(DCMs(i, :)), 'DCM'), 'DCM');
    
    if i==1
        base = spm_file(cDCM.name, 'basename');
    end
    
    % initialise with posteriors if required
    % ---------------------------------------------------------------------
    if isequal(P, 1)
        cDCM.M.P = cDCM.Ep;
    else
        cDCM.M.P = P;
    end
    
    % initialise with posteriors if required
    % ---------------------------------------------------------------------
    if isempty(pE)
        if isfield(cDCM.M,'pE')
            cDCM.M = rmfield(cDCM.M,'pE');
        end
        if isfield(cDCM.M,'pC')
            cDCM.M = rmfield(cDCM.M,'pC');
        end
    elseif ~isequal(pE, 1)
        cDCM.M.pE = pE;
        if ~isempty(pC)
            cDCM.M.pC = pC;
        end
    end
    
    for j = 1:size(DD, 1)
        DCM = cDCM;
        
        D = spm_eeg_load(deblank(DD(j, :)));
        
        [D, ok] = check(D, 'dcm');
        
        if ~ok
            if check(D, 'basic')
                warning (['The file ' D.fname ' is not ready for DCM.'...
                    'Use prepare to specify sensors and fiducials or LFP channels.']);
            else
                warning(['The meeg file ' D.fname ' is corrupt or incomplete']);
            end
            continue;
        end
        
        
        DCM.xY.Dfile  = fullfile(D);
        DCM           = spm_dcm_erp_data(DCM);
        DCM.name      = fullfile(char(job.dir), [spm_file(DCM.name, 'basename') '_' D.fname]);
        DCM.M.nograph = ~feedback;
        
        % save
        %------------------------------------------------------------------
        if isequal(job.output, 'GCM')
            GCM{j, i} = DCM;
        else
            save(DCM.name, 'DCM', spm_get_defaults('mat.format'));
            GCM{j, i} = DCM.name;
        end
    end
end

out = [];
if isequal(job.output, 'GCM')
   out.gcmmat = {fullfile(char(job.dir), ['GCM_' base '.mat'])};
   save(char(out.gcmmat), 'GCM', spm_get_defaults('mat.format')); 
else
   out.dcmmat = GCM(:);
end

%==========================================================================
function dep = vout_dcm_meeg(job)
%==========================================================================

if isequal(job.output, 'GCM')
   dep(1)            = cfg_dep;
    dep(1).sname      = 'GCM mat File(s)';
    dep(1).src_output = substruct('.','gcmmat');
    dep(1).tgt_spec   = cfg_findspec({{'filter','mat','strtype','e'}});
else
    dep(1)            = cfg_dep;
    dep(1).sname      = 'DCM mat File(s)';
    dep(1).src_output = substruct('.','dcmmat');
    dep(1).tgt_spec   = cfg_findspec({{'filter','mat','strtype','e'}});
end
function labview = spm_cfg_opm_read_lvm
% configuration file for reading lab view file
%__________________________________________________________________________
% Copyright (C) 2018 Wellcome Trust Centre for Neuroimaging

% Tim Tierney
% $Id: spm_cfg_opm_read_lvm.m 7429 2018-09-28 09:29:20Z spm $

%--------------------------------------------------------------------------
% labview file
%--------------------------------------------------------------------------
filename        = cfg_files;
filename.tag    = 'filename';
filename.name   = 'File Name';
filename.filter = '(.lvm|.zip)';
filename.num    = [1 1];
filename.help   = {'Select the (zipped) lvm file.'};


%--------------------------------------------------------------------------
% headerlength
%--------------------------------------------------------------------------
headerlength         = cfg_entry;
headerlength.tag     = 'headerlength';
headerlength.name    = 'No. of header lines';
headerlength.help    = {'The number of lines of text containing header information'};
headerlength.strtype = 'r';
headerlength.num     = [1,1];
headerlength.val     = {23};



%--------------------------------------------------------------------------
% timeind
%--------------------------------------------------------------------------
timeind         = cfg_entry;
timeind.tag     = 'timeind';
timeind.name    = 'Time Index';
timeind.help    = {'The column number of the time variable'};
timeind.strtype = 'r';
timeind.num     = [1,1];
timeind.val     = {1};

%--------------------------------------------------------------------------
% Decimal Triggers
%--------------------------------------------------------------------------
decimalTriggerInds         = cfg_entry;
decimalTriggerInds.tag     = 'decimalTriggerInds';
decimalTriggerInds.name    = 'Index of Decimal Triggers';
decimalTriggerInds.help    = {'Column numbers of channels that should be interpreted as decimal triggers. '};
decimalTriggerInds.strtype = 'r';
decimalTriggerInds.num     = [0,0];
decimalTriggerInds.val     = {74:81};

%--------------------------------------------------------------------------
% Binary Triggers
%--------------------------------------------------------------------------
binaryTriggerInds         = cfg_entry;
binaryTriggerInds.tag     = 'binaryTriggerInds';
binaryTriggerInds.name    = 'Index of Binary Triggers';
binaryTriggerInds.help    = {'Column numbers of channels that should be interpreted as binary triggers'};
binaryTriggerInds.strtype = 'r';
binaryTriggerInds.num     = [0,0];
binaryTriggerInds.val     = {[]};

%--------------------------------------------------------------------------
% Trigger Threshld
%--------------------------------------------------------------------------
trigThresh         = cfg_entry;
trigThresh.tag     = 'trigThresh';
trigThresh.name    = 'Trigger Threshold';
trigThresh.help    = {'Threshold to apply to triggers. Currently the default has units of Volts.'};
trigThresh.strtype = 'r';
trigThresh.num     = [1,1];
trigThresh.val     = {4};

%--------------------------------------------------------------------------
% read
%--------------------------------------------------------------------------
labview          = cfg_exbranch;
labview.tag      = 'labview';
labview.name     = 'Read LabView Files';
labview.val      = {filename,headerlength,timeind,decimalTriggerInds,binaryTriggerInds,trigThresh};
labview.help     = {'Reading LabView data'}';
labview.prog     = @lbv_read;
labview.vout     = @vout_lbv_read;
labview.modality = {'EEG'};


%==========================================================================
function labview = lbv_read(job)
data = spm_opm_read_lvm(job);

outfile = spm_file(job.filename{1},'ext','.mat');
save(outfile,'data');

labview.data = {outfile};


%==========================================================================
function dep = vout_lbv_read(job)
% return dependencies
dep = cfg_dep;
dep.sname = 'Prepared Labview Data';
dep.src_output = substruct('.','data');
dep.tgt_spec   = cfg_findspec({{'filter','mat'}});

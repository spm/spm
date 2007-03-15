function conf = spm_config_fmri_data
% Configuration file for specification of fMRI model
%_______________________________________________________________________
% Copyright (C) 2005 Wellcome Department of Imaging Neuroscience

% Darren Gitelman and Will Penny
% $Id: spm_config_fmri_data.m 766 2007-03-15 14:09:30Z volkmar $


% Define inline types.
%-----------------------------------------------------------------------

entry = inline(['struct(''type'',''entry'',''name'',name,'...
    '''tag'',tag,''strtype'',strtype,''num'',num,''help'',hlp)'],...
    'name','tag','strtype','num','hlp');

files = inline(['struct(''type'',''files'',''name'',name,'...
    '''tag'',tag,''filter'',fltr,''num'',num,''help'',hlp)'],...
    'name','tag','fltr','num','hlp');

mnu = inline(['struct(''type'',''menu'',''name'',name,'...
    '''tag'',tag,''labels'',{labels},''values'',{values},''help'',hlp)'],...
    'name','tag','labels','values','hlp');

branch = inline(['struct(''type'',''branch'',''name'',name,'...
    '''tag'',tag,''val'',{val},''help'',hlp)'],...
    'name','tag','val','hlp');

repeat = inline(['struct(''type'',''repeat'',''name'',name,'...
    '''tag'',tag,''values'',{values},''help'',hlp)'],...
    'name','tag','values','hlp');

choice = inline(['struct(''type'',''choice'',''name'',name,'...
    '''tag'',tag,''values'',{values},''help'',hlp)'],...
    'name','tag','values','hlp');

%-----------------------------------------------------------------------

sp_text = ['                                                      ',...
    '                                                      '];
%-----------------------------------------------------------------------

scans    = files('Scans','scans','image',[1 Inf],'Select scans');
scans.help = {[...
    'Select the fMRI scans for this session.  They must all have the same ',...
    'image dimensions, orientation, voxel size etc.']};

%-------------------------------------------------------------------------

mask = files('Explicit mask','mask','image',[0 1],'Image mask');
mask.val = {''};
p1=['Specify an image for explicitly masking the analysis. ',...
    'A sensible option here is to use a segmention of structural images ',...
    'to specify a within-brain mask. If you select that image as an ',...
    'explicit mask then only those voxels in the brain will be analysed. ',...
    'This both speeds the estimation and restricts SPMs/PPMs to within-brain ',...
    'voxels. Alternatively, if such structural images are unavailble or no ',...
    'masking is required, then leave this field empty.'];
mask.help={p1};

%-------------------------------------------------------------------------

spmmat = files('Select SPM.mat','spmmat','mat',1,'');
spmmat.help = {[...
    'Select the SPM.mat file containing the ',...
    'specified design matrix.']};

%-------------------------------------------------------------------------

conf = branch('fMRI data specification','fmri_data',...
    {scans,spmmat,mask},'fMRI data');
conf.prog   = @run_stats;
conf.vfiles = @vfiles_stats;
conf.modality = {'FMRI'};
conf.help = {['Select the data and optional explicit mask for a specified ' ...
             'design']};

return;
%-------------------------------------------------------------------------

%-------------------------------------------------------------------------
function my_cd(varargin)
% jobDir must be the actual directory to change to, NOT the job structure.
jobDir = varargin{1};
if ~isempty(jobDir)
    try
    cd(char(jobDir));
    fprintf('Changing directory to: %s\n',char(jobDir));
    catch
        error('Failed to change directory. Aborting run.')
    end
end
return;
%------------------------------------------------------------------------

%------------------------------------------------------------------------
function run_stats(job)
% Set up the design matrix and run a design.

spm_defaults;
global defaults
defaults.modality='FMRI';

original_dir = pwd;
[p n e v] = fileparts(job.spmmat{1});
my_cd(p);
load(job.spmmat{1});    
% Image filenames
%-------------------------------------------------------------
SPM.xY.P = strvcat(job.scans);

% Let SPM configure the design
%-------------------------------------------------------------
SPM = spm_fmri_spm_ui(SPM);

if ~isempty(job.mask)&&~isempty(job.mask{1})
    SPM.xM.VM         = spm_vol(job.mask{:});
    SPM.xM.xs.Masking = [SPM.xM.xs.Masking, '+explicit mask'];
end

%-Save SPM.mat
%-----------------------------------------------------------------------
fprintf('%-40s: ','Saving SPM configuration')   %-#
if str2num(version('-release'))>=14,
    save('SPM','-V6','SPM');
else
    save('SPM','SPM');
end;

fprintf('%30s\n','...SPM.mat saved')                     %-#


my_cd(original_dir); % Change back dir
fprintf('Done\n')
return
%-------------------------------------------------------------------------

%-------------------------------------------------------------------------
function vf = vfiles_stats(job)
vf    = job.spmmat;


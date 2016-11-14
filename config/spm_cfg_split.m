function split = spm_cfg_split
% SPM Configuration file for 4D to 3D volumes conversion
%__________________________________________________________________________
% Copyright (C) 2014 Wellcome Trust Centre for Neuroimaging

% Guillaume Flandin
% $Id: spm_cfg_split.m 6929 2016-11-14 13:07:31Z guillaume $


%--------------------------------------------------------------------------
% vols 4D Volume
%--------------------------------------------------------------------------
vol         = cfg_files;
vol.tag     = 'vol';
vol.name    = '4D Volume';
vol.help    = {'Select the 4D volume file to convert into a series of 3D volume files.'};
vol.filter  = 'image';
vol.ufilter = '.*';
vol.num     = [1 1];

%--------------------------------------------------------------------------
% outdir Output Directory
%--------------------------------------------------------------------------
outdir         = cfg_files;
outdir.tag     = 'outdir';
outdir.name    = 'Output Directory';
outdir.val     = {{''}};
outdir.help    = {'Specify the output directory. If no directory is given, files will be written in the same directory than the input 4D file.'};
outdir.filter  = 'dir';
outdir.ufilter = '.*';
outdir.num     = [0 1];

%--------------------------------------------------------------------------
% cat 4D to 3D File Conversion
%--------------------------------------------------------------------------
split      = cfg_exbranch;
split.tag  = 'split';
split.name = '4D to 3D File Conversion';
split.val  = {vol outdir};
split.help = {'Convert a 4D volume file into a series of 3D volume files.'};
split.prog = @(job)spm_run_split('run',job);
split.vout = @(job)spm_run_split('vout',job);


%==========================================================================
function out = spm_run_split(cmd, job)

switch lower(cmd)
    case 'run'
        V                 = char(job.vol);
        odir              = job.outdir;
        if isempty(odir{1}), odir = {}; end
        Vo                = spm_file_split(V, odir{:});
        out.splitfiles    = {Vo.fname}';
        
    case 'vout'
        out(1)            = cfg_dep;
        out(1).sname      = 'Series of 3D Volumes';
        out(1).src_output = substruct('.','splitfiles');
        out(1).tgt_spec   = cfg_findspec({{'filter','image','strtype','e'}});
end

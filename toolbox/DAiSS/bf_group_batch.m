function res = bf_group_batch(BF, S)
% Run a DAiSS batch on a group of subjects
% Copyright (C) 2015 Wellcome Trust Centre for Neuroimaging

% Vladimir Litvak
% $Id: bf_group_batch.m 7703 2019-11-22 12:06:29Z guillaume $

%--------------------------------------------------------------------------
if nargin == 0 
    batchfile = cfg_files;
    batchfile.tag = 'batchfile';
    batchfile.name = 'Batch .mat or .m file';
    batchfile.filter = '(.*.mat$)|(.*.m$)';
    batchfile.num = [1 Inf];
    batchfile.help = {'Select batch specification file.'};
    
    batch = cfg_branch;
    batch.tag = 'batch';
    batch.name = 'Batch';
    batch.val = {batchfile};
    
    res = batch;
    
    return
elseif nargin < 2
    error('Two input arguments are required');
end

batchfile = char(S.batchfile);

if isequal(spm_file(batchfile, 'ext'), 'm')
    cdir = pwd;
    if ~isempty(spm_file(batchfile, 'path'))
        try
            cd(spm_file(batchfile, 'path'));
        catch
            % This is to allow the use of pipelines saved with DAiSS
            % on different machines
            tbxdir = fileparts(mfilename('fullpath'));
            cd(tbxdir);
        end
    end
    vars = who;
    eval(spm_file(batchfile, 'basename'));
    cd(cdir);
    name = setdiff(who, [vars; {'vars'}]);
    if numel(name)~= 1
        error('Invalid batch specification');
    end
    matlabbatch = eval(char(name));
    if ~isa(matlabbatch, 'cell')
        error('Invalid batch specification');
    end
elseif isequal(spm_file(batchfile, 'ext'), 'mat')
    try
        tmp  = load(batchfile);
    catch
        tbxdir = fileparts(mfilename('fullpath'));
        tmp  = load(spm_file(batchfile, 'path', tbxdir));
    end
    name = fieldnames(tmp);
    if numel(name)~= 1
        error('Invalid batch specification');
    end
    matlabbatch = tmp.(char(name));
    if ~isa(matlabbatch, 'cell')
        error('Invalid batch specification');
    end
else
    error('Invalid batch specification');
end

try 
    matlabbatch{1}.spm.tools.beamforming;
    [~, matlabbatch] = spm_jobman('harvest', matlabbatch);
catch
    error('Invalid batch specification. DAiSS batch expected.');
end

res = cell(1, numel(BF));
for i = 1:numel(BF)
    if isfield(matlabbatch{1}.spm.tools.beamforming, 'data');
        D = spm_eeg_load(BF{i});
        matlabbatch{1}.spm.tools.beamforming.data.D = {fullfile(D)};
        dum = mkdir(D.path, [S.prefix 'BF']);
        matlabbatch{1}.spm.tools.beamforming.data.dir = {fullfile(D.path, [S.prefix 'BF'])};        
    else
        module = char(fieldnames(matlabbatch{1}.spm.tools.beamforming));
        matlabbatch{1}.spm.tools.beamforming.(module).BF = BF(i);        
    end
    out = spm_jobman('run', matlabbatch);
    res(i) = out{end}.BF;
end

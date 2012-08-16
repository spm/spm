function out = bf_sources
% Prepares source locations and lead fields for beamforming
% Copyright (C) 2012 Wellcome Trust Centre for Neuroimaging

% Vladimir Litvak
% $Id: bf_sources.m 4847 2012-08-16 17:29:23Z vladimir $

% dir Directory
% ---------------------------------------------------------------------
BF = cfg_files;
BF.tag = 'BF';
BF.name = 'BF.mat file';
BF.filter = '^BF.mat$';
BF.num = [1 1];
BF.help = {'Select BF.mat file.'};

%--------------------------------------------------------------------------
% method
%--------------------------------------------------------------------------
plugin      = cfg_choice;
plugin.tag  = 'plugin';
plugin.name = 'Source space type ';

source_funs = spm_select('List', fullfile(spm('dir'), 'toolbox', 'Beamforming'), '^bf_sources_.*\.m$');
source_funs = cellstr(source_funs );
for i = 1:numel(source_funs)
    plugin.values{i} = feval(spm_file(source_funs{i},'basename'));
end


out = cfg_exbranch;
out.tag = 'bf_sources';
out.name = 'Define sources';
out.val = {BF, plugin};
out.help = {'Define source space for beamforming'};
out.prog = @bf_source_run;
out.vout = @bf_source_vout;
out.modality = {'EEG'};
end

function  out = bf_source_run(job)

outdir = spm_file(job.BF{1}, 'fpath');

cd(outdir);

BF = bf_load('BF.mat');

plugin_name = cell2mat(fieldnames(job.plugin));

BF.sources.(plugin_name) = feval(['bf_sources_' plugin_name], BF, job.plugin.(plugin_name));
BF.sources.pos = BF.sources.(plugin_name).pos;
if isfield(BF.sources.(plugin_name), 'ori')
    BF.sources.ori = BF.sources.(plugin_name).ori;
else
    BF.sources.ori = [];
end

nvert = size(BF.sources.pos, 1);
modalities = {'MEG', 'EEG'};

for m = 1:numel(modalities)
    
    if isfield(BF.data, modalities{m})
        
        
        [vol, sens] = ft_prepare_vol_sens(BF.data.(modalities{m}).vol, BF.data.(modalities{m}).sens, 'channel', ...
            chanlabels(BF.data.D, meegchannels(BF.data.D, modalities{m})));
        
        
        spm('Pointer', 'Watch');drawnow;
        spm_progress_bar('Init', nvert, ['Computing ' modalities{m} ' leadfields']); drawnow;
        if nvert > 100, Ibar = floor(linspace(1, nvert,100));
        else Ibar = 1:nvert; end
        
        L = cell(1, nvert);
        
        for i = 1:nvert
            if ft_inside_vol(BF.sources.pos(i, :), vol)
                
                L{i}  = ft_compute_leadfield(BF.sources.pos(i, :), sens, vol);
                
                if ~isempty(BF.sources.ori)
                    L{i}  = L{i}*BF.sources.ori(i, :)';
                end
            else
                L{i} = NaN;
            end
            
            if ismember(i, Ibar)
                spm_progress_bar('Set', i); drawnow;
            end
        end
         
        spm_progress_bar('Clear');
        
        
        BF.sources.L.(modalities{m}) = L;
    end
end

bf_save(BF);

out.BF{1} = fullfile(outdir, 'BF.mat');
end

function dep = bf_source_vout(job)
% Output is always in field "D", no matter how job is structured
dep = cfg_dep;
dep.sname = 'BF.mat file';
% reference field "B" from output
dep.src_output = substruct('.','BF');
% this can be entered into any evaluated input
dep.tgt_spec   = cfg_findspec({{'filter','mat'}});
end

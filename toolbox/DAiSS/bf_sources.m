function out = bf_sources
% Prepares source locations and lead fields for beamforming
% Copyright (C) 2012 Wellcome Trust Centre for Neuroimaging

% Vladimir Litvak
% $Id: bf_sources.m 7703 2019-11-22 12:06:29Z guillaume $

% dir Directory
% ---------------------------------------------------------------------
BF = cfg_files;
BF.tag = 'BF';
BF.name = 'BF.mat file';
BF.filter = '^BF.mat$';
BF.num = [1 1];
BF.help = {'Select BF.mat file.'};

reduce_rank = cfg_entry;
reduce_rank.tag = 'reduce_rank';
reduce_rank.name = 'Reduce rank';
reduce_rank.strtype = 'r';
reduce_rank.num = [1 2];
reduce_rank.val = {[2 3]};
reduce_rank.help = {'Enter rank for MEG and EEG lead fields [MEG EEG]'};

keep3d = cfg_menu;
keep3d.tag = 'keep3d';
keep3d.name = 'Keep the original orientations';
keep3d.help = {'If not then the extra dimension is physically removed.'};
keep3d.labels = {'yes', 'no'};
keep3d.values = {1, 0};
keep3d.val = {1};

%--------------------------------------------------------------------------
% method
%--------------------------------------------------------------------------
plugin      = cfg_choice;
plugin.tag  = 'plugin';
plugin.name = 'Source space type ';

source_funs = spm_select('List', fileparts(mfilename('fullpath')), '^bf_sources_.*\.m$');
source_funs = cellstr(source_funs );
for i = 1:numel(source_funs)
    plugin.values{i} = feval(spm_file(source_funs{i},'basename'));
end

visualise = cfg_menu;
visualise.tag = 'visualise';
visualise.name = 'Visualise head model and sources';
visualise.help = {'Visualise head model and sourses to verify that everythin was done correctly'};
visualise.labels = {'yes', 'no'};
visualise.values = {1, 0};
visualise.val = {1};

out = cfg_exbranch;
out.tag = 'sources';
out.name = 'Define sources';
out.val = {BF, reduce_rank, keep3d, plugin, visualise};
out.help = {'Define source space for beamforming'};
out.prog = @bf_source_run;
out.vout = @bf_source_vout;
out.modality = {'EEG'};
end

function  out = bf_source_run(job)

outdir = spm_file(job.BF{1}, 'fpath');


BF = bf_load(fullfile(outdir, 'BF.mat'));

plugin_name = cell2mat(fieldnames(job.plugin));

field_name  = strtok(plugin_name, '_');

BF.sources = [];
BF.sources.(field_name) = feval(['bf_sources_' plugin_name], BF, job.plugin.(plugin_name));
BF.sources.pos = BF.sources.(field_name).pos;

if isfield(BF.sources.(field_name), 'ori')
    BF.sources.ori = BF.sources.(field_name).ori;
else
    BF.sources.ori = [];
end

siunits = isfield(BF.data, 'siunits') & BF.data.siunits;

nvert = size(BF.sources.pos, 1);
modalities = {'MEG', 'EEG'};
reduce_rank=job.reduce_rank; 

for m = 1:numel(modalities)
    
    if isfield(BF.data, modalities{m})
        
        if isequal(modalities{m}, 'MEG')
            chanind = indchantype(BF.data.D, {'MEG', 'MEGPLANAR'}, 'GOOD');
        elseif isequal(modalities{m}, 'EEG')
            chanind = indchantype(BF.data.D, 'EEG', 'GOOD');
        end
        
        if isempty(chanind)
            error(['No good ' modalities{m} ' channels were found.']);
        end
        
        if ischar(BF.data.(modalities{m}).vol),           
            BF.data.(modalities{m}).vol = ft_read_vol(BF.data.(modalities{m}).vol);
        end;
               
        chanunits = units(BF.data.D, chanind);
        
        
        [vol, sens] = ft_prepare_vol_sens(BF.data.(modalities{m}).vol, BF.data.(modalities{m}).sens, 'channel', ...
            chanlabels(BF.data.D, chanind));
        
        pos = BF.sources.pos;
                
        if isfield(BF.data.(modalities{m}), 'mesh_correction') && ~isempty(BF.data.(modalities{m}).mesh_correction)
            disp(['Adjusting source points for volume type ' ft_voltype(vol)]);
            cfg     = BF.data.(modalities{m}).mesh_correction;
            cfg.vol      = vol;
            cfg.grid.pos = pos;
            gridcorrect  = ft_prepare_sourcemodel(cfg);
            
            pos          = gridcorrect.pos;
        end
        
        if job.visualise
            F = spm_figure('GetWin', modalities{m});clf;
            
            if ismac
                set(F,'renderer','zbuffer');
            else
                set(F,'renderer','OpenGL');
            end
            
            ft_plot_headmodel(vol, 'edgecolor', [0 0 0], 'facealpha', 0);
            
            hold on
            
            try
                ft_plot_sens(sens, 'style', '*b', 'coil',  ft_senstype(sens, 'eeg'));
            catch
                ft_plot_sens(sens, 'style', '*b', 'coilshape', 'point', 'coil', ft_senstype(sens, 'eeg'));
            end
            
            plot3(pos(:, 1), pos(:, 2), pos(:, 3), '.r', 'MarkerSize', 10);
            
            rotate3d on;
            
            axis off
            axis vis3d
            axis equal
        end
        
        spm('Pointer', 'Watch');drawnow;
        spm_progress_bar('Init', nvert, ['Computing ' modalities{m} ' leadfields']); drawnow;
        if nvert > 100, Ibar = floor(linspace(1, nvert,100));
        else Ibar = 1:nvert; end
        
        L = cell(1, nvert);       
        
        for i = 1:nvert
            if siunits
                L{i}  = ft_compute_leadfield(BF.sources.pos(i, :), sens, vol, 'reducerank', reduce_rank(m), 'dipoleunit', 'nA*m', 'chanunit', chanunits);
            else
                L{i}  = ft_compute_leadfield(BF.sources.pos(i, :), sens, vol, 'reducerank', reduce_rank(m));
            end
            
            if ~isempty(BF.sources.ori) && any(BF.sources.ori(i, :))               
                L{i}  = L{i}*BF.sources.ori(i, :)';
            elseif ~job.keep3d &&  reduce_rank(m) < 3
                 [U_, S_, V] = svd(L{i}, 'econ');
                 L{i} = L{i}*V(:,1:reduce_rank(m));
            end
            
            if ismember(i, Ibar)
                spm_progress_bar('Set', i); drawnow;
            end
        end
         
        spm_progress_bar('Clear');
        
        BF.sources.reduce_rank.(modalities{m})=reduce_rank(m); %MWW
        BF.sources.L.(modalities{m}) = L;
        BF.sources.channels.(modalities{m}) = chanlabels(BF.data.D, chanind);
    end       
end

bf_save_path(BF,fullfile(outdir, 'BF.mat'));

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

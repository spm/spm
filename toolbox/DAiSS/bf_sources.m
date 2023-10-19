function out = bf_sources
% Prepare source locations and lead fields for beamforming
%__________________________________________________________________________

% Vladimir Litvak
% Copyright (C) 2015-2023 Wellcome Centre for Human Neuroimaging


out = cfg_exbranch;
out.tag = 'sources';
out.name = 'Define sources';
out.val = @bf_sources_cfg;
out.help = {'Define source space for beamforming'};
out.prog = @bf_source_run;
out.vout = @bf_source_vout;
out.modality = {'EEG'};


%==========================================================================
function varargout = bf_sources_cfg

persistent cfg
if ~isempty(cfg), varargout = {cfg}; return; end

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

normalise_lf = cfg_menu;
normalise_lf.tag = 'normalise_lf';
normalise_lf.name = 'Normalise lead fields';
normalise_lf.help = {'Normalise the lead fields to yield array gain beamformer'};
normalise_lf.labels = {'yes', 'no'};
normalise_lf.values = {true, false};
normalise_lf.val = {false};

visualise = cfg_menu;
visualise.tag = 'visualise';
visualise.name = 'Visualise head model and sources';
visualise.help = {'Visualise head model and sourses to verify that everythin was done correctly'};
visualise.labels = {'yes', 'no'};
visualise.values = {true, false};
visualise.val = {true};

[cfg,varargout{1}] = deal({BF, reduce_rank, keep3d, plugin, normalise_lf, visualise});


%==========================================================================
function  out = bf_source_run(job)

outdir = spm_file(job.BF{1}, 'fpath');


BF = bf_load(fullfile(outdir, 'BF.mat'));

plugin_name = cell2mat(fieldnames(job.plugin));

field_name  = strtok(plugin_name, '_');

% odd case for scalp, rename to mesh
if strcmpi(field_name,'scalp')
    field_name = 'mesh';
end

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
            
            % OPM data can have a case where a sensor is selected above, 
            % but does not have a physical location in the sens structure, 
            % find those channels and remove from the list.
            [sel1, ~] = match_str(BF.data.D.chanlabels(chanind),...
                BF.data.(modalities{m}).sens.label);
            missing = setdiff(chanind,chanind(sel1));
            if ~isempty(missing)
                warning(['Sensors below missing location information:'...
                    ' ignoring and setting to a different type'])
                missing_names = BF.data.D.chanlabels(missing);
                for ii = 1:numel(missing)
                    fprintf('\t\t%s',missing_names{ii});
                end
                fprintf('\n');
                BF.data.D = BF.data.D.chantype(missing,'OPM');
                chanind = chanind(sel1);
            end
            
        elseif isequal(modalities{m}, 'EEG')
            chanind = indchantype(BF.data.D, 'EEG', 'GOOD');
        end
        
        if isempty(chanind)
            error(['No good ' modalities{m} ' channels were found.']);
        end
        
        if ischar(BF.data.(modalities{m}).vol)
            BF.data.(modalities{m}).vol = ft_read_vol(BF.data.(modalities{m}).vol);
        end
        
        chanunits = units(BF.data.D, chanind);
        
        [vol, sens] = ft_prepare_vol_sens(BF.data.(modalities{m}).vol, BF.data.(modalities{m}).sens, 'channel', ...
            chanlabels(BF.data.D, chanind));
        
        pos = BF.sources.pos;
        
        if isfield(BF.data.(modalities{m}), 'mesh_correction') && ~isempty(BF.data.(modalities{m}).mesh_correction)
            disp(['Adjusting source points for volume type ' ft_headmodeltype(vol)]);
            cfg     = BF.data.(modalities{m}).mesh_correction;
            cfg.vol      = vol;
            cfg.grid.pos = pos;
            gridcorrect  = ft_prepare_sourcemodel(cfg);
            
            pos          = gridcorrect.pos;
        end
        
        if job.visualise
            F = spm_figure('GetWin', ['Head model for ' modalities{m}]);clf;
            
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
            
            if job.normalise_lf
                L{i} = L{i}./norm(L{i}, 'fro');
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


%==========================================================================
function dep = bf_source_vout(job)
% Output is always in field "D", no matter how job is structured
dep = cfg_dep;
dep.sname = 'BF.mat file';
% reference field "B" from output
dep.src_output = substruct('.','BF');
% this can be entered into any evaluated input
dep.tgt_spec   = cfg_findspec({{'filter','mat'}});

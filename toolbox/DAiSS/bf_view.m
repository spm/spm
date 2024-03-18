function out = bf_view
% Display the results of beamforming analysis
%__________________________________________________________________________

% George O'Neill
% Copyright (C) 2022-2023 Wellcome Centre for Human Neuroimaging


% dir Directory
% ---------------------------------------------------------------------
BF = cfg_files;
BF.tag = 'BF';
BF.name = 'BF.mat file';
BF.filter = '^BF.mat$';
BF.num = [1 1];
BF.help = {'Select BF.mat file.'};

% method
%---------------------------------------------------------------------
plugin      = cfg_choice;
plugin.tag  = 'plugin';
plugin.name = 'What to view';

view_funs = spm_select('List', fileparts(mfilename('fullpath')), '^bf_view_.*\.m$');
view_funs = cellstr(view_funs );
for i = 1:numel(view_funs)
    plugin.values{i} = feval(spm_file(view_funs{i},'basename'));
end

% save?
%---------------------------------------------------------------------
write           = cfg_menu;
write.tag       = 'write';
write.name      = 'Save as figure';
write.help      = {'Save displayed figure as .fig file?'};
write.labels    = {'yes','no'};
write.values    = {true, false};
write.val       = {false};

out = cfg_exbranch;
out.tag = 'view';
out.name = 'View';
out.val = {BF, plugin, write};
out.help = {'View results'};
out.prog = @bf_view_run;
% out.vout = @bf_view_vout;
out.modality = {'EEG'};
end

function  out = bf_view_run(job)

outdir = spm_file(job.BF{1}, 'fpath');

cd(outdir);

BF = bf_load('BF.mat', {'data', 'sources', 'features', 'output'});

plugin_name   = cell2mat(fieldnames(job.plugin));

outfield_name =  strtok(plugin_name, '_');

BF.view.(outfield_name) = feval(['bf_view_' plugin_name], BF, job.plugin.(plugin_name));

if job.write
    if ishghandle(BF.view.glass)
        fname = [BF.output.image.label '.fig'];
        saveas(BF.view.(outfield_name),fname,'fig');
    else
        fprintf('WARNING: Figure handle points to deleted figure, skipping writing to disk\n')
    end
end

end

% function dep = bf_view_vout(job)
%
%
% % % Output is always in field "D", no matter how job is structured
% % dep(1) = cfg_dep;
% % dep(1).sname = 'BF.mat file';
% % % reference field "B" from output
% % dep(1).src_output = substruct('.','BF');
% % % this can be entered into any evaluated input
% % dep(1).tgt_spec   = cfg_findspec({{'filter','mat'}});
% % dep(2).sname      = 'Output files';
% % dep(2).src_output = substruct('.','files');
% %
% % plugin_name   = cell2mat(fieldnames(job.plugin));
% %
% % switch plugin_name
% %     case 'spmeeg'
% %         dep(2).tgt_spec   = cfg_findspec({{'filter','mat','strtype','e'}});
% %     case 'nifti'
% %         dep(2).tgt_spec   = cfg_findspec({{'filter','image','strtype','e'}});
% %     case 'gifti'
% %         dep(2).tgt_spec   = cfg_findspec({{'filter','mesh','strtype','e'}});
% % end
% end
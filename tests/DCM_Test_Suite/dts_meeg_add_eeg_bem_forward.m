function outputs = dts_meeg_add_eeg_bem_forward(eeg_file)
% Build a template EEG BEM forward model and full gain matrix for one M/EEG file.
% Authored by Pranay Yadav in 2026

if exist(eeg_file, 'file') ~= 2
    error('Missing processed EEG fixture: %s', eeg_file);
end

spm('defaults', 'eeg');
spm_get_defaults('cmdline', true);
spm('CmdLine', true);

D = spm_eeg_load(eeg_file);
val = 1;
D.val = val;
D.inv = {struct('mesh', [], 'date', strvcat(date, datestr(now, 15)), 'comment', {{''}})};
D = spm_eeg_inv_mesh_ui(D, val, true, 2);

meegfid = D.fiducials;
mrifid = D.inv{val}.mesh.fid;
D = spm_eeg_inv_datareg_ui(D, val, meegfid, mrifid, 1);

D.inv{val}.forward = struct([]);
D.inv{val}.forward(1).voltype = 'EEG BEM';
D = spm_eeg_inv_forward(D);
save(D);

D = spm_eeg_load(eeg_file);
D.val = val;
Is = 1:size(D.inv{D.val}.forward(1).mesh.vert, 1);
[~, D] = spm_eeg_lgainmat(D, Is);
save(D);

outputs = struct();
outputs.final_mat = D.fullfile();
outputs.final_dat = fullfile(D.path, spm_file(D.fname, 'ext', '.dat'));
outputs.gainmat_file = fullfile(D.path, D.inv{D.val}.gainmat);
end

function fixtures = dts_meeg_sample_sensor_fixture(cfg, purpose)
% Build local IMG/ECD sensor fixtures from the shipped processed EEG data.
% Authored by Pranay Yadav in 2026

if nargin < 2 || isempty(purpose)
    purpose = 'runtime';
end

source_dir = fullfile(cfg.prep_root, 'test_05_sample_sensor_data');
source_file = fullfile(source_dir, 'EEG_MMN.mat');
if exist(source_file, 'file') ~= 2
    error(['Missing processed sensor fixture:\n%s\n' ...
        'Run tests/DCM_Test_Suite/dts_prep01_make_test_fixtures.m or ship the fixture.'], source_file);
end

local_root = fullfile(cfg.prep_root, 'local');
img_dir = fullfile(local_root, 'test_05_sample_sensor_img');
ecd_dir = fullfile(local_root, 'test_05_sample_sensor_ecd');
ready_file = fullfile(local_root, 'sample_sensor_fixture.ready');
lock_dir = fullfile(local_root, 'sample_sensor_fixture.lock');

if fixture_ready(ready_file, img_dir, ecd_dir)
    fixtures = fixture_struct(img_dir, ecd_dir, purpose);
    return
end

if ~exist(local_root, 'dir'), mkdir(local_root); end
running_on_worker = false;
try
    running_on_worker = ~isempty(getCurrentTask());
catch
end
if ~cfg.parallel && ~running_on_worker
    build_sample_sensor_fixture(source_file, img_dir, ecd_dir, ready_file);
    fixtures = fixture_struct(img_dir, ecd_dir, purpose);
    return
end

have_lock = false;
while ~have_lock
    if exist(lock_dir, 'dir')
        have_lock = false;
    else
        [created_lock, ~, lock_msgid] = mkdir(lock_dir);
        have_lock = created_lock && isempty(lock_msgid);
    end
    if have_lock
        cleanup_lock = onCleanup(@()remove_lock(lock_dir));
    elseif fixture_ready(ready_file, img_dir, ecd_dir)
        fixtures = fixture_struct(img_dir, ecd_dir, purpose);
        return
    else
        pause(2);
    end
end

if fixture_ready(ready_file, img_dir, ecd_dir)
    fixtures = fixture_struct(img_dir, ecd_dir, purpose);
    return
end

build_sample_sensor_fixture(source_file, img_dir, ecd_dir, ready_file);
fixtures = fixture_struct(img_dir, ecd_dir, purpose);
end

function build_sample_sensor_fixture(source_file, img_dir, ecd_dir, ready_file)
if exist(img_dir, 'dir'), rmdir(img_dir, 's'); end
if exist(ecd_dir, 'dir'), rmdir(ecd_dir, 's'); end
mkdir(img_dir);
mkdir(ecd_dir);

D = spm_eeg_load(source_file);
Dimg = copy(D, fullfile(img_dir, 'EEG_MMN.mat'));
Dimg = strip_forward_model(Dimg);
save(Dimg);
img_outputs = dts_meeg_add_eeg_bem_forward(Dimg.fullfile());

D = spm_eeg_load(img_outputs.final_mat);
Decd = copy(D, fullfile(ecd_dir, 'EEG_MMN.mat'));
copy_gainmatrix(D, ecd_dir);
save(Decd);

fid = fopen(ready_file, 'w');
if fid == -1
    error('Could not write sample sensor fixture ready marker: %s', ready_file);
end
fprintf(fid, '%s\n', char(datetime('now', 'Format', 'yyyyMMdd''T''HHmmss')));
fclose(fid);
end

function D = strip_forward_model(D)
try
    D.inv = {};
catch
end
end

function tf = fixture_ready(ready_file, img_dir, ecd_dir)
tf = exist(ready_file, 'file') == 2 && ...
    exist(fullfile(img_dir, 'EEG_MMN.mat'), 'file') == 2 && ...
    exist(fullfile(ecd_dir, 'EEG_MMN.mat'), 'file') == 2;
end

function fixtures = fixture_struct(img_dir, ecd_dir, purpose)
fixtures = struct();
fixtures.img_file = fullfile(img_dir, 'EEG_MMN.mat');
fixtures.ecd_file = fullfile(ecd_dir, 'EEG_MMN.mat');
fixtures.img_dir = img_dir;
fixtures.ecd_dir = ecd_dir;
fixtures.gainmat_file = '';
fixtures.purpose = purpose;
end

function remove_lock(lock_dir)
if exist(lock_dir, 'dir')
    rmdir(lock_dir, 's');
end
end

function copy_gainmatrix(D, out_dir)
try
    gain = D.inv{D.val}.gainmat;
    if isempty(gain), return, end
    [pth, nam, ext] = fileparts(gain);
    if isempty(pth), pth = D.path; end
    src = fullfile(pth, [nam ext]);
    if exist(src, 'file')
        copyfile(src, fullfile(out_dir, [nam ext]));
    end
catch
end
end

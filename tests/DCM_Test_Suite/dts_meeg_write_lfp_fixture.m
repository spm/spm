function D = dts_meeg_write_lfp_fixture(y, dt, outfile, condition_label)
% Write a one-channel LFP waveform as an SPM M/EEG object.
% Authored by Pranay Yadav in 2026

outdir = fileparts(outfile);
if isempty(outdir)
    outdir = pwd;
end
if ~exist(outdir, 'dir')
    mkdir(outdir);
end
[~, name, ext] = fileparts(outfile);
if isempty(ext)
    ext = '.mat';
end

ftdata = [];
ftdata.fsample = 1/dt;
ftdata.label   = {'A1_LFP'};
ftdata.trial   = {y(:)'};
ftdata.time    = {(0:numel(y)-1)*dt};

cwd = pwd;
cd(outdir);
cleanup = onCleanup(@()cd(cwd));

D = spm_eeg_ft2spm(ftdata, [name ext]);
D = chantype(D, ':', 'LFP');
D = units(D, ':', 'a.u.');
D = conditions(D, 1, condition_label);
save(D);
end

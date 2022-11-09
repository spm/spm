function D = spm_eeg_inv_optireg(S)
% D = spm_eeg_inv_optireg(S)
% Registers a template anatomical to SPM M/EEG dataset, using the fiducial
% information obtained from an optical scanning system at the FIL, as a
% part of MEG data collection from Apr 2021.
%
% Input:
%
% S  - input struct
% fields of S:
%
% S.D       - SPM MEEG object                           (REQUIRED)
% S.fidfile - path to .csv file with subject anatomical fidicals and coil
%                locations                              (REQUIRED)
% S.save    - logical to save registration in current dataset
%                                                       (default: TRUE)
% S.forward - calles the forward modelling ui after registration
%                                                       (default: FALSE)
%
% Output:
%
% D - Coregistered dataset.
%__________________________________________________________________________

% George O'Neill
% Copyright (C) 2021-2022 Wellcome Centre for Human Neuroimaging


spm('FnBanner', mfilename);

if ~isfield(S,'D');         error('You must specify a MEEG dataset!');  end
if ~isfield(S,'fidfile');   error('You must specify a fiducial file!'); end
if ~isfield(S,'save');      S.save = 1;                                 end
if ~isfield(S,'forward');   S.forward = 0;                              end


    
% Initialisaion of template anatomical in D
%-------------------------------------------------

[D,val] = spm_eeg_inv_check(S.D);

if val == 0
    val = 1;
end

if ~isfield(D, 'inv') || ~isfield(D.inv{val}, 'comment')
    D.inv = {struct('mesh', [])};
    D.inv{val}.date    = strvcat(date,datestr(now,15));
    D.inv{val}.comment = {''};
else
    tmp = struct('mesh', []);
    tmp.comment = D.inv{val}.comment;
    tmp.date    = D.inv{val}.date;
    D.inv{val} = tmp;
    clear inv
end

D.inv{val}.mesh = spm_eeg_inv_mesh([], 2);

% Load in fiducial file, then prepare intial registration
%----------------------------------------------------------

opti.pnt = csvread(S.fidfile,1,1);
[~, opti.label] = xlsread(S.fidfile,'A2:A7');

targets = {'anat_nas','anat_lpa','anat_rpa'};
for ii = 1:numel(targets)
    id = strmatch(targets{ii}, opti.label, 'exact');
    anatfid.pnt(ii,:) = opti.pnt(id,:);
    anatfid.label{ii} = opti.label{id};
end
anatfid.units = 'mm';
anatfid.fid.pnt = anatfid.pnt;
anatfid.fid.label = anatfid.label;

targets = {'nas','FIL_CTF_L','FIL_CTF_R'};
for ii = 1:numel(targets)
    id = strmatch(targets{ii}, D.inv{val}.mesh.fid.fid.label, 'exact');
    mrifid.pnt(ii,:) = D.inv{val}.mesh.fid.fid.pnt(id,:);
end
mrifid.label = anatfid.label;
mrifid.units = 'mm';
mrifid.fid.pnt = mrifid.pnt;
mrifid.fid.label = mrifid.label;

S1 = [];
S1.sourcefid = anatfid;
S1.targetfid = mrifid;
S1.template = 2;
S1.useheadshape = 0;

fprintf('performing initial registration...')
M1 = spm_eeg_inv_datareg(S1);
fprintf('    COMPLETE!\n')

% Append additional fiducials to template
%----------------------------------------------------------

targets = {'coil_nas','coil_lpa','coil_rpa'};
for ii = 1:numel(targets)
    id = strmatch(targets{ii}, opti.label, 'exact');
    coilfid.pnt(ii,:) = opti.pnt(id,:);
    coilfid.label{ii} = opti.label{id};
end
coilfid.units = 'mm';

coil2mri = ft_transform_geometry(M1,coilfid);

D.inv{val}.mesh.fid.fid.pnt = cat(1,D.inv{val}.mesh.fid.fid.pnt,coil2mri.pnt);
D.inv{val}.mesh.fid.fid.label{end+1} = coil2mri.label{1};
D.inv{val}.mesh.fid.fid.label{end+1} = coil2mri.label{2};
D.inv{val}.mesh.fid.fid.label{end+1} = coil2mri.label{3};

% View to check its going okay so far...
if ~spm('cmdline')
    spm_eeg_inv_checkmeshes(D,val);
end

% Perform the final coregistration between new fiducials in template
% and fidcucial coils as per CTF
%---------------------------------------------------------

ctffid = fiducials(D);
mrifid = [];
targets = {'coil_nas','coil_lpa','coil_rpa'};
for ii = 1:numel(targets)
    id = strmatch(targets{ii}, D.inv{1}.mesh.fid.fid.label, 'exact');
    mrifid.pnt(ii,:) = D.inv{val}.mesh.fid.fid.pnt(id,:);
end
mrifid.label = ctffid.fid.label;
mrifid.fid = mrifid;
mrifid.units = 'mm';

S1 = [];
S1.sourcefid = ctffid;
S1.targetfid = mrifid;
S1.template = 2;
S1.useheadshape = 0;

fprintf('performing final registration...')
M1 = spm_eeg_inv_datareg(S1);
fprintf('      COMPLETE!\n')

ind = 1;
D.inv{val}.datareg(ind).sensors = D.sensors('MEG');
D.inv{val}.datareg(ind).fid_eeg = S1.sourcefid;
D.inv{val}.datareg(ind).fid_mri = ft_transform_geometry(inv(M1), S1.targetfid);
D.inv{val}.datareg(ind).toMNI = D.inv{val}.mesh.Affine*M1;
D.inv{val}.datareg(ind).fromMNI = inv(D.inv{val}.datareg(ind).toMNI);
D.inv{val}.datareg(ind).modality = 'MEG';

% View again!
if ~spm('cmdline')
    spm_eeg_inv_checkdatareg(D,val,ind)
end

%-postable
%--------------------------------------------------------------------------

D = D.history(mfilename, S);

if S.forward
    % Intialise forward gui
    D = spm_eeg_inv_forward_ui(D,val);
end

if S.save
    save(D);
end

%-Cleanup
%--------------------------------------------------------------------------
fprintf('%-40s: %30s\n','Completed',spm('time'));                       %-#
function D = spm_eeg_inv_datareg_ui(varargin)
% Data registration user-interface routine
% commands the EEG/MEG data co-registration within original sMRI space
%
% FORMAT D = spm_eeg_inv_datareg_ui(D,[val], modality)
% Input:
% Output:
% D         - same data struct including the new required files and variables
%__________________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% Vladimir Litvak
% $Id: spm_eeg_inv_datareg_ui.m 2452 2008-11-10 18:45:32Z vladimir $

% initialise
%--------------------------------------------------------------------------
[Finter, Fgraph] = spm('FnUIsetup','MEEG/MRI coregistration', 0);

[D,val] = spm_eeg_inv_check(varargin{:});

try
    D.inv{val}.mesh.template;
catch
    D.inv{val}.mesh.template = 0;
end

if nargin > 2
    modality = varargin{3};
else
    modality = spm_eeg_modality_ui(D, 1);
end

D.inv{val}.modality = modality;

meegfid = D.fiducials;
mrifid = D.inv{val}.datareg.fid_mri;

meeglbl = meegfid.fid.label;
mrilbl = mrifid.fid.label;

newmrifid = mrifid;
newmrifid.fid.pnt = [];
newmrifid.fid.label = {};

if numel(meeglbl)> 3
    [selection ok]= listdlg('ListString', meeglbl, 'SelectionMode', 'multiple',...
        'InitialValue', spm_match_str(upper(meeglbl), upper(mrilbl)), ...
        'Name', 'Select at least 3 fiducials', 'ListSize', [400 300]);

    if ~ok || length(selection) < 3
        error('At least 3 M/EEG fiducials are required for coregistration');
    end

    meegfid.fid.pnt   = meegfid.fid.pnt(selection, :);
    meegfid.fid.label = meegfid.fid.label(selection);
    meeglbl = meeglbl(selection);
end

if numel(meeglbl)<3
    error('At least 3 M/EEG fiducials are required for coregistration');
end

if all(ismember({'spmnas', 'spmlpa', 'spmrpa'}, meegfid.fid.label))
    M1 = eye(4);

    S =[];
    S.sourcefid = meegfid;
    S.targetfid = mrifid;
    S.targetfid.fid = S.sourcefid.fid;
else
    for i = 1:length(meeglbl)
        switch spm_input(['How to specify ' meeglbl{i} ' position?'] , 1, 'select|type|click|skip')
            case 'select'
                [selection ok]= listdlg('ListString', mrilbl, 'SelectionMode', 'single',...
                    'InitialValue', strmatch(upper(meeglbl{i}), upper(mrilbl)), ...
                    'Name', ['Select matching MRI fiducial for ' meeglbl{i}], 'ListSize', [400 300]);
                if ~ok
                    continue
                end

                newmrifid.fid.pnt   = [newmrifid.fid.pnt; mrifid.fid.pnt(selection, :)];
            case 'type'
                pnt = spm_input('Input MNI coordinates', '+1', 'r', '', 3);
                newmrifid.fid.pnt   = [newmrifid.fid.pnt; pnt(:)'];
            case 'click'
                while 1
                    figure(Fgraph); clf;
                    mri = spm_vol(D.inv{val}.mesh.sMRI);
                    spm_orthviews('Reset');
                    spm_orthviews('Image', mri);
                    colormap('gray');
                    cameratoolbar('resetcamera')
                    cameratoolbar('close')
                    rotate3d off;
                    if spm_input(['Select ' meeglbl{i} ' position and click'] , 1,'OK|Retry', [1,0], 1)
                        newmrifid.fid.pnt   = [newmrifid.fid.pnt; spm_orthviews('Pos')'];
                        spm_orthviews('Reset');
                        break;
                    end
                end
            case 'skip'
                meegfid.fid.pnt(i, :) = [];
                meegfid.fid.label(i)  = [];
                continue;
        end
        newmrifid.fid.label = [newmrifid.fid.label  meeglbl{i}];
    end

    if size(newmrifid.fid.label) < 3
        error('At least 3 M/EEG fiducials are required for coregistration');
    end

    % register
    %==========================================================================
    S =[];
    S.sourcefid = meegfid;
    S.targetfid = newmrifid;

    if  D.inv{val}.mesh.template
        if strcmp(D.inv{val}.modality, 'MEG')
            S.template = 2;
        else
            S.template =1;
        end
    else
        S.template = 0;
    end

    if ~isempty(S.sourcefid.pnt)
        S.useheadshape = spm_input('Use headshape points?' , '+1','yes|no', [1,0], 1);
    else
        S.useheadshape = 0;
    end

    %--------------------------------------------------------------------------

    M1 = spm_eeg_inv_datareg(S);
end

switch D.inv{val}.modality
    case 'EEG'
        D.inv{val}.datareg.sensors = forwinv_transform_sens(M1, D.sensors(modality));
        D.inv{val}.datareg.fid_eeg = forwinv_transform_headshape(M1, S.sourcefid);
        D.inv{val}.datareg.fid_mri = S.targetfid;
        D.inv{val}.datareg.fromMNI = inv(D.inv{val}.mesh.Affine);
        D.inv{val}.datareg.toMNI = D.inv{val}.mesh.Affine;
    case 'MEG'
        D.inv{val}.forward.vol = forwinv_transform_vol(inv(M1), D.inv{val}.forward.vol);
        D.inv{val}.datareg.fid_mri = forwinv_transform_headshape(inv(M1), S.targetfid);
        D.inv{val}.mesh = spm_eeg_inv_transform_mesh(inv(M1), D.inv{val}.mesh);
        D.inv{val}.datareg.sensors = D.sensors(modality);
        D.inv{val}.datareg.fid_eeg = S.sourcefid;
        D.inv{val}.datareg.toMNI = D.inv{val}.mesh.Affine*M1;
        D.inv{val}.datareg.fromMNI = inv(D.inv{val}.datareg.toMNI);

end

% check and display registration
%--------------------------------------------------------------------------
spm_eeg_inv_checkdatareg(D);



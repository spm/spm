function D = spm_eeg_ft_datareg_manual(varargin)
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
% $Id: spm_eeg_ft_datareg_manual.m 2389 2008-10-23 11:15:23Z vladimir $

% initialise
%--------------------------------------------------------------------------
[Finter, Fgraph] = spm('FnUIsetup','Fieldtrip MEEG/MRI manual coregistration', 0);

[D,val] = spm_eeg_inv_check(varargin{:});

[ok, D] = check(D, 'sensfid');

if ~ok
    if check(D, 'basic')
        error(['The requested file is not ready for source reconstruction.'...
            'Use prep to specify sensors and fiducials.']);
    else
        error('The meeg file is corrupt or incomplete');
    end
end

try
    D.inv{val}.mesh.template;
catch
    error('Please add a head model to the file before coregistering');
end

if nargin > 2
    modality = varargin{3};
else
    modality = spm_eeg_modality_ui(D, 1);
end

D.inv{val}.modality = modality;

usepolhemus = spm_input('Use polhemus?', 1, 'yes|no', [1, 0]);

if usepolhemus
    meegfid = fileio_read_headshape(spm_select(1, '\.*', 'Select headshape file'));
    meegfid = forwinv_convert_units(meegfid, 'mm');
else
    meegfid = D.fiducials;
end

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

if numel(meeglbl)>=3
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
        
    if spm_input('Choose initial coregistration', 1, 'rigid|align', [1, 0]);
        M1 = spm_eeg_inv_headcoordinates(meegfid.fid.pnt(1, :), meegfid.fid.pnt(2, :), meegfid.fid.pnt(3, :));
        M =  spm_eeg_inv_headcoordinates(newmrifid.fid.pnt(1, :), newmrifid.fid.pnt(2, :), newmrifid.fid.pnt(3, :));
        M1 = inv(M) * M1;
    else
        M1 = spm_eeg_inv_rigidreg(newmrifid.fid.pnt', meegfid.fid.pnt');
    end
    
    meegfid = forwinv_transform_headshape(M1, meegfid);
end
%%
cfg = [];
cfg.individual.headshape = meegfid;
cfg.template.headshape = newmrifid;
cfg = ft_interactiverealign(cfg);

meegfid = forwinv_transform_headshape(cfg.m, meegfid);

M1 = cfg.m * M1;

if usepolhemus
    D.inv{val}.datareg.pol2native = M1;
    origfid = D.fiducials;
    sel1 = spm_match_str(origfid.fid.label, {'nas', 'lpa', 'rpa'});
    sel2 = spm_match_str(meegfid.fid.label, {'nas', 'lpa', 'rpa'});
    M2 = spm_eeg_inv_headcoordinates(meegfid.fid.pnt(sel2(1), :), meegfid.fid.pnt(sel2(2), :), meegfid.fid.pnt(sel2(3), :));
    M3 = spm_eeg_inv_headcoordinates(origfid.fid.pnt(sel1(1), :), origfid.fid.pnt(sel1(2), :), origfid.fid.pnt(sel1(3), :));
    M1 = inv(M2) * M3; 
end
    
switch D.inv{val}.modality
    case 'EEG'
        D.inv{val}.datareg.sensors = forwinv_transform_sens(M1, D.sensors(modality));
        D.inv{val}.datareg.fid_eeg = forwinv_transform_headshape(M1, sensors(D, 'EEG'));
        D.inv{val}.datareg.fid_mri = newmrifid;
        D.inv{val}.datareg.fromMNI = inv(D.inv{val}.mesh.Affine);
        D.inv{val}.datareg.toMNI = D.inv{val}.mesh.Affine;
    case 'MEG'
        D.inv{val}.forward.vol = forwinv_transform_vol(inv(M1), D.inv{val}.forward.vol);
        D.inv{val}.datareg.fid_mri = forwinv_transform_headshape(inv(M1), newmrifid);
        D.inv{val}.mesh = spm_eeg_inv_transform_mesh(inv(M1), D.inv{val}.mesh);
        D.inv{val}.datareg.sensors = D.sensors(modality);
        D.inv{val}.datareg.fid_eeg = forwinv_transform_headshape(inv(M1), meegfid);
        D.inv{val}.datareg.toMNI = D.inv{val}.mesh.Affine*M1;
        D.inv{val}.datareg.fromMNI = inv(D.inv{val}.datareg.toMNI);
end

% check and display registration
%--------------------------------------------------------------------------
spm_eeg_inv_checkdatareg(D);

save(D);
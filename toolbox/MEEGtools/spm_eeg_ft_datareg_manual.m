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
% $Id: spm_eeg_ft_datareg_manual.m 4315 2011-04-26 13:56:07Z vladimir $

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

usepolhemus = spm_input('Use polhemus?', 1, 'yes|no', [1, 0]);

if usepolhemus
    meegfid = ft_read_headshape(spm_select(1, '\.*', 'Select polhemus file'));
    meegfid = ft_convert_units(meegfid, 'mm');
else
    meegfid = D.fiducials;
end

mrifid = D.inv{val}.mesh.fid;

if spm_input('Use individual head surface?', 1, 'yes|no', [1, 0])
    isurf = export(gifti(spm_select(1, 'mesh', 'Select head surface file')), 'ft');
    mrifid.pnt = isurf.pnt;
    mrifid.tri = isurf.tri;
end

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

    switch spm_input('Choose initial coregistration', 1, 'rigid|align1|align2|spm');
        case 'align1'
            M1 = spm_eeg_inv_headcoordinates(meegfid.fid.pnt(1, :), meegfid.fid.pnt(2, :), meegfid.fid.pnt(3, :));
            M =  spm_eeg_inv_headcoordinates(newmrifid.fid.pnt(1, :), newmrifid.fid.pnt(2, :), newmrifid.fid.pnt(3, :));
            M1 = inv(M) * M1;
        case 'align2'
            M1 = spm_eeg_inv_rigidreg(newmrifid.fid.pnt', meegfid.fid.pnt');
            tempfid = ft_transform_headshape(M1, meegfid);
            tempfid.fid.pnt(:, 2) = tempfid.fid.pnt(:, 2)- tempfid.fid.pnt(1, 2)+ newmrifid.fid.pnt(1, 2);
            tempfid.fid.pnt(:, 3) = tempfid.fid.pnt(:, 3)- mean(tempfid.fid.pnt(2:3, 3))+ mean(newmrifid.fid.pnt(2:3, 3));
            M1 = spm_eeg_inv_rigidreg(tempfid.fid.pnt', meegfid.fid.pnt');
        case 'rigid'
            M1 = spm_eeg_inv_rigidreg(newmrifid.fid.pnt', meegfid.fid.pnt');            
        case 'spm'
            S =[];
            S.sourcefid = meegfid;
            S.targetfid = newmrifid;
            S.template  = 0;
            S.useheadshape = ~isempty(S.sourcefid.pnt);
            M1 = spm_eeg_inv_datareg(S);
    end    
    meegfid = ft_transform_headshape(M1, meegfid);
end
%%
cfg = [];
cfg.individual.headshape = meegfid;
cfg.template.headshape = newmrifid;
cfg = ft_interactiverealign(cfg);

meegfid = ft_transform_headshape(cfg.m, meegfid);

M1 = cfg.m * M1;

if usepolhemus
    origfid = D.fiducials;
    sel1 = spm_match_str(origfid.fid.label, {'nas', 'lpa', 'rpa'});
    sel2 = spm_match_str(meegfid.fid.label, {'nas', 'lpa', 'rpa'});
    M2 = spm_eeg_inv_headcoordinates(meegfid.fid.pnt(sel2(1), :), meegfid.fid.pnt(sel2(2), :), meegfid.fid.pnt(sel2(3), :));
    M3 = spm_eeg_inv_headcoordinates(origfid.fid.pnt(sel1(1), :), origfid.fid.pnt(sel1(2), :), origfid.fid.pnt(sel1(3), :));
    M1 = inv(M2) * M3;
end

ind = 1;
D.inv{val}.datareg = struct([]);

if ~isempty(D.sensors('EEG'))
    D.inv{val}.datareg(ind).sensors = ft_transform_sens(M1, D.sensors('EEG'));
    D.inv{val}.datareg(ind).fid_eeg = ft_transform_headshape(M1, D.fiducials);
    D.inv{val}.datareg(ind).fid_mri = newmrifid;
    D.inv{val}.datareg(ind).toMNI = D.inv{val}.mesh.Affine;
    D.inv{val}.datareg(ind).fromMNI = inv(D.inv{val}.datareg(ind).toMNI);
    D.inv{val}.datareg(ind).modality = 'EEG';
    ind = ind+1;
end

if ~isempty(D.sensors('MEG'))
    D.inv{val}.datareg(ind).sensors = D.sensors('MEG');
    D.inv{val}.datareg(ind).fid_eeg = ft_transform_headshape(inv(M1), meegfid);
    D.inv{val}.datareg(ind).fid_mri = ft_transform_headshape(inv(M1), newmrifid);
    D.inv{val}.datareg(ind).toMNI = D.inv{val}.mesh.Affine*M1;
    D.inv{val}.datareg(ind).fromMNI = inv(D.inv{val}.datareg(ind).toMNI);
    D.inv{val}.datareg(ind).modality = 'MEG';
end


% check and display registration
%--------------------------------------------------------------------------
spm_eeg_inv_checkdatareg(D);

save(D);

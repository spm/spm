function D = spm_eeg_inv_datareg_ui(varargin)
% Data registration user-interface routine
% commands the EEG/MEG data co-registration within original sMRI space
%
% FORMAT D = spm_eeg_inv_datareg_ui(D,[val], [meegfid, newmrifid, useheadshape])
% Input:
% D         - M/EEG dataset
%
% meegfid   - M/EEG fiducials
% mrifid    - MRI fiducials
% useheadshape - use headshape points (1)
%
% Output:
% D         - same data struct including the new required files and variables
%__________________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% Vladimir Litvak
% $Id: spm_eeg_inv_datareg_ui.m 4116 2010-11-10 14:23:54Z vladimir $

% initialise
%--------------------------------------------------------------------------
[Finter, Fgraph] = spm('FnUIsetup','MEEG/MRI coregistration', 0);

[D,val] = spm_eeg_inv_check(varargin{:});

mrifid = D.inv{val}.mesh.fid;
mrilbl = mrifid.fid.label;

if nargin>=3
    meegfid = varargin{3};
    interactive = 0;
else
    interactive = 1;
    meegfid = D.fiducials;
    meeglbl = meegfid.fid.label;
    
    if numel(meeglbl)> 3
        [selection ok]= listdlg('ListString', meeglbl, 'SelectionMode', 'multiple',...
            'InitialValue', spm_match_str(upper(meeglbl), upper(mrilbl)), ...
            'Name', 'Select at least 3 fiducials', 'ListSize', [400 300]);
        
        if ~ok || length(selection) < 3
            error('At least 3 M/EEG fiducials are required for coregistration');
        end
        
        meegfid.fid.pnt   = meegfid.fid.pnt(selection, :);
        meegfid.fid.label = meegfid.fid.label(selection);
    end
end

meeglbl = meegfid.fid.label;

if numel(meeglbl)<3
    error('At least 3 M/EEG fiducials are required for coregistration');
end

if all(ismember({'spmnas', 'spmlpa', 'spmrpa'}, meegfid.fid.label)) && isempty(D.sensors('MEG'))
    S =[];
    S.sourcefid = meegfid;
    S.targetfid = mrifid;
    
    if D.inv{val}.mesh.template
        M1 = eye(4);
        S.targetfid.fid = S.sourcefid.fid;
        S.useheadshape = 0;
    else
        M1 = [];
        S.sourcefid.fid.label{strmatch('spmnas', S.sourcefid.fid.label, 'exact')} = 'nas';
        S.sourcefid.fid.label{strmatch('spmlpa', S.sourcefid.fid.label, 'exact')} = 'lpa';
        S.sourcefid.fid.label{strmatch('spmrpa', S.sourcefid.fid.label, 'exact')} = 'rpa';
        S.targetfid.fid.pnt = S.targetfid.fid.pnt(1:3, :);
        S.targetfid.fid.label = S.targetfid.fid.label(1:3, :);
        S.useheadshape = 1;
    end
elseif nargin>=4
    M1 = [];
    S.sourcefid = meegfid;
    S.targetfid = varargin{4};
    
    if nargin >= 5
        S.useheadshape = varargin{5};
    end
else    
    newmrifid = mrifid;
    newmrifid.fid.pnt = [];
    newmrifid.fid.label = {};    
    M1 = [];
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
end

if ~isempty(S.sourcefid.pnt)
    if ~isfield(S, 'useheadshape')
        S.useheadshape = spm_input('Use headshape points?' , '+1','yes|no', [1,0], 1);
    end
else
    S.useheadshape = 0;
end

ind = 1;
D.inv{val}.datareg = struct([]);

[junk, modalities] = modality(D);

if ismember('EEG', modalities) && ~isempty(D.sensors('EEG'))
    if isempty(M1)
        S.template = (D.inv{val}.mesh.template | S.useheadshape);
        M1 = spm_eeg_inv_datareg(S);
    end
    
    D.inv{val}.datareg(ind).sensors = ft_transform_sens(M1, D.sensors('EEG'));
    D.inv{val}.datareg(ind).fid_eeg = ft_transform_headshape(M1, S.sourcefid);
    D.inv{val}.datareg(ind).fid_mri = S.targetfid;
    D.inv{val}.datareg(ind).toMNI = D.inv{val}.mesh.Affine;
    D.inv{val}.datareg(ind).fromMNI = inv(D.inv{val}.datareg(ind).toMNI);
    D.inv{val}.datareg(ind).modality = 'EEG';
    ind = ind+1;
end

if ismember('MEG', modalities) && ~isempty(D.sensors('MEG'))
    if  D.inv{val}.mesh.template
        S.template = 2;
    else
        S.template = 0;
    end
    
    M1 = spm_eeg_inv_datareg(S);
    
    D.inv{val}.datareg(ind).sensors = D.sensors('MEG');
    D.inv{val}.datareg(ind).fid_eeg = S.sourcefid;
    D.inv{val}.datareg(ind).fid_mri = ft_transform_headshape(inv(M1), S.targetfid);
    D.inv{val}.datareg(ind).toMNI = D.inv{val}.mesh.Affine*M1;
    D.inv{val}.datareg(ind).fromMNI = inv(D.inv{val}.datareg(ind).toMNI);
    D.inv{val}.datareg(ind).modality = 'MEG';
end

% check and display registration
%--------------------------------------------------------------------------
if interactive
    spm_eeg_inv_checkdatareg(D);
else
   for i = 1:numel(D.inv{val}.datareg)
       spm_eeg_inv_checkdatareg(D, val, i);
   end
end



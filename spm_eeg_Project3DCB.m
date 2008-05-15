function [xy,label] = spm_eeg_Project3DCB(D, datatype)
% Wrapper function to a fieldtrip function to project 3D locations 
% onto a 2D plane. 
% FORMAT [xy,label] = spm_eeg_Project3DCB(D, datatype)
% _______________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% Stefan Kiebel, Vladimir Litvak
% $Id: spm_eeg_Project3DCB.m 1668 2008-05-15 17:50:41Z vladimir $

cfg = [];

switch datatype
    case 'EEG'
        cfg.elec = D.sensors('EEG');
        label = cfg.elec.label;
        cfg.rotate = 0;
    case 'MEG'
       cfg.grad = D.sensors('MEG');
       label = cfg.grad.label;
    otherwise
        error('Unknown data type');
end

lay = ft_prepare_layout(cfg);
[sel1, sel2] = spm_match_str(label, lay.label);

label =lay.label(sel2)';
xy = lay.pos(sel2, :);

nchan = size(xy, 1);

xy =(xy-repmat(min(xy), nchan, 1));
xy = xy./repmat(max(xy), nchan, 1);
xy = xy*0.9+0.05;
xy = xy';



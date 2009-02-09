function [xy,label] = spm_eeg_layout3D(sens, modality)
% Wrapper function to a fieldtrip function to project 3D locations 
% onto a 2D plane. 
% FORMAT [xy,label] = spm_eeg_project3D(sens, modality)
% _______________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% Vladimir Litvak
% $Id: spm_eeg_layout3D.m 2720 2009-02-09 19:50:46Z vladimir $


switch modality
    case 'EEG'
        xy = sens.pnt;
        label = sens.label;
    case 'MEG'
        cfg = [];
        cfg.style = '3d';
        cfg.rotate = 0;
        cfg.grad = sens;

        lay = ft_prepare_layout(cfg);

        [sel1, sel2] = spm_match_str(sens.label, lay.label);

        label = lay.label(sel2, 1);
        xy = lay.pos(sel2, :);
end


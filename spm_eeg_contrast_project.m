function c = spm_eeg_contrast_project(SPM, i, c_in)
% function to generate a contrast matrix for factor i.
% FORMAT c = spm_eeg_contrast_project(SPM, i, c_in)
%
% SPM  - SPM struct
% i    - index of factor
% c_in - vector/matrix in observation space
%
% c    - contrast vector/matrix 
%_______________________________________________________________________
%
% spm_eeg_contrast_project is an internally used function that takes
% vector/matrices in observation space and transforms them to contrast
% vector/matrices. See also [1], Eq. ??. 
% To generate contrast weights, the function has to, if necessary, iterate
% through subordinate design matrices of a hierarchical model.
%
% [1] S.J. Kiebel and K.J. Friston (2004). Statistical Parametric Mapping
% for event-related potentials (II): a hierarchical temporal model.
% Neuroimage, 22: 503-520
%_______________________________________________________________________
% Copyright (C) 2005 Wellcome Department of Imaging Neuroscience

% Stefan Kiebel
% $Id: spm_eeg_contrast_project.m 112 2005-05-04 18:20:52Z john $


Fname = SPM.eeg.factor{i};
Ilevel = SPM.eeg.Ilevel;

if strcmpi(SPM.xBF.name_d(1,i), 'Identity') & Ilevel > 1
    % These are parameters taken up from a lower level and we have 
    % to get at the design matrix component at a lower level. 
    % We assume that all parameters have been taken up, not only some
    % contrast(s) of them. (This assumption might be relaxed in future
    % versions.)
    
    
    SPM1 = SPM;
    ind = [];
    while Ilevel > 1 & isempty(ind)
        SPM1 = load(SPM1.eeg.subordinate);
        SPM1 = SPM1.SPM;
        for k = 1:SPM1.eeg.Ncomp_d
            for j = 1:SPM1.eeg.Nfactors
                if strcmpi(SPM1.eeg.factor{j}, Fname) & ~strcmpi(SPM1.xBF.name_d{k, j}, 'Identity') & ~strcmpi(SPM1.xBF.name_d{k, j}, 'Constant')
                    ind = [k j];
                end
            end
        end
        Ilevel = Ilevel-1;
    end
    
    % first-level time-design matrix
    if isempty(ind)
        error('Couldn''t get design matrix component for factor %s', Fname);
    end
    
    X = SPM1.eeg.X_d{ind(1), ind(2)};
else
    X = SPM.eeg.X_d{1, i};
end

% pad with zeros if necessary
if size(c_in, 2) < size(X, 1)
    c_in = [c_in zeros(size(c_in,1), size(X, 1) - size(c_in, 2))];
elseif size(c_in, 2) > size(X, 1)
    warning('Wrong vector length for factor %s', SPM.eeg.factor{i});
    w = 1;
end

c = (pinv(c_in')*X);


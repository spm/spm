function SPM = spm_eeg_get_vc(SPM);
% generate variance components for EEG design
% FORMAT SPM = spm_eeg_get_vc(SPM);
%
% SPM   - SPM struct
% in SPM.xVi:
%    Qidentical   - cell vector of 0/1-values that determines whether
%                   factors have identical variances or not
%                   
%    Qindependent - cell vector of 0/1-values that determines whether
%                   factors have correlated error or not
%
% SPM.xVi.V       - cell vector of covariance components
%_______________________________________________________________________
% 
% spm_eeg_get_vc generates variance components for a given design. For each
% factor, the user specifies whether its levels have identical variances
% and are uncorrelated. The individual components for each factor are
% combined into covariance components by using the kronecker tensor
% product, if one has the same number of observations under all levels and
% factors. If there are unequal number of observations at different levels,
% the function uses concatenation of block matrices instead.
% 
% The functionality of spm_eeg_get_vc is similar to that of
% spm_non_sphericity. The difference is that spm_eeg_get_vc can accomodate
% any number of factors and is more general, because it can cope with
% different number of observations under different levels of a factor.
%_______________________________________________________________________
% Copyright (C) 2005 Wellcome Department of Imaging Neuroscience

% Stefan Kiebel
% $Id: spm_eeg_get_vc.m 213 2005-08-22 12:43:29Z stefan $



% Cycle through design components
% There are only 2 design decisions to make for each design component:
% (i) Identical variances yes/no and (ii) independence yes/no

for i = 1:SPM.eeg.Nfactors
    try
        Qidentical{i} = SPM.xVi.Qidentical{i};
    catch
        Qidentical{i} = spm_input(sprintf('Identical variance for factor %s?', SPM.eeg.factor{i}), '+1', 'yes|no', [1, 0]);
        SPM.xVi.Qidentical{i} = Qidentical{i};
    end
    try
        Qindependent{i} = SPM.xVi.Qindependent{i};
    catch
        Qindependent{i} = spm_input(sprintf('Independence for factor %s?', SPM.eeg.factor{i}), '+1', 'yes|no', [1, 0]);
        SPM.xVi.Qindependent{i} = Qindependent{i};
        
    end
end

% create components for each factor
%-----------------------------------------------------------------------
Q = {};

% for each factor
for i = 1:SPM.eeg.Nfactors

    % identical variances
    if Qidentical{i}
        
        % if there is the same number of levels over all levels of the
        % supraordinate factor
        if length(SPM.eeg.Nlevels{i}) == 1
            Lmax = SPM.eeg.Nlevels{i};    
            % otherwise 
        else
            Lmax = max(SPM.eeg.Nlevels{i});
        end
    
        Q{i}{1} = speye(Lmax, Lmax);
        
        % non-identical variances
    else
        
        % if there is the same number of levels over all levels of the
        % supraordinate factor        
        if length(SPM.eeg.Nlevels{i}) == 1
            Lmax = SPM.eeg.Nlevels{i};
        else
            Lmax = max(SPM.eeg.Nlevels{i});
        end
        
        for j = 1:Lmax
            Q{i}{j} = sparse(j, j, 1, Lmax, Lmax);
        end
    end
    
    if ~Qindependent{i}
        if length(SPM.eeg.Nlevels{i}) == 1
            Lmax = SPM.eeg.Nlevels{i};
        else
            Lmax = max(SPM.eeg.Nlevels{i});
        end
        
        for j = 1:Lmax
            for k = j+1:Lmax
                Q{i}{end+1} = sparse([j k], [k j], 1, Lmax, Lmax);
            end
        end
    end
end

% generate variance components
%-----------------------------------------------------------------------
Xind = SPM.eeg.Xind;

tmp = {1}; 

% loop over factors
for i = 1:SPM.eeg.Nfactors
    V = {};
    
    % loop over number of variance components per factor
    for j = 1:length(Q{i})
        
        % loop over number of already existing variance components
        for k = 1:length(tmp)
            
            % if all levels equal at this level
            if length(SPM.eeg.Nlevels{i}) == 1
                V{end+1} = kron(tmp{k}, Q{i}{j});
                
                % otherwise
            else
                Mlevel = max(SPM.eeg.Nlevels{i});
                
                tmp2 = [];
                Nlevels = SPM.eeg.Nlevels{i};
                for l = 1:length(Nlevels)
                    
                    % is there one or more missing levels?
                    if SPM.eeg.Nlevels{i}(l) < Mlevel
                        % find the existing indices from Xind
                        L = Xind{i-1}(l, :);
                        ind = find(all(Xind{i}(:, 1:end-1) == kron(ones(size(Xind{i},1), 1), L), 2));
                        ind_1 = Xind{i}(ind, end);
                    else
                        ind_1 = [1:Mlevel];
                    end
                    
                    tmp3 = [];
                    for m = 1:length(Nlevels)
                        if l == m
                            % variance
                            tmp3 = [tmp3 tmp{k}(l, m)*Q{i}{j}(ind_1, ind_1)];
                        else
                            % covariance
                            % is there one or more missing levels?
                            if SPM.eeg.Nlevels{i}(m) < Mlevel
                                % find the existing indices from Xind
                                L = Xind{i-1}(m, :);
                                ind = find(all(Xind{i}(:, 1:end-1) == kron(ones(size(Xind{i},1), 1), L), 2));
                                ind_2 = Xind{i}(ind, end);
                            else
                                ind_2 = [1:Mlevel];
                            end
                            
                            tmp3 = [tmp3 tmp{k}(l, m)*Q{i}{j}(ind_1, ind_2)];
                        end
                    end
                    tmp2 = [tmp2; tmp3];
                end
                V{end+1} = tmp2;
            end
        end    
    end
    tmp = V; 
end

% remove all-zero variance components
ind = [];
for i = 1:length(tmp)
    if sum(tmp{i}(:)) == 0
        ind = [ind i];
    end
end
tmp(ind) = [];

SPM.xVi.Vi = tmp;



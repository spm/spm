function [SPM] = spm_eeg_design(SPM)
% Assembles hierarchical linear model for EEG/MEG data
% FORMAT [SPM] = spm_eeg_design(SPM)
%_______________________________________________________________________
% 
% spm_eeg_design is a comprehensive and general routine to assmemble all
% information needed to specify a hierarchical linear model. Each level of
% the hierarchy is specified by a separate call to spm_eeg_design. At each
% level, one can specify a model with an arbitrary number of levels. The
% different levels of the hierarchy are linked by an entry in the SPM
% struct (SPM.eeg.subordinate) that stores the location of the
% 'subordinate' SPM file.
% At each level, one specifies the model as components, where each
% component is associated with a factor. Examples for factors are group,
% subject, trial type, single trial, time. Each model component is fully
% specified by providing a design component and a variance component. At
% the inference stage, contrast weights are formed in a similar fashion.
%_______________________________________________________________________
% Stefan Kiebel & Karl Friston $Id$

%-GUI setup
%-----------------------------------------------------------------------
% SPMid = spm('SFnBanner',mfilename,SCCSid);
[Finter,Fgraph,CmdLine] = spm('FnUIsetup', 'EEG stats model setup',0);
spm_help('!ContextHelp', mfilename)

% is this the first level or a subordinate level?
%-----------------------------------------------------------------------
try
    Ilevel        = SPM.eeg.Ilevel;
catch
    Ilevel        = spm_input('Is this a first-level design?', '+1', 'yes|no', [1 0]);    
end

if Ilevel ~= 1
    try
        SPM.eeg.subordinate;
    catch
        SPM.eeg.subordinate = spm_get(1, 'SPM.mat', 'Select SPM.mat of subordinate level');
    end
    
    S = load(SPM.eeg.subordinate);
    SPM.eeg.Ilevel = S.SPM.eeg.Ilevel+1;
else
    SPM.eeg.Ilevel = Ilevel;
end

% construct Design matrix {X}
%=======================================================================

% get sampling frequency
%-----------------------------------------------------------------------
if Ilevel == 1
    try
        SPM.xY.RT;
    catch
        spm_input('Basic parameters...',1,'d',mfilename)
        SPM.xY.RT = 1/spm_input('Sampling frequency (Hz)','+1','r',[],1);
    end
else
    try
        SPM.eeg.RT = S.SPM.eeg.RT;
    catch
        SPM.eeg.RT = S.SPM.xY.RT;
    end
end

% get peri-stimulus times of epoch
%-----------------------------------------------------------------------
if Ilevel == 1
    try
        SPM.eeg.pt;
    catch
        pt        = spm_input('Peri-stimulus times (ms)','+1', 'r', '-95 1180', 2);
        SPM.eeg.pt = pt;
    end
else
    SPM.eeg.pt = S.SPM.eeg.pt;
end

% get number of factors
try
	SPM.eeg.Nfactors;
catch
	SPM.eeg.Nfactors = spm_input(['How many factors?'], '+1');
end

for i = 1:SPM.eeg.Nfactors
    
    % get factor name
    try
        SPM.eeg.factor{i};
    catch
        SPM.eeg.factor{i} = spm_input(sprintf('Factor name %d', i),'+1', 's');
    end
    
    % get number of levels
    try
        SPM.eeg.Nlevels{i};
    catch
        % number of level can be either scalar or vector. If scalar, all
        % levels of the supraordinate factor have the same number of
        % levels, otherwise the elements of the vector list all number of
        % levels explicitly.
        
        % recursive computation of number of levels that we need at this
        % level
        if i == 1
            Ntotal(1) = 1;
        else
            if length(SPM.eeg.Nlevels{i-1}) == 1
                Ntotal(i) = Ntotal(i-1)*SPM.eeg.Nlevels{i-1};
            else
                Ntotal(i) = sum(SPM.eeg.Nlevels{i-1});
            end
        end
        
        str = sprintf('How many levels for factor %s', SPM.eeg.factor{i});
        Ypos = -1;
        while 1
            if Ypos == -1   
                [SPM.eeg.Nlevels{i}, Ypos] = spm_input(str, '+1', 'n');
            else
                SPM.eeg.Nlevels{i} = spm_input(sprintf(str, SPM.eeg.factor{i}), Ypos, 'n');
            end

            if length(SPM.eeg.Nlevels{i}) == 1 | length(SPM.eeg.Nlevels{i}) == Ntotal(i), break, end
            str = sprintf('enter a scalar or [%d] vector', Ntotal(i));
        end
    end
end

% time units are always ms!
%-----------------------------------------------------------------------
SPM.xBF.UNITS = 'millisecs';

% separate specifications for non-replicated sessions
%-----------------------------------------------------------------------
% assume for now that we have the same design matrix for each ERP
% likely exceptions: different number of single trials per trial type and
% different number of subjects per group.

% How many design components
%-----------------------------------------------------------------------
try
	SPM.eeg.Ncomp_d;
catch
	SPM.eeg.Ncomp_d = spm_input(['How many design partitions?'], '+1');
end

try
    SPM.xBF.name_d;
catch
    for i = 1:SPM.eeg.Ncomp_d
        for j = 1:SPM.eeg.Nfactors
            spm_input(sprintf('Partition %d: design component %s', i, SPM.eeg.factor{j}), 1, 'd');
            Ctype = {
                'Identity',...
                    'Constant',...
                    'Wavelets Db4',...
                    'User defined'};
            str   = 'Select design component';
            Sel   = spm_input(str, 2, 'm', Ctype);
            SPM.xBF.name_d{i, j} = Ctype{Sel};
        end        
    end
end

% get parameters for design components
% We assume the same parameters for each design partition
%------------------------------------------------------
for j = 1:SPM.eeg.Nfactors
    switch SPM.xBF.name_d{1, j}
        case {'Wavelets Db4'}
            try
                SPM.xBF.order_d(j);
            catch
                if length(unique(SPM.eeg.Nlevels{j})) > 1
                    error('Factor %s must have equal number of levels!', SPM.eeg.factor{j});
                else
                    SPM.xBF.order_d(j) = spm_input('order', '+1', 'n', log2(SPM.eeg.Nlevels{j}));
                end
            end
    end
end

for j = 1:SPM.eeg.Nfactors
    switch SPM.xBF.name_d{1, j}
        % get parameters for design components
        case {'Wavelets Db4'}
            try	
                SPM.xBF.truncation_d(j);
            catch
                SPM.xBF.truncation_d(j) = spm_input('truncation', '+1', 'n', 2);
            end
    end
end

% design 'anomalies', i.e. missing cells
%------------------------------------------------------
Xind = {};
for i = 1:SPM.eeg.Nfactors
    
    tmp = [];
    
    % build matrix Xind with level indices
    %-------------------------------------
    if i > 1 & length(SPM.eeg.Nlevels{i}) == 1
        Ilevels = SPM.eeg.Nlevels{i}*ones(size(Xind{i-1}, 1), 1);
    else
        Ilevels = SPM.eeg.Nlevels{i};
    end
    
    Mlevel = max(Ilevels);

    % loop through factors and ask for level indices of missing cells
    %----------------------------------------
    for j = 1:length(Ilevels)
        if i == 1
            tmp = [1:Ilevels(j)]';
        else
            try
                m = SPM.eeg.Imissing{i}(j);
            catch
                if Ilevels(j) < Mlevel
                    str = '';
                    for k = 1:i-1
                        str = [str sprintf('%s (%d), ', SPM.eeg.factor{k}, Xind{i-1}(j, k))];
                    end
                    if Mlevel-Ilevels(j) > 1                        
                        str = [str 'missing levels'];
                    else
                        str = [str 'missing level'];
                    end
                    m = spm_input(str, '+1', 'n', num2str([Ilevels(j)+1:Mlevel]), Mlevel-Ilevels(j), Mlevel);              
                else
                    m = [];
                end
            end
            
            I = [1:Mlevel]';
            I(m) = [];
            tmp = [tmp; [kron(Xind{i-1}(j,:), ones(Ilevels(j),1)) I]];
        end
    end
    Xind{i} = tmp;
end

SPM.eeg.Xind = Xind;


% generate components
%---------------------------------------------------------------
Xn = {};
for j = 1:SPM.eeg.Nfactors
    for i = 1:SPM.eeg.Ncomp_d
        switch SPM.xBF.name_d{i, j}
            case{'Identity'}
                SPM.eeg.X_d{i,j} = speye(max(SPM.eeg.Nlevels{j}));
                for k = 1:max(SPM.eeg.Nlevels{j})
                    Xn{i, j}{k} = sprintf('%s %d', SPM.eeg.factor{j}, k);
                end
                
            case{'Constant'}
                SPM.eeg.X_d{i,j} = ones(max(SPM.eeg.Nlevels{j}), 1);
                Xn{i, j}{1} = sprintf('avg over %s', SPM.eeg.factor{j});
            case{'Wavelets Db4'}
                [Xt, Nres] = make_wavelets(SPM.eeg.Nlevels{j}, 'db4', SPM.xBF.order_d(j));
                % truncation
                Xt = Xt(:, 1:Nres(end-SPM.xBF.truncation_d(j)));
                SPM.xBF.Nres{i, j} = Nres(1:end-SPM.xBF.truncation_d(j));
                SPM.eeg.X_d{i, j} = Xt;
                for k = 1:size(SPM.eeg.X_d{i,j}, 2)
                    Xn{i, j}{k} = sprintf('bf_W %d', k);
                end
            case{'User defined'}
                % User defined regressors have different partition for each
                % combination of supra-ordinate factors
                try 
                    SPM.eeg.X_d{i, j};
                catch
                    % GUI not functional yet!!
                    for k = 1:2
                        % get design partitions, each input one matrix
                        SPM.eeg.X_d{i,j} = spm_input('UD for ', '+1', 'w1', 0);
                    end    
                end
                
                % still to do: build in error checking of design partitions
               
                for k = 1:size(SPM.eeg.X_d{i,j}, 2)
                    Xn{i, j}{k} = sprintf('bf_U %d', k);
                end
                
            otherwise
                error('Unknown design component');
        end
    end        
end

% generate design matrix
X = [];
for i = 1:SPM.eeg.Ncomp_d
    tmp = 1;
    for j = 1:SPM.eeg.Nfactors
        if length(SPM.eeg.Nlevels{j}) == 1
            tmp = kron(tmp, SPM.eeg.X_d{i, j});
        else
            tmp2 = [];

            for k = 1:length(SPM.eeg.Nlevels{j})
                
                % extract indices
                if ~strcmpi(SPM.xBF.name_d{i, j}, 'User defined')
                    Lmax = max(SPM.eeg.Nlevels{j});
                    if SPM.eeg.Nlevels{j}(k) < Lmax
                        L = Xind{j-1}(k,:);
                        ind = find(all(Xind{j}(:, 1:end-1) == kron(ones(size(Xind{j},1), 1), L), 2));
                        ind_1 = Xind{j}(ind, end);
                    else
                        ind_1 = 1:Lmax;
                    end
                end
          
                tmp3 = [];
                for l = 1:size(tmp,2)

                    % use switch 
                    switch SPM.xBF.name_d{i, j}
                        case{'Identity'}
                            % identity special case, because one regressor per
                            % level
                            Lmax = max(SPM.eeg.Nlevels{j});
                            if SPM.eeg.Nlevels{j}(l) < Lmax
                                
                                % indices of different levels in Xind
                                L = Xind{j-1}(l,:);
                                ind = find(all(Xind{j}(:, 1:end-1) == kron(ones(size(Xind{j},1), 1), L), 2));
                                
                                % are there other levels?
                                if j < 3
                                    ind2 = [1:size(Xind{j}, j)]';
                                else
                                    ind2 = Xind{j}(find(Xind{j}(:, j-2)==Xind{j}(j, 1)), j);
                                end
                                ind_2 = unique([ind2; Xind{j}(ind, end)]);
                                
                            else
                                ind_2 = 1:Lmax;
                            end
                            Xpart = SPM.eeg.X_d{i, j}(ind_1, ind_2);
                        case{'Constant', 'Wavelets Db4'}
                            ind_2 = 1:size(SPM.eeg.X_d{i, j}, 2);
                            Xpart = SPM.eeg.X_d{i, j}(ind_1, ind_2);
                        case{'User defined'}
                            Xpart = SPM.eeg.X_d{i, j}{k};
                        otherwise
                            error('Unknown design partition type');                            
                    end
                    tmp3 = [tmp3 tmp(k,l)*Xpart];
                    
                end
                tmp2 = [tmp2; tmp3];
            end
            tmp = tmp2;
        end
    end
    X = [X tmp];
end

SPM.xX.X = X;

% construct Xname cell vector
Xname = {};
if SPM.eeg.Nfactors == 1
    for i = 1:SPM.eeg.Ncomp_d
        Xname = [Xname Xn{1,1}];
    end
else
    for i = 1:SPM.eeg.Ncomp_d    
        Xname_d = {};
        
        for j = 1:SPM.eeg.Nfactors
            tmp = {};
            
            for k = 1:max(1, length(Xname_d))
                for l = 1:length(Xn{i, j})
                    if isempty(Xname_d)
                        tmp = [tmp Xn{i,j}(l)];
                    else
                        tmp = [tmp {[Xname_d{k} ' x ' Xn{i,j}{l}]}];
                    end
                end
            end
            Xname_d = tmp;
        end    
        Xname = [Xname Xname_d];
    end
end

% remove all-zero regressors and save indices for
% contrast vector/matrix generation
ind = find(all(~SPM.xX.X));
SPM.xX.X(:, ind) = [];
Xname(ind) = [];
SPM.eeg.remove = ind;

% finished
%-----------------------------------------------------------------------
SPM.xX.iH     = [];
SPM.xX.iC     = [1:size(X,2)];
SPM.xX.iB     = [];
SPM.xX.iG     = [];
SPM.xX.name   = {Xname{:}};

% get variance components
%===============================================================
SPM = spm_eeg_get_vc(SPM);


%-End: Save SPM.mat
%-----------------------------------------------------------------------
fprintf('\t%-32s: ','Saving ERP design')                            %-#
fprintf('%30s\n','...done')                                 %-#
save SPM SPM;

spm_input('!DeleteInputObj')


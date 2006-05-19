function [SPM] = spm_eeg_design(SPM)
% Assembles hierarchical linear model for EEG/MEG data
% FORMAT [SPM] = spm_eeg_design(SPM)
%_______________________________________________________________________
% 
% spm_eeg_design is a general routine to assmemble all
% information needed to specify a hierarchical linear model. Each level of
% the hierarchy is specified by a separate call to spm_eeg_design.
% At each level, one specifies the model as components, where each
% component is associated with a factor. Examples for factors are group,
% subject, trial type, single trial, time. Each model component is fully
% specified by providing design components and variance components. At
% the inference stage, contrast weights are formed in a similar fashion.
%_______________________________________________________________________
% Copyright (C) 2005 Wellcome Department of Imaging Neuroscience

% Stefan Kiebel
% $Id: spm_eeg_design.m 539 2006-05-19 17:59:30Z Darren $

%-GUI setup
%-----------------------------------------------------------------------
[Finter,Fgraph,CmdLine] = spm('FnUIsetup', 'EEG stats model setup',0);
spm_help('!ContextHelp', mfilename)

Oanalysis = SPM.eeg.Oanalysis;

% set all options if shortcut
if Oanalysis == 1
    % ERP analysis
    SPM.eeg.Ilevel = 1;
    SPM.eeg.Nfactors = 2;
    SPM.eeg.factor{1} = 'conditions';
    SPM.eeg.factor{2} = 'time';
    SPM.eeg.Ncomp_d = 1;
    SPM.xBF.name_d{1, 1} = 'Identity';
    SPM.xBF.name_d{1, 2} = 'Identity';    
    SPM.xVi.Qidentical{1} = 1;
    SPM.xVi.Qindependent{1} = 1;
    SPM.xVi.Qidentical{2} = 1;
    SPM.xVi.Qindependent{2} = 1;
    SPM.xX.fullrank = 1;
end

% is this the first level?
%-----------------------------------------------------------------------
try
    Ilevel = SPM.eeg.Ilevel;
catch
    Ilevel = spm_input('Is this a first-level design?', '+1', 'yes|no', [1 0]);
    SPM.eeg.Ilevel = Ilevel;
end

if Ilevel == 1
    try
        D = spm_eeg_ldata(SPM.eeg.D);
    catch
        % ask for generating ERP file to get at timing parameters
        D = spm_eeg_ldata(spm_select(1, 'mat', 'Select one of the original M/EEG-mat file(s)'));
    end
    SPM.xY.RT = 1/D.Radc;
    SPM.eeg.pt = 1000/D.Radc*[-D.events.start D.events.stop];
end

% construct Design matrix {X}
%=======================================================================

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
        
        str = sprintf('#levels for factor %s', SPM.eeg.factor{i});
        Ypos = -1;
        while 1
            if Ypos == -1
                if strcmpi(SPM.eeg.factor{i}, 'time')
                    if ~isempty(who('D'))
                        defstr = sprintf('%d', D.Nsamples); 
                    end
                else
                    defstr = '';
                end
                [SPM.eeg.Nlevels{i}, Ypos] = spm_input(str, '+1', 'n', defstr);
            else
                SPM.eeg.Nlevels{i} = spm_input(sprintf(str, SPM.eeg.factor{i}), Ypos, 'n', defstr);
            end
            
            if length(SPM.eeg.Nlevels{i}) == 1 | length(SPM.eeg.Nlevels{i}) == Ntotal(i), break, end
            str = sprintf('enter a scalar or [%d] vector', Ntotal(i));
        end
    end
end

% time units are always ms!
%-----------------------------------------------------------------------
SPM.xBF.UNITS = 'millisecs';

% How many design components
%-----------------------------------------------------------------------
try
	SPM.eeg.Ncomp_d;
catch
    if Ilevel ~= 1
        SPM.eeg.Ncomp_d = spm_input(['How many design partitions?'], '+1');
    else
        SPM.eeg.Ncomp_d = 1;
    end
end

try
    SPM.xBF.name_d;
catch
    for i = 1:SPM.eeg.Ncomp_d
        for j = 1:SPM.eeg.Nfactors
            spm_input(sprintf('Partition %d: design component %s', i, SPM.eeg.factor{j}), 1, 'd');
            Ctype = {
                'Identity',...
                'Constant'};
            str = 'Select design component';
            Sel = spm_input(str, 2, 'm', Ctype);
            SPM.xBF.name_d{i, j} = Ctype{Sel};
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

    % loop over factors and ask for level indices of missing cells
    %-------------------------------------------------------------
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

% special case of a single factor (e.g. 1-sample test)
if length(Xind) == 1 
    Xind{2} = 1;
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
                            Xpart = SPM.eeg.X_d{i, j}(ind_2, ind_2);
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
ind = find(~any(SPM.xX.X));
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

if size(SPM.xX.X, 1) == size(SPM.xX.X, 2)
    % i.e. we cannot estimate any error
    SPM.xX.fullrank = 1;
    for i = 1 : SPM.eeg.Nfactors
        SPM.xVi.Qidentical{i} = 1;
        SPM.xVi.Qindependent{i} = 1;
    end
end

% get variance components
%===============================================================
SPM = spm_eeg_get_vc(SPM);


%-End: Save SPM.mat
%-----------------------------------------------------------------------
fprintf('\t%-32s: ','Saving ERP design')                            %-#
if spm_matlab_version_chk('7') >= 0
    save('SPM', 'SPM', '-V6');
else
    save('SPM', 'SPM');
end;
fprintf('%30s\n','...done')                                 %-#

spm_input('!DeleteInputObj')


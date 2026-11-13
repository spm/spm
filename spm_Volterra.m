function [X,Xname,Fc] = spm_Volterra(U,bf,V)
% Generalized convolution of inputs (U) with basis set (bf)
% FORMAT [X,Xname,Fc] = spm_Volterra(U,bf,V)
% U          -  input structure array (see spm_get_ons.m)
% bf         -  Basis functions (see spm_get_bf.m)
% V          -  [1 or 2] order of Volterra expansion [default = 1]
%
% X          -  Design Matrix
% Xname      -  names of regressors [columns] in X
% Fc(i).i    -  indices pertaining to input i (and interactions)
% Fc(i).name -  names pertaining to input i   (and interactions)
% Fc(i).p    -  grouping of regressors per parameter
%__________________________________________________________________________
%
% For first order expansions spm_Volterra simply convolves the causes (e.g.
% stick functions) in U.u by the basis functions in bf to create a design
% matrix X.  For second order expansions new entries appear in X, Xname and
% Fc that correspond to the interaction among the original causes. The
% basis functions for these effects are two dimensional and are used to
% assemble the second order kernel in spm_graph.m. Second order effects are
% computed for only the first column of U.u.
%__________________________________________________________________________

% Karl Friston
% Copyright (C) 1999-2022 Wellcome Centre for Human Neuroimaging


%-Order of Volterra expansion
%--------------------------------------------------------------------------
if nargin == 2, V = 1; end

%-Construct X
%--------------------------------------------------------------------------

%-1st order terms
%==========================================================================
X     = [];
Xname = {};
Fc    = [];
for i = 1:numel(U)
    ind   = [];
    ip    = [];
    for k = 1:size(U(i).u,2)
        for p = 1:size(bf,2)
            x = U(i).u(:,k);
            d = 1:length(x);
            x = conv(full(x),bf(:,p));
            x = x(d);
            X = [X x];
            
            %-Indices and regressor names
            %------------------------------------------------------------------
            Xname{end + 1} = sprintf('%s*bf(%i)',U(i).name{k},p);
            ind(end + 1)   = size(X,2);
            ip(end + 1)    = k;
        end
    end
    Fc(end + 1).i = ind;
    Fc(end).name  = U(i).name{1};
    Fc(end).p     = ip;
end

%-Return if first order
%--------------------------------------------------------------------------
if V == 1, return, end

%-2nd order terms
%==========================================================================

% List all permutations of conditions and basis functions
K = [];
for i = 1:numel(U)
    for p = 1:size(bf,2)
       K(end+1,:) = [i, p]; % condition, bf
    end
end

% Initialize array to store indices of design matrix columns involved in 
% each 2-way interaction
Xcond = cell(numel(U),numel(U));

% Compute all 2-way interactions
for i = 1:size(K,1)
    for j = i:size(K,1)        
        % Condition indices
        condx = K(i,1);
        condy = K(j,1);
        
        % Basis function indices
        bfx = K(i,2);
        bfy = K(j,2);
        
        % Convolve HRF then compute interaction
        x = U(condx).u(:,1);
        y = U(condy).u(:,1);
        x = conv(full(x),bf(:,bfx));
        y = conv(full(y),bf(:,bfy));
        x = x(d);
        y = y(d);
        X = [X x.*y];
        
        % Regressor names
        Xname{end + 1} = sprintf('%s*bf(%i)x%s*bf(%i)',...
            U(condx).name{1}, bfx,...
            U(condy).name{1}, bfy);
        
        % Design matrix columns per [condition x condition] permutation
        if isempty(Xcond{condx,condy})
            Xcond{condx,condy} = size(X,2);
        else
            Xcond{condx,condy}(end+1) = size(X,2);
        end
    end
end

% Set metadata about each new interaction condition
for i = 1:numel(U)
    for j = 1:numel(U)
        if ~isempty(Xcond{i,j})
            Fc(end+1).i  = Xcond{i,j};
            Fc(end).name = [U(i).name{1} 'x' U(j).name{1}];
            Fc(end).p    = ones(1,numel(Xcond{i,j}));
        end
    end
end
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
% Copyright (C) 1999-2012 Wellcome Trust Centre for Neuroimaging

% Karl Friston
% $Id: spm_Volterra.m 4925 2012-09-14 11:17:01Z guillaume $


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
for i = 1:numel(U) 
for j = i:numel(U)
    ind   = [];
    ip    = [];
    for p = 1:size(bf,2)
    for q = 1:size(bf,2)
        x = U(i).u(:,1);
        y = U(j).u(:,1);
        x = conv(full(x),bf(:,p));
        y = conv(full(y),bf(:,q));
        x = x(d);
        y = y(d);
        X = [X x.*y];

        %-Indices and regressor names
        %------------------------------------------------------------------
        Xname{end + 1} = sprintf('%s*bf(%i)x%s*bf(%i)',...
                            U(i).name{1}, p,...
                            U(j).name{1}, q);
        ind(end + 1)   = size(X,2);
        ip(end + 1)    = 1;
    end
    end
    Fc(end + 1).i = ind;
    Fc(end).name  = [U(i).name{1} 'x' U(j).name{1}];
    Fc(end).p     = ip;
end
end

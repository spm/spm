function [X,Xname,Fc] = spm_Volterra(U,bf,V)
% generalized convolution of inputs (U) with basis set (bf)
% FORMAT [X,Xname,Fc] = spm_Volterra(U,bf,V);
% U          -  input structure array
% bf         -  Basis functions
% V          -  [1 or 2] order of Volterra expansion [default = 1]
%
% X          -  Design Matrix
% Xname      -  names of regressors [columns] in X
% Fc(j).i    -  indices pertaining to input i (and interactions)
% Fc(j).name -  names pertaining to input i   (and interactions)
%___________________________________________________________________________
%
% For first order expansions spm_Volterra simply convolves the causes
% (e.g. stick functions) in U.u by the basis functions in bf to create
% a design matrix X.  For second order expansions new entries appear
% in ind, bf and name that correspond to the interaction among the
% original causes. The basis functions for these effects are two dimensional
% and are used to assemble the second order kernel in spm_graph.m.
% Second order effects are computed for only the first column of U.u.
%___________________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% Karl Friston
% $Id: spm_Volterra.m 2021 2008-08-27 10:05:32Z volkmar $



% 1st order terms
%---------------------------------------------------------------------------
if nargin == 2, V = 1; end

% Construct X
%===========================================================================

% 1st order terms
%---------------------------------------------------------------------------
X     = [];
Xname = {};
ind   = {};
Uname = {};
Fc    = {};
for i = 1:length(U)
    ind   = [];
    for k = 1:size(U(i).u,2)
    for p = 1:size(bf,2)
        x      = U(i).u(:,k);
        d      = 1:length(x);
        x      = conv(full(x),bf(:,p));
        x      = x(d);
        X      = [X x];

        % indices and regressor names
        %-----------------------------------------------------------
        str            = sprintf('%s*bf(%i)',U(i).name{k},p);
        Xname{end + 1} = str;
        ind(end + 1)   = size(X,2);
    end
    end
    Fc(end + 1).i = ind;
    Fc(end).name  = U(i).name{1};
end

% return if first order
%---------------------------------------------------------------------------
if V == 1, return, end

% 2nd order terms
%---------------------------------------------------------------------------
for i = 1:length(U) 
for j = i:length(U)
    ind   = [];
    for p = 1:size(bf,2)
    for q = 1:size(bf,2)
        x      = U(i).u(:,1);
        y      = U(j).u(:,1);
        x      = conv(full(x),bf(:,p));
        y      = conv(full(y),bf(:,q));
        x      = x(d);
        y      = y(d);
        X      = [X x.*y];

        % indices and regressor names
        %-----------------------------------------------------------    
        str            = sprintf('%s*bf(%i)x%s*bf(%i)',...
                          U(i).name{1},p,...
                          U(j).name{1},q);
        Xname{end + 1} = str;
        ind(end + 1)   = size(X,2);
    end
    end
    Fc(end + 1).i = ind;
    Fc(end).name  = [U(i).name{1} 'x' U(j).name{1}];
end
end

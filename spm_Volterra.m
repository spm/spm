function [X,Xname,ind,Uname] = spm_Volterra(U,bf,V)
% generalized convolution of inputs (U) with basis set (bf)
% FORMAT [X,Xname,ind,Uname] = spm_Volterra(U,bf,V);
% U{i}     -  input structure
% bf       -  Basis functions
% V        -  [1 or 2] order of Volterra expansion [default = 1]
%
% X        -  Design Matrix
% Xname    -  names of regressors [columns] in X
% ind{i}   -  indices pertaining to input i (and interactions)
% Canme{i} -  names pertaining to input i (and interactions)
%___________________________________________________________________________
%
% For first order expansions spm_Volterra simply convolves the causes
% (e.g. stick functions) in U.u by the basis functions in bf to create
% a design matrix X.  For second order expansions new entries appear
% in ind, bf and name that correspond to the interaction among the
% orginal causes. The basis functions for these efects are two dimensional
% and are used to assemble the second order kernel in spm_graph.m.
% Second order effects are computed for only the first column of U.u.
%___________________________________________________________________________
% %W% Karl Friston %E%


% 1st order terms
%---------------------------------------------------------------------------
if nargin == 2, V == 1; end

% Construct X
%===========================================================================

% 1st order terms
%---------------------------------------------------------------------------
X     = [];
Xname = {};
ind   = {};
Uname = {};
for i = 1:length(U)
    ind     = {ind{:} {}};
    for k = 1:size(U{i}.u,2)
	for p = 1:size(bf,2)
		x      = U{i}.u(:,k);
		d      = 1:length(x);
		x      = conv(full(x),bf(:,p));
		x      = x(d);
		X      = [X x];

		% indices and regressor names
		%-----------------------------------------------------------
		str            = sprintf('%s*bf(%i)',U{i}.Uname{k},p);
		Xname{end + 1} = str;
		ind{end}       = [ind{end} size(X,2)];
	end
    end
    Uname{end + 1} = U{i}.Uname{1};
end

% return if first order
%---------------------------------------------------------------------------
if V == 1, return, end

% 2nd order terms
%---------------------------------------------------------------------------
for i = 1:length(U) 
for j = i:length(U)
    	ind   = {ind{:} {}};
	for p = 1:size(bf,2)
	for q = 1:size(bf,2)
		x      = U{i}.u(:,1);
		y      = U{j}.u(:,1);
		x      = conv(full(x),bf(:,p));
		y      = conv(full(y),bf(:,q));
		x      = x(d);
		y      = y(d);
		X      = [X x.*y];

		% indices and regressor names
		%-----------------------------------------------------------	
		str            = sprintf('%s*bf(%i)x%s*bf(%i)',...
					      U{i}.Uname{1},p,...
					      U{j}.Uname{1},q);
		Xname{end + 1} = str;
		ind{end}       = [ind{end} size(X,2)];
	end
	end
	Uname{end + 1} = [U{i}.Uname{1} 'x' U{j}.Uname{1}];
end
end

function [X,Xn,IND,BF,name] = spm_Volterra(SF,BF,name,N)
% returns [design] matrix of explanatory variables
% FORMAT [X Xn IND BF name] = spm_Volterra(SF,BF,name,N);
% SF{i}   -  multivariate causes: SF{i}(:,j) = casue i,  expansion j
% BF{i}   -  Basis functions:     BF{i}      = basis set for cause i
% name{i} -  name of cause i
% N       -  [1 or 2] order of Volterra expansion
%
% X       -  Design Matrix
% Xn{i}   -  name of cause i (now including interactions among causes)
% IND{i}  -  indices pertaining to cause i (and interactions)
%___________________________________________________________________________
%
% For first order expansions spm_Volterra simply convolves the causes
% (e.g. stick functions) in SF by the basis functions in BF to create
% a design matrix X.  For second order expansions new entries appear
% in IND, BF and name that correspond to the interaction among the
% orginal causes (if the events are sufficiently close in time).
% The basis functions for these are two dimensional and are used to
% assemble the second order kernel in spm_graph.m.  Second order effects
% are computed for only the first column of SF.
%___________________________________________________________________________
% %W% Karl Friston %E%


% Construct X
%===========================================================================

% 1st order terms
%---------------------------------------------------------------------------
X     = [];
Xn    = {};
IND   = cell(1,size(SF,2));
for i = 1:size(SF,2)
	for j = 1:size(BF{i},2)
		for k = 1:size(SF{i},2)
			x      = SF{i}(:,k);
			d      = 1:length(x);
			x      = conv(full(x),BF{i}(:,j));
			x      = x(d);
			X      = [X x];
			IND{i} = [IND{i} size(X,2)];
			if size(SF{i},2) > 1
				str = [name{i} sprintf('(%i)[%i]',j,(k - 1))];
			else
				str = [name{i} sprintf('(%i)',j)];
			end
			Xn{end + 1} = str;
		end
	end
end

% return if first order
%---------------------------------------------------------------------------
if N == 1, return, end

% 2nd order terms
%---------------------------------------------------------------------------
k     = length(name);
for i = 1:size(SF,2) 
	for j = i:size(SF,2)

		% ensure events can interact
		%-----------------------------------------------------------
		skip = 0;
		if i == j
			p    = diff(find(SF{i}(:,1)));
			skip = (size(BF{i},1) <= min(p)) | ~any(diff(p));
		end

		if ~skip

			k       = k + 1;
			ind     = [];
			bf      = {};
			for p = 1:size(BF{i},2)
				for q = 1:size(BF{j},2)
					ni     = [name{i} sprintf('(%i)',p)];
					nj     = [name{j} sprintf('(%i)',q)];
					x      = SF{i}(:,1);
					y      = SF{j}(:,1);
					x      = conv(full(x),BF{i}(:,p));
					y      = conv(full(y),BF{j}(:,q));
					x      = x(d);
					y      = y(d);
					X      = [X x.*y];		
					ind    = [ind size(X,2)];
					Xn{end + 1} = [ni ' x ' nj];
					bf     = [bf {BF{i}(:,p)*BF{j}(:,q)'}];
				end
			end
			name{k} = [name{i} ' x ' name{j}];
			IND{k}  = ind;
			BF{k}   = bf;
		end
	end
end

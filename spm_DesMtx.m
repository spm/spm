function [X,Enames,Index]=spm_DesMtx(P1,P2,P3,P4,P5,P6,P7,P8,P9,P10,P11,P12);
% Design matrix construction from factor level and covariate vectors
% FORMAT [X,Enames]=spm_DesMtx(<FCLevels-constraint-name> list);
%
% <FCLevels-constraints-name>
%	- Set of arguments specifying a portion of design matrix.
%	- name parameter, or constraint and name parameters, can be omitted 
%	  from the triple required for maximum specification.
%_______________________________________________________________________
%
% Returns design matrix corresponding to given vectors containing
% levels of a factor; two way interactions; covariates (n vectors);
% ready-made sections of design matrix; and factor by covariate
% interactions.
%
% The specification for the design matrix is passed in sets of arguments,
% each set corresponding to a particular Factor/Covariate/&c., specifying
% a section of the design matrix. The set of arguments consists of the 
% FCLevels matrix (Factor/Covariate levels), an optional constraint string,
% and an optional (string) name matrix containing the names of the 
% Factor/Covariate/&c.
%
% MAIN EFFECTS: For a main effect, or single factor, the FCLevels
% matrix is an integer vector whose values represent the levels of the
% factor. The integer factor levels need not be positive, nor in order.
% Effects for the factor levels are entered into the design matrix *in
% increasing order* of the factor levels. Check Enames to find out which
% columns correspond to which levels of the factor.
%
% TWO WAY INTERACTIONS: For a two way interaction effect between two
% factors, the FCLevels matrix is an nx2 integer matrix whose columns
% indicate the levels of the two factors. An effect is included for each
% unique combination of the levels of the two factors.
%
% CONSTRAINTS: Each FactorLevels vector/matrix may be followed by an 
% (optional) ConstraintString. Constraints are applied to the last level
% of the factor.
%  ConstraintStrings for main effects are:
%                  '-'   - No Constraints
%                  '+0'  - sum-to-zero constraints
%                  '+0m' - Implicit sum-to-zero constraints
%                  '.'   - CornerPoint constraints
%  Constraints for two way interaction effects are
%    '-'                 - No Constraints
%    '+i0','+j0','+ij0'  - sum-to-zero constraints
%    '.i', '.j', '.ij'   - CornerPoint constraints
%    '+0im', '+0jm'      - Implicit sum-to-zero constraints
%
% The implicit sum-to-zero constraints "mean correct" appropriate rows
% of the relevant design matrix block. For a main effect, constraint
% '+0m' "mean corrects" the main effect block across columns,
% corresponding to factor effects B_i, where B_i = B'_i - mean(B'_i) :
% The B'_i are the fitted parameters, effectively *relative* factor
% parameters, relative to their mean. This leads to a rank deficient
% design matrix block. If pinv is used to solve the least squares
% problem, then the solution with smallest L2 norm is found, which has
% mean(B'_i)=0 provided the remainder of the design is unique (design
% matrix blocks of full rank).  In this case therefore the B_i are
% identically the B'_i - the mean correction imposes the constraint.
%      
% COVARIATES: The FCLevels matrix here is an nxc matrix whose columns
% contain the covariate values. An effect is included for each covariate.
% Covariates are identified by ConstraintString 'C'.
%
% PRE-SPECIFIED DESIGN BLOCKS: ConstraintString 'X' identifies a
% ready-made bit of design matrix - the effect is the same as 'C'.
%
% FACTOR BY COVARIATE INTERACTIONS: are identified by ConstraintString
% 'FxC'. The FCLevels matrix must contain two columns, one of integers
% (indicating the factor levels) and one containing the covariate. (If
% the covariate has only integer values, it *must* be in the second
% column!)
%
% Each Factor/Covariate can be 'named', by passing a name string.
% Pass a string matrix, with rows naming the factors/covariates in the
% respective columns of the FCLevels matrix.  These names default to
% <Fac> &c., and are used in the construction of the Enames effect names
% matrix.
% E.g. for an interaction, spm_DesMtx([F1,F2],'+ij0',['Fac1';'Fac2'])
%
% Enames returns a string matrix whose successive rows describe the
% effects parameterised in the corresponding columns of the design
% matrix. The characters `_' and `_' are used in the names, for example
% `Fac1_2-Fac2_3' would efer to the effect for the interaction of the
% two factors Fac1 & Fac2, at the 2nd level of the former and the 3rd
% level of the latter. The `_' & `-' characters are used by other
% programs (spm_DesMtxSca), and should therefore be avoided when
% naming effects and ovariates.
%
% The function uses recursion to break up the problem.
%
% Integer Index matrix is returned if only a single block of design
% matrix is being computed (single set of parameters). It indexes the
% actual order of the effects in the design matrix block.  (Factor levels
% are introduced in order, regardless of order of appearence in the
% factor index matrices, so that the parameters vector has a sensible
% order.) This is used to aid recursion.
%
%__________________________________________________________________________
% %W% Andrew Holmes %E%


%-Parse arguments for recursive construction of design matrices
%=======================================================================
if nargin<3 P3=[]; end
if nargin<2 P2=[]; end
if nargin==0 error('Insufficient arguments'), end

if ~isstr(P2)
	if nargin>=2
		Args='P2'; for t=3:nargin Args=[Args, ',P' num2str(t)]; end
		[X1,Enames1]=spm_DesMtx(P1);
		[X2,Enames2]=eval(['spm_DesMtx(',Args ')']);
		X=[X1,X2]; Enames=str2mat(Enames1,Enames2);
		return
	end
elseif ~isstr(P3)
	if nargin>=3
		Args='P3'; for t=4:nargin Args=[Args, ',P' num2str(t)]; end
		[X1,Enames1]=spm_DesMtx(P1,P2);
		[X2,Enames2]=eval(['spm_DesMtx(',Args ')']);
		X=[X1,X2]; Enames=str2mat(Enames1,Enames2);
		return
	end
elseif nargin>=4
		Args='P4'; for t=5:nargin Args=[Args, ',P' num2str(t)]; end
		[X1,Enames1]=spm_DesMtx(P1,P2,P3);
		[X2,Enames2]=eval(['spm_DesMtx(',Args ')']);
		X=[X1,X2]; Enames=str2mat(Enames1,Enames2);
		return
end


%-Computation
%=======================================================================

%-Sort out arguments P1, P2 & P3
%-----------------------------------------------------------------------
I=P1;
%-If I is a vector, make it a row vector
if size(I,1)==1 I=I'; end
%-Sort out constraint and Factor/Covariate name parameters from P2 & P3
if isempty(P2), Constraint='-'; else, Constraint=P2; end
FCnames=P3;

%-Empty I case
%-----------------------------------------------------------------------
if isempty(I)
	X=[];
	Enames=[];
	Index=[];

%-Factor dependent covariate effect
%-----------------------------------------------------------------------
elseif strcmp(Constraint,'FxC')
	if isempty(FCnames) FCnames = [' ';' ']; end
	if size(FCnames,1)~=2
		error('Two row FCnames required for FxC interaction'), end

	if size(I,2)~=2 error('Two column I for ''FxC'''), I, end
	%-Factor by Covariate interaction
	if all(floor(I(:,1))==ceil(I(:,1)))
		C=I(:,2); I=I(:,1);
		Fname=FCnames(1,:); Cname=FCnames(2,:);
	elseif all(floor(I(:,2))==ceil(I(:,2)))
		C=I(:,1); I=I(:,2);
		Fname=FCnames(2,:); Cname=FCnames(1,:);
	else
		error('Neither column of a ''FxC'' contains integer indicies')
	end % if (all integer in one column)
	if strcmp(Fname,' ') Fname='Fac'; Cname='Cov'; end

	[X,Enames,Index]=spm_DesMtx(I,'-',deblank(Fname));
	X=X.*meshgrid(C,1:size(X,2))';
	Enames=[Enames, setstr(ones(length(Index),1)*['-',deblank(Cname)])];

%-Covariate effect, or ready-made design matrix
%-----------------------------------------------------------------------
elseif strcmp(Constraint,'C')|strcmp(Constraint,'X')
	%-I contains a covariate (C), or is to be inserted "as is" (X)
	X=I;

	if isempty(FCnames)
		if strcmp(Constraint,'C') FCnames='<Cov>'; end
		if strcmp(Constraint,'X') FCnames='<X>'; end
	end % if
	if (size(FCnames,1)>1) & (size(FCnames,1)~=size(X,2))
		error('FCnames doesn''t match covariate/X matrix'), end
	if size(FCnames,1)==1
		Enames=setstr(ones(size(X,2),1)*FCnames);
	else
		Enames=FCnames;
	end

%-Simple main effect
%-----------------------------------------------------------------------
elseif size(I,2)==1
	%-Sort out unique factor levels
	if ~all(floor(I)==ceil(I)) error('Non-integer indicator vector'), end
	temp=sort(I);
	Index=temp(logical([1;diff(temp(:))>0]));
	clear temp
	
	%-Determine sizes of design matrix X
	nXrows=size(I,1);
	nXcols=length(Index);
	
	%-Set up unconstrained X matrix
	%-Columns in ascending order of corresponding factor level
	X=zeros(nXrows,nXcols);
	for p_i=1:nXcols		%-p_i is position of i in Index
		X(:,p_i)=I==Index(p_i);
		%-Can't use for i=Index X(:,i)=I==i in case
		% Index has holes &/or doesn't start at 1!
	end % for i
	
	%-Impose constraint if specified & if more than one effect
	%-Apply uniqueness constraints ('+0' & '.')  to last effect, which is
	% in last column, since column i corresponds to level Index(i)
	if strcmp(Constraint,'+0') & nXcols>1
		X(X(:,nXcols),:)=-1*ones(sum(X(:,nXcols)),nXcols);
		X(:,nXcols)=[];
		Index(nXcols)=[];
	elseif strcmp(Constraint,'+0m') & nXcols>1
		X = X - 1/nXcols;
	elseif strcmp(Constraint,'.') & nXcols>1
		X(:,nXcols)=[];
		Index(nXcols)=[];
	elseif strcmp(Constraint,'-')
	elseif ~strcmp(Constraint,'-') & nXcols==1
		error('Can''t constrain constant effect')
	else
		error('unknown constraint type for main effect')
	end % if SumToZero

	%-Construct effect name index
	if isempty(FCnames) FCnames = '<Fac>'; end
	if size(FCnames,1)>1 error('Too many FCnames in matrix'), end
	if ~isempty(Index)
		Enames=[FCnames,'_',int2str(Index(1))];
		for p_i=2:length(Index)
			Enames = str2mat(Enames,...
				[FCnames,'_',int2str(Index(p_i))]);
		end % (for)
	end % (if)


%-Two way interaction effects
%-----------------------------------------------------------------------
elseif size(I,2)==2
	if ~all(floor(I(:))==ceil(I(:)))
		error('Non-integer indicator vector')
	end
	IJ=I; I=IJ(:,1); J=IJ(:,2);

	%-Create design matrix for all interaction effects
	%-Make "raw" index to IJ pairs
	rIndex=(max(J)-min(J)+1)*(I-min(I))+(J-min(J)+1);
	[X,null,sIndex]=spm_DesMtx(rIndex);	%-make design matrix
	%-Build Index matrix
	Index=[];
	for c=1:length(sIndex)
		rp=min(find(rIndex==sIndex(c)));
		Index=[Index, [I(rp);J(rp)]];
	end

	%-Impose Constraints
	if (Constraint(length(Constraint))=='m')
		%-Implicit sum to zero constraints
		SumIToZero=strcmp(Constraint,'+i0m');
		SumJToZero=strcmp(Constraint,'+j0m');
		if (~(SumIToZero|SumJToZero))
		    error('Illegal SumToZero constraint type for interaction')
		end % if (error)
		if (SumIToZero)
			%-impose SumIToZero constraints
			Js = sort(Index(2,:)); Js = Js([1,diff(Js)>0]);
			for j = Js
				rows = find(J==j);
				cols = find(Index(2,:)==j);
				if length(cols)==1
				   error('Only one level: Can''t constrain')
				end
				X(rows,cols) = X(rows,cols) - 1/length(cols);
			end
		end % if SumIToZero

		if (SumJToZero)
			Is = sort(Index(1,:)); Is = Is([1,diff(Is)>0]);
			for i = Is
				rows = find(I==i);
				cols = find(Index(1,:)==i);
				if length(cols)==1
				   error('Only one level: Can''t constrain')
				end
				X(rows,cols) = X(rows,cols) - 1/length(cols);
			end
		end % if SumJToZero

	elseif (Constraint(1)=='+')
		%-Explicit sum to zero constraints
		SumIToZero=strcmp(Constraint,'+i0')|strcmp(Constraint,'+ij0');
		SumJToZero=strcmp(Constraint,'+j0')|strcmp(Constraint,'+ij0');
		if (~(SumIToZero|SumJToZero))
		    error('Illegal SumToZero constraint type for interaction')
		end % if (error)
		if (SumIToZero)
			%-impose SumIToZero constraints
			i=max(Index(1,:));
			if i==min(Index(1,:))
				error('Only one i level: Can''t constrain')
			end
			cols=find(Index(1,:)==i); % columns to delete
			for c=cols
				j=Index(2,c);
				t_cols=find(Index(2,:)==j);
				t_rows=find(X(:,c));
				%-This ij equals -sum(ij) over other i
				% (j fixed for this col c).
				%-So subtract weight of this ij factor frow
				% weights for all other ij factors for this j
				% to impose the constraint.
				X(t_rows,t_cols) = X(t_rows,t_cols)...
					-1*meshgrid(...
					X(t_rows,c),1:length(t_cols))';
%-(This next line would do it, but only first time round, when all weights
% (are 1, and only one weight per row for this j.
% X(t_rows,t_cols)=-1*ones(length(t_rows),length(t_cols));
			end % for c
			%-delete columns
			X(:,cols)=[]; Index(:,cols)=[];
		end % if SumIToZero

		if (SumJToZero)
			%-impose SumJToZero constraints
			j=max(Index(2,:));
			if j==min(Index(2,:))
				error('Only one j level: Can''t constrain')
			end
			cols=find(Index(2,:)==j); % columns to delete
			for c=cols
				i=Index(1,c);
				t_cols=find(Index(1,:)==i);
				t_rows=find(X(:,c));
				X(t_rows,t_cols) = X(t_rows,t_cols)...
					-1*meshgrid(...
					X(t_rows,c),1:length(t_cols))';
			end % for c
			%-delete columns
			X(:,cols)=[]; Index(:,cols)=[];
		end % if SumJToZero

	elseif (Constraint(1)=='.')
		%-Corner point constraints
		CornerPointI=strcmp(Constraint,'.i')|strcmp(Constraint,'.ij');
		CornerPointJ=strcmp(Constraint,'.j')|strcmp(Constraint,'.ij');
		if (~(CornerPointI|CornerPointJ))
		error('Illegal CornerPoint constraint type for interaction')
		end % if (error)
		if (CornerPointI)
			%-impose CornerPointI constraints
			i=max(Index(1,:));
			if i==min(Index(1,:))
				error('Only one i level: Can''t constrain')
			end
			cols=find(Index(1,:)==i); % columns to delete
			%-delete columns
			X(:,cols)=[]; Index(:,cols)=[];
		end % if CornerPointI

		if (CornerPointJ)
			%-impose CornerPointJ constraints
			j=max(Index(2,:));
			if j==min(Index(2,:))
				error('Only one j level: Can''t constrain')
			end
			cols=find(Index(2,:)==j); % columns to delete
			%-delete columns
			X(:,cols)=[]; Index(:,cols)=[];
		end % if CornerPointJ

	elseif strcmp(Constraint,'-')
	else
		error('unknown constraint type for interaction effects')
	end % if SumToZero type constraints
	
	%-Construct effect name index
	if isempty(FCnames) FCnames = ['<Fac1>';'<Fac2>']; end
	if size(FCnames,1)~=2
		error('Two row FCnames required for interaction'), end
	Iname=deblank(FCnames(1,:)); Jname=deblank(FCnames(2,:));
	Enames=[Iname,' x ',Jname,'_',...
		int2str(Index(1,1)),'x',int2str(Index(2,1))];
	for Xcol=2:size(Index,2)
		Enames=str2mat(Enames,...
			[Iname,' x ',Jname,'_',...
			int2str(Index(1,Xcol)),'x',int2str(Index(2,Xcol))]);
	end % for


%-Mis-specified arguments - ERROR
%-----------------------------------------------------------------------
else
	error('Can only do main effects and two way interactions!')
end % if size (Switch for main effect/interactions matrix)

function [X,Enames,Index]=spm_DesMtx(varargin);
% Design matrix construction from factor level and covariate vectors
% FORMAT [X,Enames] = spm_DesMtx(<FCLevels-Constraint-FCnames> list)
% FORMAT [X,Enames,Index] = spm_DesMtx(FCLevels,Constraint,FCnames)
%
% <FCLevels-Constraints-FCnames>
%        - Set of arguments specifying a portion of design matrix (see below)
%        - FCnames parameter, or Constraint and FCnames parameters, are optional
%        - A list of multiple <FCLevels-Constraint-FCnames> triples can be
% 	   specified, where FCnames or Constraint-FCnames may be omitted
% 	   within any triple. The program then works recursively.
%
% X      - Design matrix
% Enames - Effect names (constructed from FCnames)
% Index  - Integer index of factor levels
%
%                           ----------------
%
% FORMAT [nX,nEnames] = spm_DesMtx('sca',X1,Enames1,X2,Enames2,...)
% Produces a scaled design matrix nX with max(abs(nX(:))<=1, suitable
% for imaging with: image((nX+1)*32)
% X1,X2,...             - Design matrix partitions
% Enames1, Enames2,...  - Corresponding effect name string matrices (optional)
% nX                    - Scaled design matrix
% nEnames               - COncatenated effect names for columns of nX
%
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
%    '+i0m', '+j0m'      - Implicit sum-to-zero constraints
%
% The implicit sum-to-zero constraints "mean correct" appropriate rows
% of the relevant design matrix block. For a main effect, constraint
% '+0m' "mean corrects" the main effect block across columns,
% corresponding to factor effects B_i, where B_i = B'_i - mean(B'_i) :
% The B'_i are the fitted parameters, effectively *relative* factor
% parameters, relative to their mean. This leads to a rank deficient
% design matrix block. If Matlab's pinv, which uses a Moore-Penrose
% pseudoinverse, is used to solve the least squares problem, then the
% solution with smallest L2 norm is found, which has mean(B'_i)=0
% provided the remainder of the design is unique (design matrix blocks
% of full rank). In this case therefore the B_i are identically the
% B'_i - the mean correction imposes the constraint.
%      
% COVARIATES: The FCLevels matrix here is an nxc matrix whose columns
% contain the covariate values. An effect is included for each covariate.
% Covariates are identified by ConstraintString 'C'.
%
% PRE-SPECIFIED DESIGN BLOCKS: ConstraintString 'X' identifies a
% ready-made bit of design matrix - the effect is the same as 'C'.
%
% FACTOR BY COVARIATE INTERACTIONS: are identified by ConstraintString
% 'FxC'. The last column is understood to contain the covariate. Other
% columns are taken to contain integer FactorLevels vectors. The
% (unconstrained) interaction of the factors is interacted with the
% covariate.
%
% NAMES: Each Factor/Covariate can be 'named', by passing a name
% string.  Pass a string matrix, with rows naming the
% factors/covariates in the respective columns of the FCLevels matrix.
% These names default to <Fac>, <Cov>, <Fac1>, <Fac2> &c., and are used
% in the construction of the Enames effect names matrix.
% E.g. for an interaction, spm_DesMtx([F1,F2],'+ij0',['subj';'cond'])
% giving effect names such as subj*cond_{1,2} etc...
%
% Enames returns a string matrix whose successive rows describe the
% effects parameterised in the corresponding columns of the design
% matrix. `Fac1*Fac2_{2,3}' would refer to the effect for the
% interaction of the two factors Fac1 & Fac2, at the 2nd level of the
% former and the 3rd level of the latter. Other forms are
%  - Simple main effect (level 1)        : <Fac>_{1}
%  - Three way interaction (level 1,2,3) : <Fac1>*<Fac2>*<Fac3>_{1,2,3}
%  - Two way factor interaction by covariate interaction :
%                                        : <Fac1>*<Fac2>_{1,1}*<Cov>
%  - Column 3 of prespecified DesMtx block (if unnamed)
%                                        : <X> (1)
% The special characters `_*{}()' are recognised by the scaling
% function (spm_DesMtxSca('Xsca',...), and should therefore be avoided
% when naming effects and covariates.
%
% INDEX: An Integer Index matrix is returned if only a single block of
% design matrix is being computed (single set of parameters). It
% indexes the actual order of the effects in the design matrix block.
% (Factor levels are introduced in order, regardless of order of
% appearence in the factor index matrices, so that the parameters
% vector has a sensible order.) This is used to aid recursion.
%
%                           ----------------
%
% The design matrix scaling feature is designed to return a scaled
% version of a design matrix, with values in [-1,1], suitable for
% visualisation. Special care is taken to apply the same normalisation
% to blocks of design matrix reflecting a single effect, to preserve
% appropriate relationships between columns. Identification of effects
% corresponding to columns of design matrix portions is via the effect
% names matrices. The design matrix may be passed in any number of
% parts, provided the corresponding effect names are given. It is
% assummed that the block representing an effect is contained within a
% single partition. Partitions supplied without corresponding effect
% names are scaled on a column by column basis, the effects labelled as
% <UnSpec> in the returned nEnames matrix.
% 
% Effects are identified using the special characters `_*{}()' used in
% effect naming as follows: (here ? is a wildcard)
%       - ?_{?}         - main effect or interaction of main effects
%       - ?_{?}*?       - factor by covariate interaction
%       - ?(?)          - part of a pre-specified block
% Blocks are identified by looking for runs of effects of the same type
% with the same names: E.g. a block of main effects for factor 'Fac1'
% would have names like Fac1_{?}.
% 
% Scaling is as follows:
% 	* No scaling is carried out if max(abs(tX(:))) is in [.4,1]
% 	  This protects dummy variables from normalisation, even if
% 	  using implicit sum-to-zero constraints.  * If the block has
% 	a single value, it's replaced by 1's * FxC blocks are
% 	normalised so the covariate values cover [-1,1]
% 	  but leaving zeros as zero.  * Otherwise, block is scaled to
% 	cover [-1,1].
%
%_______________________________________________________________________
% %W% Andrew Holmes %E%


%-Parse arguments for recursive construction of design matrices
%=======================================================================
if nargin==0 error('Insufficient arguments'), end

if ischar(varargin{1})
	%-Non-recursive action string usage
	Constraint=varargin{1};
elseif nargin>=2 & ~ischar(varargin{2})
	[X1,Enames1]=spm_DesMtx(varargin{1});
	[X2,Enames2]=spm_DesMtx(varargin{2:end});
	X=[X1,X2]; Enames=strvcat(Enames1,Enames2);
	return
elseif nargin>=3 & ~ischar(varargin{3})
	[X1,Enames1]=spm_DesMtx(varargin{1:2});
	[X2,Enames2]=spm_DesMtx(varargin{3:end});
	X=[X1,X2]; Enames=strvcat(Enames1,Enames2);
	return
elseif nargin>=4
	[X1,Enames1]=spm_DesMtx(varargin{1:3});
	[X2,Enames2]=spm_DesMtx(varargin{4:end});
	X=[X1,X2]; Enames=strvcat(Enames1,Enames2);
	return
else
	%-If I is a vector, make it a column vector
	I=varargin{1}; if size(I,1)==1, I=I'; end
	%-Sort out constraint and Factor/Covariate name parameters
	if nargin<2, Constraint='-'; else, Constraint=varargin{2}; end
	if isempty(I), Constraint='mt'; end
	if nargin<3, FCnames=''; else, FCnames=varargin{3}; end
end



switch Constraint, case 'mt'                              %-Empty I case
%=======================================================================
X=[];
Enames=[];
Index=[];



case {'C','X'}           %-Covariate effect, or ready-made design matrix
%=======================================================================
%-I contains a covariate (C), or is to be inserted "as is" (X)
X = I;

%-Construct effect name index
%-----------------------------------------------------------------------
if isempty(FCnames)
	if strcmp(Constraint,'C'), FCnames='<Cov>';
		else, FCnames='<X>'; end
end

if size(FCnames,1)==1 & size(X,2)>1
	Enames = '';
	for i=1:size(X,2)
		Enames = strvcat(Enames,sprintf('%s (%d)',FCnames,i));
	end
elseif size(FCnames,1)~=size(X,2)
	error('FCnames doesn''t match covariate/X matrix')
else
	Enames = FCnames;
end



case '-(1)'                                         %-Simple main effect
%=======================================================================
%-Sort out arguments & unique factor levels
%-----------------------------------------------------------------------
if size(I,2)>1, error('Simple main effect requires vector index'), end
if any(I~=floor(I)), error('Non-integer indicator vector'), end
if isempty(FCnames), FCnames = '<Fac>';
elseif size(FCnames,1)>1, error('Too many FCnames in matrix'), end

tmp   = sort(I');
Index = tmp([1,1+find(diff(tmp))]);
clear tmp

%-Determine sizes of design matrix X
nXrows = size(I,1);
nXcols = length(Index);

%-Set up unconstrained X matrix & construct effect name index
%-----------------------------------------------------------------------
%-Columns in ascending order of corresponding factor level
X      = zeros(nXrows,nXcols);
Enames = '';
for ii=1:nXcols			%-ii indexes i in Index
	X(:,ii) = I==Index(ii);
	%-Can't use: for i=Index, X(:,i) = I==i; end
	% in case Index has holes &/or doesn't start at 1!
	Enames = strvcat(Enames,sprintf('%s_{%d}',FCnames,Index(ii)));
end



case {'-','-(*)'}             %-Main effect or unconstrained interaction
%=======================================================================
if size(I,2)==1
	[X,Enames,Index] = spm_DesMtx(I,'-(1)',FCnames);
	return
end

if any((I(:))~=floor(I(:))), error('Non-integer indicator vector'), end

%-Create design matrix for all (interaction) effects
%-----------------------------------------------------------------------
%-Make "raw" index to unique effects
nI     = I - ones(size(I,1),1)*min(I);
tmp    = max(I)-min(I)+1;
tmp    = [fliplr(cumprod(tmp(end:-1:2))),1];
rIndex = sum(nI.*(ones(size(I,1),1)*tmp),2)+1;

[X,null,sIndex]=spm_DesMtx(rIndex,'-(1)');

%-Build Index matrix
%-----------------------------------------------------------------------
Index = zeros(size(I,2),length(sIndex));
for c = 1:length(sIndex)
	Index(:,c) = I(min(find(rIndex==sIndex(c))),:)';
end

%-Construct effect name index
%-----------------------------------------------------------------------
if isempty(FCnames)
	tmp = ['<Fac1>',sprintf('*<Fac%d>',2:size(I,2))];
elseif size(FCnames,1)==size(I,2)
	tmp = [deblank(FCnames(1,:))];
	for i = 2:size(I,2), tmp=[tmp,'*',deblank(FCnames(i,:))]; end
else
	error('#FCnames mismatches #Factors in interaction')
end

Enames = '';
for c = 1:size(Index,2)
    Enames = strvcat(Enames,...
	[sprintf('%s_{%d',tmp,Index(1,c)),...
	 sprintf(',%d',Index(2:end,c)),'}']);
end



case 'FxC'      %-Factor dependent covariate effect - Cov in last column
%=======================================================================
%-Check
%-----------------------------------------------------------------------
if size(I,2)==1, error('FxC requires multi-column I'), end

F = I(:,1:end-1);
C = I(:,end);

if ~all(all(F==floor(F),1),2)
	error('non-integer indicies in F partition of FxC'), end

if isempty(FCnames)
	Fnames = '';
	Cnames = '<Cov>';
elseif size(FCnames,1)==size(I,2)
	Fnames = FCnames(1:end-1,:);
	Cnames = deblank(FCnames(end,:));
else
	error('#FCnames mismatches #Factors+#Cov in FxC')
end

%-Set up design matrix X & names matrix
%-----------------------------------------------------------------------
[X,Enames,Index] = spm_DesMtx(F,'-',Fnames);
X = X.*(C*ones(1,size(X,2)));
Enames = [Enames,repmat(['*',Cnames],length(Index),1)];



case {'.','+0','+0m'}                   %-Constrained simple main effect
%=======================================================================

if size(I,2)~=1, error('Simple main effect requires vector index'), end

[X,Enames,Index] = spm_DesMtx(I,'-(1)',FCnames);

%-Impose constraint if more than one effect
%-----------------------------------------------------------------------
%-Apply uniqueness constraints ('.' & '+0')  to last effect, which is
% in last column, since column i corresponds to level Index(i)
nXcols = size(X,2);
if nXcols==1
	error('Can''t constrain a constant effect!')
elseif strcmp(Constraint,'.')
	X(:,nXcols)=[];	Enames(nXcols,:)=''; Index(nXcols)=[];
elseif strcmp(Constraint,'+0')
	X(find(X(:,nXcols)),:)=-1;
	X(:,nXcols)=[];	Enames(nXcols,:)=''; Index(nXcols)=[];
elseif strcmp(Constraint,'+0m')
	X = X - 1/nXcols;
end



case {'.i','.j','.ij','+i0','+j0','+ij0','+i0m','+j0m'}
                                           %-Two way interaction effects
%=======================================================================
if size(I,2)~=2, error('Two way interaction requires Nx2 index'), end

[X,Enames,Index] = spm_DesMtx(I,'-',FCnames);

%-Implicit sum to zero
%-----------------------------------------------------------------------
if any(strcmp(Constraint,{'+i0m','+j0m'}))
	SumIToZero = strcmp(Constraint,'+i0m');
	SumJToZero = strcmp(Constraint,'+j0m');

	if SumIToZero	%-impose implicit SumIToZero constraints
		Js = sort(Index(2,:)); Js = Js([1,1+find(diff(Js))]);
		for j = Js
			rows = find(I(:,2)==j);
			cols = find(Index(2,:)==j);
			if length(cols)==1
			   error('Only one level: Can''t constrain')
			end
			X(rows,cols) = X(rows,cols) - 1/length(cols);
		end
	end

	if SumJToZero	%-impose implicit SumJToZero constraints
		Is = sort(Index(1,:)); Is = Is([1,1+find(diff(Is))]);
		for i = Is
			rows = find(I(:,1)==i);
			cols = find(Index(1,:)==i);
			if length(cols)==1
			   error('Only one level: Can''t constrain')
			end
			X(rows,cols) = X(rows,cols) - 1/length(cols);
		end
	end

%-Explicit sum to zero
%-----------------------------------------------------------------------
elseif any(strcmp(Constraint,{'+i0','+j0','+ij0'}))
	SumIToZero = any(strcmp(Constraint,{'+i0','+ij0'}));
	SumJToZero = any(strcmp(Constraint,{'+j0','+ij0'}));

	if SumIToZero	%-impose explicit SumIToZero constraints
		i = max(Index(1,:));
		if i==min(Index(1,:))
			error('Only one i level: Can''t constrain'), end
		cols = find(Index(1,:)==i); %-columns to delete
		for c=cols
			j=Index(2,c);
			t_cols=find(Index(2,:)==j);
			t_rows=find(X(:,c));
			%-This ij equals -sum(ij) over other i
			% (j fixed for this col c).
			%-So subtract weight of this ij factor from
			% weights for all other ij factors for this j
			% to impose the constraint.
			X(t_rows,t_cols) = X(t_rows,t_cols)...
			    -X(t_rows,c)*ones(1,length(t_cols));
%-( Next line would do it, but only first time round, when all          )
% ( weights are 1, and only one weight per row for this j.              )
% X(t_rows,t_cols)=-1*ones(length(t_rows),length(t_cols));
		end
		%-delete columns
		X(:,cols)=[]; Enames(cols,:)=''; Index(:,cols)=[];
	end

	if SumJToZero	%-impose explicit SumJToZero constraints
		j = max(Index(2,:));
		if j==min(Index(2,:))
			error('Only one j level: Can''t constrain'), end
		cols=find(Index(2,:)==j);
		for c=cols
			i=Index(1,c);
			t_cols=find(Index(1,:)==i);
			t_rows=find(X(:,c));
			X(t_rows,t_cols) = X(t_rows,t_cols)...
			    -X(t_rows,c)*ones(1,length(t_cols));
		end
		%-delete columns
		X(:,cols)=[]; Enames(cols,:)=''; Index(:,cols)=[];
	end

%-Corner point constraints
%-----------------------------------------------------------------------
elseif any(strcmp(Constraint,{'.i','.j','.ij'}))
	CornerPointI = any(strcmp(Constraint,{'.i','.ij'}));
	CornerPointJ = any(strcmp(Constraint,{'.j','.ij'}));

	if CornerPointI	%-impose CornerPointI constraints
		i=max(Index(1,:));
		if i==min(Index(1,:))
			error('Only one i level: Can''t constrain')
		end
		cols=find(Index(1,:)==i); %-columns to delete
		%-delete columns
		X(:,cols)=[]; Enames(cols,:)=''; Index(:,cols)=[];
	end

	if CornerPointJ	%-impose CornerPointJ constraints
		j=max(Index(2,:));
		if j==min(Index(2,:))
			error('Only one j level: Can''t constrain')
		end
		cols=find(Index(2,:)==j);
		X(:,cols)=[]; Enames(cols,:)=''; Index(:,cols)=[];
	end
end



case 'sca'                                   %-Scale DesMtx for imaging
%=======================================================================
nX   = []; nEnames = ''; Carg = 2;

%-Loop through the arguments accumulating scaled design matrix nX
%-----------------------------------------------------------------------
while(Carg <= nargin)
    rX = varargin{Carg}; Carg=Carg+1;
    if Carg<=nargin & isstr(varargin{Carg})
	rEnames = varargin{Carg}; Carg=Carg+1;
    else	%-No names to work out blocks from - normalise by column
	rEnames = repmat('<UnSpec>',size(rX,2),1);
    end
    %-Pad out rEnames with 20 spaces to permit looking past line ends
    rEnames = [rEnames,repmat(' ',size(rEnames,1),20)];


    while(~isempty(rX))
	if size(rX,2)>1 & max(1,findstr(rEnames(1,:),'_{')) < ...
					max(0,find(rEnames(1,:)=='}'))
	%-Factor, interaction of factors, or FxC: find the rest
	%===============================================================
		c1 = max(findstr(rEnames(1,:),'_{'));
		d  = any(diff(abs(rEnames(:,1:c1+1))),2)...
			| ~any(rEnames(2:end,c1+2:end)=='}',2);
		t  = min(find([d;1]));
		if t>1
			%-Guard against following terms involving same factors
			[null,c2] = find(rEnames(1:t,c1+2:end)=='}');
			c2 = min(c2,[],2);
			d  = any(diff(abs(rEnames(1:t,c1+1+c2:c1+c2+21))),2);
			t  = min(find([d;1]));
		end

		%-Normalise block
		%-------------------------------------------------------
		tX = rX(:,1:t);
		if any(findstr(rEnames(1,c1:end),'}*'))
			%-Factor by covariate interaction
			C         = tX(tX~=0);
			tX(tX~=0) = 2*(C-min(C))/max(C-min(C))-1;
			nX        = [nX,tX];
		else
			nX = [nX,tXsca(tX)];
		end
		nEnames  = strvcat(nEnames,rEnames(1:t,:));
		rX(:,1:t) = []; rEnames(1:t,:)=[];


	elseif size(rX,2)>1 & max(1,find(rEnames(1,:)=='(')) < ...
					max(0,find(rEnames(1,:)==')'))
	%-Block: find the rest & normalise together
	%===============================================================
		c1 = max(find(rEnames(1,:)=='('));
		d  = any(diff(abs(rEnames(:,1:c1))),2)...
			| ~any(rEnames(2:end,c1+1:end)==')',2);
		t  = min(find([d;1]));

		%-Normalise block
		%-------------------------------------------------------
		nX = [nX,tXsca(rX(:,1:t))];
		nEnames  = strvcat(nEnames,rEnames(1:t,:));
		rX(:,1:t) = []; rEnames(1:t,:)=[];


	else                              %-Dunno! Just column normalise
	%===============================================================
		nX = [nX,tXsca(rX(:,1))];
		nEnames  = strvcat(nEnames,rEnames(1,:));
		rX(:,1) = []; rEnames(1,:)=[];

	end
    end
end

%-Trim off redundent trailing spaces from nEnames
X      = nX;
Enames = nEnames(:,1:...
		min(find([fliplr(cumprod(fliplr(all(nEnames==' ',1)))),1]))-1 );



otherwise                              %-Mis-specified arguments - ERROR
%=======================================================================
if ischar(varargin{1})
	error('unrecognised action string')
else
	error('unrecognised constraint type')
end

%=======================================================================
end



%=======================================================================
% - S U B F U N C T I O N S
%=======================================================================
function nX = tXsca(tX)
if nargin==0, nX=[]; return, end
if abs(max(abs(tX(:)))-0.7)<(.3+eps)
	nX = tX;
elseif all(tX(:)==tX(1))
	nX = ones(size(tX));
else
	nX = 2*(tX-min(tX(:)))/max(tX(:)-min(tX(:)))-1;
end

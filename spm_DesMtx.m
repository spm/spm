function [X,Pnames,Index,idx,jdx,kdx]=spm_DesMtx(varargin)
% Design matrix construction from factor level and covariate vectors
% FORMAT [X,Pnames] = spm_DesMtx(<FCLevels-Constraint-FCnames> list)
% FORMAT [X,Pnames,Index,idx,jdx,kdx] = spm_DesMtx(FCLevels,Constraint,FCnames)
%
% <FCLevels-Constraints-FCnames>
%        - set of arguments specifying a portion of design matrix (see below)
%        - FCnames parameter, or Constraint and FCnames parameters, are optional
%        - a list of multiple <FCLevels-Constraint-FCnames> triples can be
%      specified, where FCnames or Constraint-FCnames may be omitted
%      within any triple. The program then works recursively.
%
% X      - design matrix
% Pnames - paramater names as (constructed from FCnames) - a cellstr
% Index  - integer index of factor levels
%        - only returned when computing a single design matrix partition
%
% idx,jdx,kdx - reference vectors mapping I & Index (described below)
%             - only returned when computing a single design matrix partition
%               for unconstrained factor effects ('-' or '~')
%
%                           ----------------
% - Utilities:
%
% FORMAT i = spm_DesMtx('pds',v,m,n)
% Patterned data setting function - inspired by MINITAB's "SET" command
% v - base pattern vector
% m - (scalar natural number) #replications of elements of v [default 1]
% n - (scalar natural number) #repeats of pattern [default 1]
% i - resultant pattern vector, with v's elements replicated m times,
%     the resulting vector repeated n times.
%
% FORMAT [nX,nPnames] = spm_DesMtx('sca',X1,Pnames1,X2,Pnames2,...)
% Produces a scaled design matrix nX with max(abs(nX(:))<=1, suitable
% for imaging with: image((nX+1)*32)
% X1,X2,...             - Design matrix partitions
% Pnames1, Pnames2,...  - Corresponding parameter name string mtx/cellstr (opt)
% nX                    - Scaled design matrix
% nPnames               - Concatenated parameter names for columns of nX
%
% FORMAT Fnames = spm_DesMtx('Fnames',Pnames)
% Converts parameter names into suitable filenames
% Pnames - string mtx/cellstr containing parameter names
% Fnames - filenames derived from Pnames. (cellstr)
%
% FORMAT TPnames = spm_DesMtx('TeXnames',Pnames)
% Removes '*'s and '@'s from Pnames, so TPnames suitable for TeX interpretation
% Pnames  - string mtx/cellstr containing parameter names
% TPnames - TeX-ified parameter names
%
% FORMAT Map = spm_DesMtx('ParMap',aMap)
% Returns Nx2 cellstr mapping (greek TeX) parameters to English names,
% using the notation established in the SPMcourse notes.
% aMap - (optional) Mx2 cellstr of additional or over-ride mappings
% Map  - cellstr of parameter names (col1) and corresponding English names (col2)
%
% FORMAT EPnames = spm_DesMtx('ETeXnames',Pnames,aMap)
% Translates greek (TeX) parameter names into English using mapping given by
% spm_DesMtx('ParMap',aMap)
% Pnames  - string mtx/cellstr containing parameter names
% aMap    - (optional) Mx2 cellstr of additional or over-ride mappings
% EPnames - cellstr of converted parameter names
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
% factor. The integer factor levels need not be positive, nor in
% order.  In the '~' constraint types (below), a factor level of zero
% is ignored (treated as no effect), and no corresponding column of
% design matrix is created.  Effects for the factor levels are entered
% into the design matrix *in increasing order* of the factor levels.
% Check Pnames to find out which columns correspond to which levels of
% the factor.
%
% TWO WAY INTERACTIONS: For a two way interaction effect between two
% factors, the FCLevels matrix is an nx2 integer matrix whose columns
% indicate the levels of the two factors. An effect is included for
% each unique combination of the levels of the two factors. Again,
% factor levels must be integer, though not necessarily positive.
% Zero levels are ignored in the '~' constraint types described below.
%
% CONSTRAINTS: Each FactorLevels vector/matrix may be followed by an 
% (optional) ConstraintString.
%
% ConstraintStrings for main effects are:
%                  '-'   - No Constraint
%                  '~'   - Ignore zero level of factor
%                          (I.e. cornerPoint constraint on zero level,
%                          (same as '.0', except zero level is always ignored,
%                          (even if factor only has zero level, in which case
%                          (an empty DesMtx results and a warning is given
%                  '+0'  - sum-to-zero constraint
%                  '+0m' - Implicit sum-to-zero constraint
%                  '.'   - CornerPoint constraint
%                  '.0'  - CornerPoint constraint applied to zero factor level
%                          (warns if there is no zero factor level)
%  Constraints for two way interaction effects are
%    '-'                 - No Constraints
%                  '~'   - Ignore zero level of any factor
%                          (I.e. cornerPoint constraint on zero level,
%                          (same as '.ij0', except zero levels are always ignored
%    '+i0','+j0','+ij0'  - sum-to-zero constraints
%    '.i', '.j', '.ij'   - CornerPoint constraints
%    '.i0','.j0','.ij0'  - CornerPoint constraints applied to zero factor level
%                          (warns if there is no zero factor level)
%    '+i0m', '+j0m'      - Implicit sum-to-zero constraints
%
% With the exception of the "ignore zero" '~' constraint, constraints
% are only applied if there are sufficient factor levels. CornerPoint
% and explicit sum-to-zero Constraints are applied to the last level of
% the factor.
%
% The implicit sum-to-zero constraints "mean correct" appropriate rows
% of the relevant design matrix block. For a main effect, constraint
% '+0m' "mean corrects" the main effect block across columns,
% corresponding to factor effects B_i, where B_i = B'_i - mean(B'_i) :
% The B'_i are the fitted parameters, effectively *relative* factor
% parameters, relative to their mean. This leads to a rank deficient
% design matrix block. If Matlab's pinv, which implements a
% Moore-Penrose pseudoinverse, is used to solve the least squares
% problem, then the solution with smallest L2 norm is found, which has
% mean(B'_i)=0 provided the remainder of the design is unique (design
% matrix blocks of full rank). In this case therefore the B_i are
% identically the B'_i - the mean correction imposes the constraint.
%      
%
% COVARIATES: The FCLevels matrix here is an nxc matrix whose columns
% contain the covariate values. An effect is included for each covariate.
% Covariates are identified by ConstraintString 'C'.
%
%
% PRE-SPECIFIED DESIGN BLOCKS: ConstraintString 'X' identifies a
% ready-made bit of design matrix - the effect is the same as 'C'.
%
%
% FACTOR BY COVARIATE INTERACTIONS: are identified by ConstraintString
% 'FxC'. The last column is understood to contain the covariate. Other
% columns are taken to contain integer FactorLevels vectors. The
% (unconstrained) interaction of the factors is interacted with the
% covariate. Zero factor levels are ignored if ConstraintString '~FxC'
% is used.
%
%
% NAMES: Each Factor/Covariate can be 'named', by passing a name
% string.  Pass a string matrix, or cell array (vector) of strings,
% with rows (cells) naming the factors/covariates in the respective
% columns of the FCLevels matrix.  These names default to <Fac>, <Cov>,
% <Fac1>, <Fac2> &c., and are used in the construction of the Pnames
% parameter names.
% E.g. for an interaction, spm_DesMtx([F1,F2],'+ij0',['subj';'cond'])
% giving parameter names such as subj*cond_{1,2} etc...
%
% Pnames returns a string matrix whose successive rows describe the
% effects parameterised in the corresponding columns of the design
% matrix. `Fac1*Fac2_{2,3}' would refer to the parameter for the
% interaction of the two factors Fac1 & Fac2, at the 2nd level of the
% former and the 3rd level of the latter. Other forms are
%  - Simple main effect (level 1)        : <Fac>_{1}
%  - Three way interaction (level 1,2,3) : <Fac1>*<Fac2>*<Fac3>_{1,2,3}
%  - Two way factor interaction by covariate interaction :
%                                        : <Cov>@<Fac1>*<Fac2>_{1,1}
%  - Column 3 of prespecified DesMtx block (if unnamed)
%                                        : <X> [1]
% The special characters `_*()[]{}' are recognised by the scaling
% function (spm_DesMtx('sca',...), and should therefore be avoided
% when naming effects and covariates.
%
%
% INDEX: An Integer Index matrix is returned if only a single block of
% design matrix is being computed (single set of parameters). It
% indexes the actual order of the effect levels in the design matrix block.
% (Factor levels are introduced in order, regardless of order of
% appearence in the factor index matrices, so that the parameters
% vector has a sensible order.) This is used to aid recursion.
%
% Similarly idx,jdx & kdx are indexes returned for a single block of
% design matrix consisting of unconstrained factor effects ('-' or '~').
% These indexes map I and Index (in a similar fashion to the `unique`
% function) as follows:
%  - idx & jdx are such that I = Index(:,jdx)' and  Index = I(idx,:)'
%    where vector I is given as a column vector
%  - If the "ignore zeros" constraint '~' is used, then kdx indexes the
%    non-zero (combinations) of factor levels, such that
%                     I(kdx,:) = Index(:,jdx)' and  Index == I(kdx(idx),:)'
%
%                           ----------------
%
% The "patterned data setting" (spm_DesMtx('pds'...) is a simple
% utility for setting patterned indicator vectors, inspired by
% MINITAB's "SET" command.
%
% The vector v has it's elements replicated m times, and the resulting
% vector is itself repeated n times, giving a resultant vector i of
% length n*m*length(v)
%
% Examples:
%     spm_DesMtx('pds',1:3)       % = [1,2,3]
%     spm_DesMtx('pds',1:3,2)     % = [1,1,2,2,3,3]
%     spm_DesMtx('pds',1:3,2,3)   % = [1,1,2,2,3,3,1,1,2,2,3,3,1,1,2,2,3,3]
% NB: MINITAB's "SET" command has syntax n(v)m:
%
%                           ----------------
%
% The design matrix scaling feature is designed to return a scaled
% version of a design matrix, with values in [-1,1], suitable for
% visualisation. Special care is taken to apply the same normalisation
% to blocks of design matrix reflecting a single effect, to preserve
% appropriate relationships between columns. Identification of effects
% corresponding to columns of design matrix portions is via the parameter
% names matrices. The design matrix may be passed in any number of
% parts, provided the corresponding parameter names are given. It is
% assummed that the block representing an effect is contained within a
% single partition. Partitions supplied without corresponding parameter
% names are scaled on a column by column basis, the parameters labelled as
% <UnSpec> in the returned nPnames matrix.
% 
% Effects are identified using the special characters `_*()[]{}' used in
% parameter naming as follows: (here ? is a wildcard)
%       - ?(?)          - general  block (column normalised)
%       - ?[?]          - specific block (block normalised)
%       - ?_{?}         - main effect or interaction of main effects
%       - ?@?_{?}       - factor by covariate interaction
% Blocks are identified by looking for runs of parameters of the same type
% with the same names: E.g. a block of main effects for factor 'Fac1'
% would have names like Fac1_{?}.
% 
% Scaling is as follows:
%       * fMRI blocks are scaled around zero to lie in [-1,1]
%       * No scaling is carried out if max(abs(tX(:))) is in [.4,1]
%         This protects dummy variables from normalisation, even if
%         using implicit sum-to-zero constraints.
%       * If the block has a single value, it's replaced by 1's
%       * FxC blocks are normalised so the covariate values cover [-1,1]
%         but leaving zeros as zero.
%       * Otherwise, block is scaled to cover [-1,1].
%
%_______________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% Andrew Holmes
% $Id: spm_DesMtx.m 4137 2010-12-15 17:18:32Z guillaume $



%-Parse arguments for recursive construction of design matrices
%=======================================================================
if nargin==0 error('Insufficient arguments'), end

if ischar(varargin{1})
    %-Non-recursive action string usage
    Constraint=varargin{1};
elseif nargin>=2 && ~(ischar(varargin{2}) || iscell(varargin{2}))
    [X1,Pnames1]=spm_DesMtx(varargin{1});
    [X2,Pnames2]=spm_DesMtx(varargin{2:end});
    X=[X1,X2]; Pnames=[Pnames1;Pnames2];
    return
elseif nargin>=3 && ~(ischar(varargin{3}) || iscell(varargin{3}))
    [X1,Pnames1]=spm_DesMtx(varargin{1:2});
    [X2,Pnames2]=spm_DesMtx(varargin{3:end});
    X=[X1,X2]; Pnames=[Pnames1;Pnames2];
    return
elseif nargin>=4
    [X1,Pnames1]=spm_DesMtx(varargin{1:3});
    [X2,Pnames2]=spm_DesMtx(varargin{4:end});
    X=[X1,X2]; Pnames=[Pnames1;Pnames2];
    return
else
    %-If I is a vector, make it a column vector
    I=varargin{1}; if size(I,1)==1, I=I'; end
    %-Sort out constraint and Factor/Covariate name parameters
    if nargin<2, Constraint='-'; else Constraint=varargin{2}; end
    if isempty(I), Constraint='mt'; end
    if nargin<3, FCnames={}; else FCnames=varargin{3}; end
    if char(FCnames), FCnames=cellstr(FCnames); end
end



switch Constraint, case 'mt'                              %-Empty I case
%=======================================================================
X      = [];
Pnames = {};
Index  = [];



case {'C','X'}           %-Covariate effect, or ready-made design matrix
%=======================================================================
%-I contains a covariate (C), or is to be inserted "as is" (X)
X = I;

%-Construct parameter name index
%-----------------------------------------------------------------------
if isempty(FCnames)
    if strcmp(Constraint,'C'), FCnames={'<Cov>'}; else FCnames={'<X>'}; end
end

if length(FCnames)==1 && size(X,2)>1
    Pnames = cell(size(X,2),1);
    for i=1:size(X,2)
        Pnames{i} = sprintf('%s [%d]',FCnames{1},i);
    end
elseif length(FCnames)~=size(X,2)
    error('FCnames doesn''t match covariate/X matrix')
else
    Pnames = FCnames;
end



case {'-(1)','~(1)'}      %-Simple main effect ('~' ignores zero levels)
%=======================================================================
%-Sort out arguments
%-----------------------------------------------------------------------
if size(I,2)>1, error('Simple main effect requires vector index'), end
if any(I~=floor(I)), error('Non-integer indicator vector'), end
if isempty(FCnames), FCnames = {'<Fac>'};
elseif length(FCnames)>1, error('Too many FCnames'), end

nXrows = size(I,1);

% Sort out unique factor levels - ignore zero level in '~(1)' usage
%-----------------------------------------------------------------------
if Constraint(1)~='~'
    [Index,idx,jdx] = unique(I');
    kdx = [1:nXrows];
else
    [Index,idx,jdx] = unique(I(I~=0)');
    kdx             = find(I~=0)';
    if isempty(Index)
        X=[]; Pnames={}; Index=[];
        warning(['factor has only zero level - ',...
            'returning empty DesMtx partition'])
        return
    end
end

%-Set up unconstrained X matrix & construct parameter name index
%-----------------------------------------------------------------------
nXcols = length(Index);

%-Columns in ascending order of corresponding factor level
X      = zeros(nXrows,nXcols);
Pnames = cell(nXcols,1);
for ii=1:nXcols         %-ii indexes i in Index
    X(:,ii) = I==Index(ii);
    %-Can't use: for i=Index, X(:,i) = I==i; end
    % in case Index has holes &/or doesn't start at 1!
    Pnames{ii} = sprintf('%s_{%d}',FCnames{1},Index(ii));
end
%-Don't append effect level if only one level
if nXcols==1, Pnames=FCnames; end


case {'-','~'}     %-Main effect / interaction ('~' ignores zero levels)
%=======================================================================
if size(I,2)==1
    %-Main effect - process directly
    [X,Pnames,Index,idx,jdx,kdx] = spm_DesMtx(I,[Constraint,'(1)'],FCnames);
    return
end

if any((I(:))~=floor(I(:))), error('Non-integer indicator vector'), end

% Sort out unique factor level combinations & build design matrix
%-----------------------------------------------------------------------
%-Make "raw" index to unique effects
nI     = I - ones(size(I,1),1)*min(I);
tmp    = max(I)-min(I)+1;
tmp    = [fliplr(cumprod(tmp(end:-1:2))),1];
rIndex = sum(nI.*(ones(size(I,1),1)*tmp),2)+1;

%-Ignore combinations where any factor has level zero in '~' usage
if Constraint(1)=='~'
    rIndex(any(I==0,2))=0;
    if all(rIndex==0)
        X=[]; Pnames={}; Index=[];
        warning(['no non-zero factor level combinations - ',...
            'returning empty DesMtx partition'])
        return
    end
end

%-Build design matrix based on unique factor combinations
[X,null,sIndex,idx,jdx,kdx]=spm_DesMtx(rIndex,[Constraint,'(1)']);

%-Sort out Index matrix
Index = I(kdx(idx),:)';

%-Construct parameter name index
%-----------------------------------------------------------------------
if isempty(FCnames)
    tmp = ['<Fac1>',sprintf('*<Fac%d>',2:size(I,2))];
elseif length(FCnames)==size(I,2)
    tmp = [FCnames{1},sprintf('*%s',FCnames{2:end})];
else
    error('#FCnames mismatches #Factors in interaction')
end

Pnames = cell(size(Index,2),1);
for c = 1:size(Index,2)
    Pnames{c} = ...
    [sprintf('%s_{%d',tmp,Index(1,c)),sprintf(',%d',Index(2:end,c)),'}'];
end



case {'FxC','-FxC','~FxC'}           %-Factor dependent covariate effect
%                                       ('~' ignores zero factor levels)
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
elseif length(FCnames)==size(I,2)
    Fnames = FCnames(1:end-1);
    Cnames = FCnames{end};
else
    error('#FCnames mismatches #Factors+#Cov in FxC')
end

%-Set up design matrix X & names matrix - ignore zero levels if '~FxC' use
%-----------------------------------------------------------------------
if Constraint(1)~='~',  [X,Pnames,Index] = spm_DesMtx(F,'-',Fnames);
    else       [X,Pnames,Index] = spm_DesMtx(F,'~',Fnames); end
X = X.*(C*ones(1,size(X,2)));
Pnames = cellstr([repmat([Cnames,'@'],size(Index,2),1),char(Pnames)]);



case {'.','.0','+0','+0m'}              %-Constrained simple main effect
%=======================================================================

if size(I,2)~=1, error('Simple main effect requires vector index'), end

[X,Pnames,Index] = spm_DesMtx(I,'-(1)',FCnames);

%-Impose constraint if more than one effect
%-----------------------------------------------------------------------
%-Apply uniqueness constraints ('.' & '+0')  to last effect, which is
% in last column, since column i corresponds to level Index(i)
%-'.0' corner point constraint is applied to zero factor level only
nXcols = size(X,2);
zCol   = find(Index==0);
if nXcols==1 && ~strcmp(Constraint,'.0')
    error('only one level: can''t constrain')
elseif strcmp(Constraint,'.')
    X(:,nXcols)=[]; Pnames(nXcols)=[]; Index(nXcols)=[];
elseif strcmp(Constraint,'.0')
    zCol = find(Index==0);
    if isempty(zCol),   warning('no zero level to constrain')
    elseif nXcols==1,   error('only one level: can''t constrain'), end
    X(:,zCol)=[];   Pnames(zCol)=[]; Index(zCol)=[];
elseif strcmp(Constraint,'+0')
    X(find(X(:,nXcols)),:)=-1;
    X(:,nXcols)=[]; Pnames(nXcols)=[]; Index(nXcols)=[];
elseif strcmp(Constraint,'+0m')
    X = X - 1/nXcols;
end



case {'.i','.i0','.j','.j0','.ij','.ij0','+i0','+j0','+ij0','+i0m','+j0m'}
                                           %-Two way interaction effects
%=======================================================================
if size(I,2)~=2, error('Two way interaction requires Nx2 index'), end

[X,Pnames,Index] = spm_DesMtx(I,'-',FCnames);

%-Implicit sum to zero
%-----------------------------------------------------------------------
if any(strcmp(Constraint,{'+i0m','+j0m'}))
    SumIToZero = strcmp(Constraint,'+i0m');
    SumJToZero = strcmp(Constraint,'+j0m');

    if SumIToZero   %-impose implicit SumIToZero constraints
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

    if SumJToZero   %-impose implicit SumJToZero constraints
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

    if SumIToZero   %-impose explicit SumIToZero constraints
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
        X(:,cols)=[]; Pnames(cols)=[]; Index(:,cols)=[];
    end

    if SumJToZero   %-impose explicit SumJToZero constraints
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
        X(:,cols)=[]; Pnames(cols)=[]; Index(:,cols)=[];
    end

%-Corner point constraints
%-----------------------------------------------------------------------
elseif any(strcmp(Constraint,{'.i','.i0','.j','.j0','.ij','.ij0'}))
    CornerPointI = any(strcmp(Constraint,{'.i','.i0','.ij','.ij0'}));
    CornerPointJ = any(strcmp(Constraint,{'.j','.j0','.ij','.ij0'}));

    if CornerPointI %-impose CornerPointI constraints
        if Constraint(end)~='0',    i = max(Index(1,:));
            else           i = 0; end
        cols=find(Index(1,:)==i); %-columns to delete
        if isempty(cols)
            warning('no zero i level to constrain')
        elseif all(Index(1,:)==i)
            error('only one i level: can''t constrain')
        end
        %-delete columns
        X(:,cols)=[]; Pnames(cols)=[]; Index(:,cols)=[];
    end

    if CornerPointJ %-impose CornerPointJ constraints
        if Constraint(end)~='0',    j = max(Index(2,:));
            else           j = 0; end
        cols=find(Index(2,:)==j);
        if isempty(cols)
            warning('no zero j level to constrain')
        elseif all(Index(2,:)==j)
            error('only one j level: can''t constrain')
        end
        X(:,cols)=[]; Pnames(cols)=[]; Index(:,cols)=[];
    end
end


case {'PDS','pds'}                          %-Patterned data set utility
%=======================================================================
% i = spm_DesMtx('pds',v,m,n)
if nargin<4, n=1; else n=varargin{4}; end
if nargin<3, m=1; else m=varargin{3}; end
if nargin<2, varargout={[]}; return, else v=varargin{2}; end
if any([size(n),size(m)])>1, error('n & m must be scalars'), end
if any(([m,n]~=floor([m,n]))|([m,n]<1))
    error('n & m must be natural numbers'), end
if sum(size(v)>1)>1, error('v must be a vector'), end

%-Computation
%-----------------------------------------------------------------------
si = ones(1,ndims(v)); si(find(size(v)>1))=n*m*length(v);
X = reshape(repmat(v(:)',m,n),si);



case {'Sca','sca'}                            %-Scale DesMtx for imaging
%=======================================================================
nX   = []; nPnames = {}; Carg = 2;

%-Loop through the arguments accumulating scaled design matrix nX
%-----------------------------------------------------------------------
while(Carg <= nargin)
    rX = varargin{Carg}; Carg=Carg+1;
    if Carg<=nargin && ~isempty(varargin{Carg}) && ...
            (ischar(varargin{Carg}) || iscellstr(varargin{Carg}))
    rPnames = char(varargin{Carg}); Carg=Carg+1;
    else    %-No names to work out blocks from - normalise by column
    rPnames = repmat('<UnSpec>',size(rX,2),1);
    end
    %-Pad out rPnames with 20 spaces to permit looking past line ends
    rPnames = [rPnames,repmat(' ',size(rPnames,1),20)];


    while(~isempty(rX))
    if size(rX,2)>1 && max([1,find(rPnames(1,:)=='(')]) < ...
                    max([0,find(rPnames(1,:)==')')])
    %-Non-specific block: find the rest & column normalise round zero
    %===============================================================
        c1 = max(find(rPnames(1,:)=='('));
        d  = any(diff(abs(rPnames(:,1:c1))),2)...
            | ~any(rPnames(2:end,c1+1:end)==')',2);
        t  = min(find([d;1]));

        %-Normalise columns of block around zero
        %-------------------------------------------------------
        tmp = size(nX,2);
        nX  = [nX, zeros(size(rX,1),t)];
        for i=1:t
            if ~any(rX(:,i))
                nX(:,tmp+i) = 0;
            else
                nX(:,tmp+i) = rX(:,i)/max(abs(rX(:,i)));
            end
        end
        nPnames   = [nPnames; cellstr(rPnames(1:t,:))];
        rX(:,1:t) = []; rPnames(1:t,:)=[];


    elseif size(rX,2)>1 && max([1,find(rPnames(1,:)=='[')]) < ...
                    max([0,find(rPnames(1,:)==']')])
    %-Block: find the rest & normalise together
    %===============================================================
        c1 = max(find(rPnames(1,:)=='['));
        d  = any(diff(abs(rPnames(:,1:c1))),2)...
            | ~any(rPnames(2:end,c1+1:end)==']',2);
        t  = min(find([d;1]));

        %-Normalise block
        %-------------------------------------------------------
        nX = [nX,sf_tXsca(rX(:,1:t))];
        nPnames  = [nPnames; cellstr(rPnames(1:t,:))];
        rX(:,1:t) = []; rPnames(1:t,:)=[];


    elseif size(rX,2)>1 && max([1,strfind(rPnames(1,:),'_{')]) < ...
                    max([0,find(rPnames(1,:)=='}')])
    %-Factor, interaction of factors, or FxC: find the rest...
    %===============================================================
        c1 = max(strfind(rPnames(1,:),'_{'));
        d  = any(diff(abs(rPnames(:,1:c1+1))),2)...
            | ~any(rPnames(2:end,c1+2:end)=='}',2);
        t  = min(find([d;1]));

        %-Normalise block
        %-------------------------------------------------------
        tX = rX(:,1:t);
        if any(rPnames(1,1:c1)=='@')    %-FxC interaction
            C         = tX(tX~=0);
            tX(tX~=0) = 2*(C-min(C))/max(C-min(C))-1;
            nX        = [nX,tX];
        else                %-Straight interaction
            nX = [nX,sf_tXsca(tX)];
        end
        nPnames  = [nPnames; cellstr(rPnames(1:t,:))];
        rX(:,1:t) = []; rPnames(1:t,:)=[];


    else                              %-Dunno! Just column normalise
    %===============================================================
        nX       = [nX,sf_tXsca(rX(:,1))];
        nPnames  = [nPnames; cellstr(rPnames(1,:))];
        rX(:,1)  = []; rPnames(1,:)=[];

    end
    end
end

X      = nX;
Pnames = nPnames;


case {'Fnames','fnames'}     %-Turn parameter names into valid filenames
%=======================================================================
% Fnames = spm_DesMtx('FNames',Pnames)
if nargin<2, varargout={''}; return, end
Fnames = varargin{2};
for i=1:numel(Fnames)
    str = Fnames{i};
    str(str==',')='x';          %-',' to 'x'
    str(str=='*')='-';          %-'*' to '-'
    str(str=='@')='-';          %-'@' to '-'
    str(str==' ')='_';          %-' ' to '_'
    str(str=='/')='';           %- delete '/'
    str(str=='.')='';           %- delete '.'
    Fnames{i} = str;
end
Fnames = spm_str_manip(Fnames,'v');     %- retain only legal characters
X = Fnames;


case {'TeXnames','texnames'}   %-Remove '@' & '*' for TeX interpretation
%=======================================================================
% TPnames = spm_DesMtx('TeXnames',Pnames)
if nargin<2, varargout={''}; return, end
TPnames = varargin{2};
for i=1:prod(size(TPnames))
    str = TPnames{i};
    str(str=='*')='';           %- delete '*'
    str(str=='@')='';           %- delete '@'
    TPnames{i} = str;
end
X = TPnames;


case {'ParMap','parmap'}          %-Parameter mappings: greek to english
%=======================================================================
% Map = spm_DesMtx('ParMap',aMap)
Map = { '\mu',      'const';...
    '\theta',   'repl';...
    '\alpha',   'cond';...
    '\gamma',   'subj';...
    '\rho',     'covint';...
    '\zeta',    'global';...
    '\epsilon', 'error'};
if nargin<2, aMap={}; else aMap = varargin{2}; end
if isempty(aMap), X=Map; return, end
if ~(iscellstr(aMap) && ndims(aMap)==2), error('aMap must be an nx2 cellstr'), end
for i=1:size(aMap,1)
    j = find(strcmp(aMap{i,1},Map(:,1)));
    if isempty(j)
        Map=[aMap(i,:); Map];
    else
        Map(j,2) = aMap(i,2);
    end
end
X = Map;


case {'ETeXNames','etexnames'} %-Search & replace: for Englishifying TeX
%=======================================================================
% EPnames = spm_DesMtx('TeXnames',Pnames,aMap)
if nargin<2, varargout={''}; return, end
if nargin<3, aMap={}; else aMap = varargin{3}; end
Map = spm_DesMtx('ParMap',aMap);
EPnames = varargin{2};
for i=1:size(Map,1)
    EPnames = strrep(EPnames,Map{i,1},Map{i,2});
end
X = EPnames;


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
function nX = sf_tXsca(tX)
if nargin==0, nX=[]; return, end
if abs(max(abs(tX(:)))-0.7)<(.3+1e-10)
    nX = tX;
elseif all(tX(:)==tX(1))
    nX = ones(size(tX));
elseif max(abs(tX(:)))<1e-10
    nX = zeros(size(tX));
else
    nX = 2*(tX-min(tX(:)))/max(tX(:)-min(tX(:)))-1;
end

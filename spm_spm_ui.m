function varargout=spm_spm_ui(varargin)
%_______________________________________________________________________
% %W% Andrew Holmes, Karl Friston %E%
SCCSid  = '%I%';

%-Condition arguments
%-----------------------------------------------------------------------
if (nargin==0), Action = 'CFG'; else, Action = varargin{1}; end


switch lower(Action), case 'cfg'
%=======================================================================
% - C O N F I G U R E   D E S I G N
%=======================================================================
% spm_spm_ui('CFG',D)
if nargin<2, D=[]; else, D=varargin{2}; end

%-GUI setup
%-----------------------------------------------------------------------
SPMid = spm('FnBanner',mfilename,SCCSid);
[Finter,Fgraph,CmdLine] = spm('FnUIsetup','Stats: Setup analysis',0);
spm_help('!ContextHelp',mfilename)


%-Option definitions
%-----------------------------------------------------------------------
%-Generic factor names
sF = {'sF1','sF2','sF3','sF4'};

%-Covariate by factor interaction options
sCFI = {'<none>';...							%-1
	'with sF1';'with sF2';'with sF3';'with sF4';...			%-2:5
	'with sF2 (within sF4)';'with sF3 (within sF4)'};		%-6,7

%-DesMtx argument components for covariate by factor interaction options
% (Used for CFI's Covariate Centering (CC), GMscale & Global normalisation)
CFIforms = {	'[]',		'C',	'{}';...			%-1
		'I(:,1)',	'FxC',	'{D.sF{1}}';...			%-2
		'I(:,2)',	'FxC',	'{D.sF{2}}';...			%-3
		'I(:,3)',	'FxC',	'{D.sF{3}}';...			%-4
		'I(:,4)',	'FxC',	'{D.sF{4}}';...			%-5
		'I(:,[4,2])',	'FxC',	'{D.sF{4},D.sF{2}}';...		%-6
		'I(:,[4,3])',	'FxC',	'{D.sF{4},D.sF{3}}'	};	%-7

%-Centre (mean correction) options for covariates & globals            (CC)
% (options 9-12 are for centering of global when using AnCova GloNorm) (GC)
sCC = {		'around overall mean';...				%-1
		'around sF1 means';...					%-2
		'around sF2 means';...					%-3
		'around sF3 means';...					%-4
		'around sF4 means';...					%-5
		'around sF2 (within sF4) means';...			%-6
		'around sF3 (within sF4) means';...			%-7
		'<no centering>';...					%-8
		'around user specified value';...			%-9
		'(as implied by AnCova)';...				%-10
		'GM';...						%-11
		'(redundant: not doing AnCova)'}';			%-12
%-DesMtx I forms for covariate centering options
CCforms = {'ones(nScan,1)',CFIforms{2:end,1},''}';


%-Global normalization options (options 1-7 match CFIforms)       (GloNorm)
sGloNorm = {	'AnCova';...						%-1
		'AnCova by sF1';...					%-2
		'AnCova by sF2';...					%-3
		'AnCova by sF3';...					%-4
		'AnCova by sF4';...					%-5
		'AnCova by sF2 (within sF4)';...			%-6
		'AnCova by sF3 (within sF4)';...			%-7
		'proportional scaling';...				%-8
		'<no global normalisation>'};				%-9

%-Grand mean scaling options                                        (GMsca)
sGMsca = {	'scaling of overall grand mean';...			%-1
		'scaling of sF1 grand means';...			%-2
		'scaling of sF2 grand means';...			%-3
		'scaling of sF3 grand means';...			%-4
		'scaling of sF4 grand means';...			%-5
		'scaling of sF2 (within sF4) grand means';...		%-6
		'scaling of sF3 (within sF4) grand means';...		%-7
		'(implicit in PropSca global normalisation)';...	%-8
		'<no grand Mean scaling>'	};			%-9
%-NB: Grand mean scaling by subject is redundent for proportional scaling


%-Global calculation options                                       (GXcalc)
sGXcalc  = {	'omit';...						%-1
		'user specified';...					%-2
		'mean voxel value (within per image fullmean/8 mask)'};	%-3


%-Variable "decoder"
%-----------------------------------------------------------------------
%-****


%=======================================================================
%-D E S I G N   P A R A M E T E R S
%=======================================================================
%-Get design type
%-----------------------------------------------------------------------
if isempty(D)
	tmp = spm_input('Select design class...',1,'m',...
		{'Basic stats','SPM98 PET designs','SPM96 PET designs'});
	switch tmp
	case 1,	D = spm_spm_ui('DesDefs_Stats');
	case 2,	D = spm_spm_ui('DesDefs_SPM98PET');
	case 3,	D = spm_spm_ui('DesDefs_SPM96PET');
	otherwise, error('Don''t know that one!')
	end
end

D = D(spm_input('Select design type...',1,'m',{D.DesName}'));

%-Set factor names for this design
%-----------------------------------------------------------------------
sCC       = sf_estrrep(sCC,[sF',D.sF']);
sCFI      = sf_estrrep(sCFI,[sF',D.sF']);
sGloNorm  = sf_estrrep(sGloNorm,[sF',D.sF']);
sGMsca    = sf_estrrep(sGMsca,[sF',D.sF']);

%-Get filenames & factor indicies
%-----------------------------------------------------------------------
[P,I] = spm_spm_ui('Files&Indices',D.sF,D.n,D.b.aTime);	%-Files & indices
nScan = length(P);					%-#observations

%-Additional design parameters
%-----------------------------------------------------------------------
bL  = any(diff(I,1),1); 	%-Multiple factor levels?
	% NB: bL(2) might be thrown by user specified f1 levels
	%     (D.b.aTime & D.n(2)>1) - assumme user is consistent?
bFI = [bL(1),bL(2:3)&~bL(4),bL(4),bL([2,3])&bL(4)];
	%-Allowable interactions for covariates
	%-Only offer interactions with multi-level factors, and
	% don't offer by F2|F3 if bL(4)!

%-Build Condition (H) and Block (B) partitions
%=======================================================================
eval(['[H,Hnames] = spm_DesMtx(',D.Hform,');'])
if rank(H)==nScan, error('unestimable condition effects'), end
eval(['[B,Bnames] = spm_DesMtx(',D.Bform,');'])
if rank(B)==nScan, error('unestimable block effects'), end

%-Drop a constant H partition if B partition can model constant
if size(H,2)>0 & all(H(:)==1) & (rank([H B])==rank(B))
	H = []; Hnames = {};
	warning('Dropping redundant constant H partition')
end


%-Covariate partition(s): interest (C) & nuisance (G) excluding global
%=======================================================================
nC = D.nC;			%-Default #covariates
C  = {[],[]}; Cnames = {{},{}};	%-Covariate DesMtx partitions & names
xC = [];			%-Struct array to hold raw covariates
				% Fields:{'c','cname','iCC','iCFI','type','cols'}

dcname = {'CovInt','NusCov'};	%-Default root names for covariates
dstr   = {'covariate','nuisance variable'};

GUIpos = spm_input('!NextPos');
nc     = [0,0];
for i=1:2			% 1:covariates of interest, 2:nuisance variables
    if isinf(nC(i)), nC(i)=spm_input(['# ',dstr{i},'s'],GUIpos,'w1'); end

    while nc(i) < nC(i)

	%-Create prompt, get covariate, get covariate name
        %---------------------------------------------------------------
	if nC(i)==1, str=dstr{i}; else, str=sprintf('%s %d',dstr{i},nc(i)+1); end
        c = spm_input(str,GUIpos,'r',[],[nScan,Inf]);
        if any(isnan(c(:))), break, end		%-NaN is dummy value to exit
	nc(i)  = nc(i)+1;			%-#Covariates (so far)
	if nC(i)>1,	tstr = sprintf('%s^{%d}',dcname{i},nc(i));
	else,		tstr = dcname{i}; end
       	cname  = spm_input([str,' name?'],'+1','s',tstr);
       	rc     = c;				%-Save covariate value
	rcname = cname;				%-Save covariate name

        %-Interaction option? (if single covariate vector entered)?
        %---------------------------------------------------------------
        if size(c,2) == 1
       	    if length(D.iCFI{i})>1
       		%-User choice of interaction options, default is negative
       		%-Only offer interactions for appropriate factor combinations
		iCFI = intersect(abs(D.iCFI{i}),find([1,bFI]));
		dCFI = max([1,intersect(iCFI,-D.iCFI{i}(D.iCFI{i}<0))]);
        	iCFI = spm_input([str,': interaction?'],'+1','m',...
			sCFI(iCFI),iCFI,find(iCFI==dCFI));
	    else
		iCFI = abs(D.iCFI{i});		%-AutoSelect default option
	    end
	else
	    iCFI = 1;
	end

        %-Centre covariate(s)? (Default centring to correspond to CFI)
        % Always offer "no centering" as default for design matrix blocks
        %---------------------------------------------------------------
	DiCC = D.iCC{i};
	if size(c,2)>1, DiCC = union(DiCC,-8); end
        if length(DiCC)>1
        	%-User has a choice of centering options
		%-Only offer factor specific for appropriate factor combinations
		iCC = intersect(abs(DiCC),find([1,bFI,1]) );
        	%-Default is max -ve option in D, overridden by iCFI if CFI
		if iCFI==1, dCC=-DiCC(DiCC<0); else, dCC=iCFI; end
		dCC = max([1,intersect(iCC,dCC)]);
		iCC = spm_input([str,': centre?'],'+1','m',...
			sCC(iCC),iCC,find(iCC==dCC));
        else
        	iCC = abs(DiCC);	%-AutoSelect default option
        end
	%-Centre within factor levels as appropriate
        if any(iCC==[1:7]), c = c - spm_meanby(c,eval(CCforms{iCC})); end

        %-Do any interaction (only for single covariate vectors)
        %---------------------------------------------------------------
        if iCFI>1				%-(NB:iCFI=1 if size(c,2)>1)
       		tI        = [eval(CFIforms{iCFI,1}),c];
		tConst    = CFIforms{iCFI,2};
		tFnames   = [eval(CFIforms{iCFI,3}),{cname}];
		[c,cname] = spm_DesMtx(tI,tConst,tFnames);
	elseif size(c,2)>1			%-Design matrix block
		[null,cname] = spm_DesMtx(c,'X',cname);
	else
		cname = {cname};
	end

	%-Store raw covariate details in xC struct for reference
	%-Pack c into appropriate DesMtx partition
        %---------------------------------------------------------------
	%-Construct description string for covariate
	str = {sprintf('%s: %s',str,rcname)};
	if size(rc,2)>1, str = {sprintf('%s (block of %d covariates)',...
		str{:},size(rc,2))}; end
	if iCC<8, str=[str;{['used centered ',sCC{iCC}]}]; end
	if iCFI>1, str=[str;{['fitted as interaction ',sCFI{iCFI}]}]; end

	tmp       = struct(	'rc',rc,	'rcname',rcname,...
				'c',c,		'cname',{cname},...
				'iCC',iCC,	'iCFI',iCFI,...
				'type',i,...
				'cols',[1:size(c,2)] + ...
						size([H,C{1}],2) +  ...
						size([B,C{2}],2)*(i-1),...
				'descrip',{str}				);
	if isempty(xC), xC=tmp; else, xC=[xC,tmp]; end
	C{i}      = [C{i},c];
	Cnames{i} = [Cnames{i}; cname];

    end	% (while)

end % (for)
clear tI tConst tFnames
spm_input('!SetNextPos',GUIpos);

%-Unpack into C & G design matrix sub-partitions
G = C{2}; Gnames = Cnames{2};
C = C{1}; Cnames = Cnames{1};


%-Options...
%=======================================================================
%-Global normalization options                                 (GloNorm)
%-----------------------------------------------------------------------
if length(D.iGloNorm)>1
	%-User choice of global normalisation options, default is negative
	%-Only offer factor specific for appropriate factor combinations
	iGloNorm = intersect(abs(D.iGloNorm),find([1,bFI,1,1]));
	dGloNorm = max([0,intersect(iGloNorm,-D.iGloNorm(D.iGloNorm<0))]);
	iGloNorm = spm_input('GloNorm: Select global normalisation','+1','m',...
	    	sGloNorm(iGloNorm),iGloNorm,find(iGloNorm==dGloNorm));
else
	iGloNorm = abs(D.iGloNorm);
end


%-Grand mean scaling options                                     (GMsca)
%-----------------------------------------------------------------------
if iGloNorm==8
	iGMsca=8;	%-grand mean scaling implicit in PropSca GloNorm
elseif length(D.iGMsca)==1
	iGMsca = abs(D.iGMsca);
else
	%-User choice of grand mean scaling options
	%-Only offer factor specific for appropriate factor combinations
	iGMsca = intersect(abs(D.iGMsca),find([1,bFI,0,1]));
        %-Default is max -ve option in D, overridden by iGloNorm if AnCova
        if iGloNorm==9, dGMsca=-D.iGMsca(D.iGMsca<0); else, dGMsca=iGloNorm; end
	dGMsca = max([0,intersect(iGMsca,dGMsca)]);
	iGMsca = spm_input('GMsca: grand mean scaling','+1','m',...
	    	sGMsca(iGMsca),iGMsca,find(iGMsca==dGMsca));
end


%-Value for PropSca / GMsca                                         (GM)
%-----------------------------------------------------------------------
if iGMsca==9                            %-Not scaling (GMsca or PropSca)
    GM=0;                               %-Set GM to zero when not scaling
else                                    %-Ask user value of GM
	if iGloNorm==8
		str='PropSca global mean to';
	else
		str=[strrep(sGMsca{iGMsca},'scaling of','scale'),' to'];
	end
	GM = spm_input(str,'+1','r',D.GM,1);
	%-If GM is zero then don't GMsca! or PropSca GloNorm
	if GM==0, iGMsca=9; if iGloNorm==8, iGloNorm=9; end, end
end

%-Sort out description strings for GloNorm and GMsca
sGloNorm = sGloNorm{iGloNorm};
sGMsca   = sGMsca{iGMsca};
if iGloNorm==8
	sGloNorm = sprintf('%s to %-4g',sGloNorm,GM);
elseif iGMsca<8
	sGMsca = sprintf('%s to %-4g',sGMsca,GM);
end


%-Global centering (for AnCova GloNorm)                             (GC)
%-----------------------------------------------------------------------
%-Specify the centering option for the global covariate for AnCova
%-Basically, if 'GMsca'ling then should centre to GM (iGC=11). Otherwise,
% should centre in similar fashion to AnCova (i.e. by the same factor(s)),
% such that models are seperable (iGC=10). This is particularly important
% for subject specific condition effects if then passed on to a second-level
% model. (See also spm_adjmean_ui.m) SPM96 (& earlier) used to just centre
% GX around its (overall) mean (iGC=1).

%-This code allows more general options to be specified (but is a bit complex)
%-Setting D.iGC=[-10,-11] gives the standard choices above

%-If not doing AnCova then GC is irrelevant
if ~any(iGloNorm==[1:7])
	iGC = 12;
	gc  = [];
else
	%-Annotate options 10 & 11 with specific details
	%---------------------------------------------------------------
	%-Tag '(as implied by AnCova)' with actual AnCova situation
	sCC{10} = [sCC{iGloNorm},' (<= ',sGloNorm,')'];
	%-Tag 'GM' case with actual GM & GMsca case
	sCC{11} = sprintf('around GM=%g (i.e. %s after grand mean scaling)',...
		GM,strrep(sCC{iGMsca},'around ',''));

	%-Constuct vector of allowable iGC
	%---------------------------------------------------------------
	%-Weed out redundent factor combinations from pre-set allowable options
	iGC = intersect(abs(D.iGC),find([1,bFI,1,1,1,1]));
	%-Omit 'GM' option if didn't GMsca (iGMsca~=8 'cos doing AnCova)
	if any(iGMsca==[8,9]), iGC = setdiff(iGC,11); end
	%-Omit 'GM' option if same as '(as implied by AnCova)'
	if iGloNorm==iGMsca, iGC = setdiff(iGC,11); end

	%-If there's a choice, set defaults (if any), & get answer
	%---------------------------------------------------------------
	if length(iGC)>1
		dGC = max([0,intersect(iGC,-D.iGC(D.iGC<0))]);
		str = 'Centre global covariate';
		if iGMsca<8, str = [str,' (after grand mean scaling)']; end
		iGC = spm_input(str,'+1','m',sCC(iGC),iGC,find(iGC==dGC));
	elseif isempty(iGC)
		error('Configuration error: empty iGC')
	end

	%-If 'user specified' then get value
	%---------------------------------------------------------------
	if iGC==9
		gc     = spm_input('Centre globals around','+0','r',D.GM,1);
		sCC{9} = sprintf('%s of %g',sCC{iGC},gc);
	else
		gc  = 0;
	end
end


%-Thresholds & masks defining voxels to analyse                   (MASK)
%=======================================================================
GUIpos = spm_input('!NextPos');

%-Analysis threshold mask
%-----------------------------------------------------------------------
%-Work out available options:
% -Inf=>None, complex=>proportional, real=>absolute (i.e. times global)
M_T = D.M_.T; if isempty(M_T), M_T = [-Inf, 100, 0.8*sqrt(-1)]; end
M_T = {	'none',		M_T(min(find(isinf(M_T))));...
	'absolute',	M_T(min(find(isfinite(M_T)&(M_T==real(M_T)))));...
	'prop''nal',	M_T(min(find(isfinite(M_T)&(M_T~=real(M_T)))))	};

%-Work out available options
q = ~[isempty(M_T{1,2}), isempty(M_T{2,2}), isempty(M_T{3,2})];

%-If there's a choice between proportional and absolute then ask
if all(q(2:3))
	tmp = spm_input('Threshold masking',GUIpos,'b',M_T(q,1),find(q));
	q(setdiff([1:3],tmp))=0;
end

%-Get mask value - note that at most one of q(2:3) is true
if ~any(q)				%-Oops - nothing specified!
	M_T = -Inf;
elseif all(q==[1,0,0])			%-no threshold masking
	M_T = -Inf;
else					%-get mask value
	if q(1),	args = {'br1','None',-Inf,abs(M_T{1+find(q(2:3)),2})};
	else,		args = {'r',abs(M_T{1+find(q(2:3)),2})}; end
	if q(2)
		M_T = spm_input('analysis threshold',GUIpos,args{:});
	elseif q(3)
		M_T = spm_input('analysis thresh  (prop''n of global)',GUIpos,...
								args{:});
		if isfinite(M_T) & isreal(M_T), M_T=M_T*sqrt(-1); end
	else
		error('Shouldn''t get here!')
	end
end

%-Make a description string
if isinf(M_T)
	xsM.Analysis_threshold = 'None (-Inf)';
elseif isreal(M_T)
	xsM.Analysis_threshold = sprintf('images thresholded at %6g',M_T);
else
	xsM.Analysis_threshold = sprintf(['images thresholded at %6g ',...
		'times global'],imag(M_T));
end


%-Implicit masking: Ignore zero voxels in low data-types?
%-----------------------------------------------------------------------
% (Implicit mask is NaN in higher data-types.)
type = getfield(spm_vol(P{1}),'dim')*[0,0,0,1]';
if ~spm_type(type,'nanrep')
	switch D.M_.I
	case Inf,    M_I = spm_input('Implicit mask (ignore zero''s)?',...
			'+1','y/n',[1,0],1);		%-Ask
	case {0,1}, M_I = D.M_.I;			%-Pre-specified
	otherwise,  error('unrecognised D.M_.I type')
	end

	if M_I, xsM.Implicit_masking = 'Yes: zero''s treated as missing';
	else,   xsm.Implicit_masking = 'No'; end
else
	M_I = 1;
	xsM.Implicit_masking = 'Yes: NaN''s treated as missing';
end


%-Explicit mask images (map them later...)
%-----------------------------------------------------------------------
switch(D.M_.X)
case Inf,    M_X = spm_input('explicit mask images?','+1','y/n',[1,0],2);
case {0,1}, M_X = D.M_.X;
otherwise,  error('unrecognised D.M_.X type')
end
if M_X, M_P = spm_get(Inf,'*.img',{'select mask images'}); else, M_P = {}; end


%-Global calculation                                            (GXcalc)
%=======================================================================
iGXcalc = abs(D.iGXcalc);
%-Only offer "omit" option if not doing any GloNorm, GMsca or PropTHRESH
if ~(iGloNorm==9 & iGMsca==9 & (isinf(M_T)|isreal(M_T)))
	iGXcalc = intersect(iGXcalc,[2:size(sGXcalc,1)]);
end
if isempty(iGXcalc)
	error('no GXcalc options')
elseif length(iGXcalc)>1
	%-User choice of global calculation options, default is negative
	dGXcalc = max([1,intersect(iGXcalc,-D.iGXcalc(D.iGXcalc<0))]);
	iGXcalc = spm_input('Global calculation','+1','m',...
	    	sGXcalc(iGXcalc),iGXcalc,find(iGXcalc==dGXcalc));
else
	iGXcalc = abs(D.iGXcalc);
end

if iGXcalc==2				%-Get user specified globals
	g = spm_input('globals','+0','r',[],[nScan,1]);
end
sGXcalc = sGXcalc{iGXcalc};


%=======================================================================
% - C O N F I G U R E   D E S I G N
%=======================================================================
spm('FigName','Stats: configuring',Finter,CmdLine);
spm('Pointer','Watch');


%-Images & image info: Map Y image files and check consistency of
% dimensions and orientation / voxel size
%=======================================================================
VY = spm_vol(char(P));

if any(any(diff(cat(1,VY.dim),1,1),1)&[1,1,1,0]) %NB: Bombs for single image
	error('images do not all have the same dimensions'), end
if any(any(any(diff(cat(3,VY.mat),1,3),3)))
	error('images do not all have same orientation & voxel size'), end


%-Global values, scaling and global normalisation
%=======================================================================
%-Compute global values
%-----------------------------------------------------------------------
switch iGXcalc, case 1
	%-Don't compute => no GMsca (iGMsca==9) or GloNorm (iGloNorm==9)
	g = [];
case 2
	%-User specified globals
case 3
	%-Compute as mean voxel value (within per image fullmean/8 mask)
	g = zeros(nScan,1);
	for i = 1:nScan, g(i) = spm_global(VY(i)); end
otherwise
	error('illegal iGXcalc')
end
rg = g;


%-Scaling: compute global scaling factors gSF required to implement proportional
% scaling global normalisation (PropSca) or grand mean scaling (GMsca),
% as specified by iGMsca (& iGloNorm)
%-----------------------------------------------------------------------
switch iGMsca, case 8
	%-Proportional scaling global normalisation
	if iGloNorm~=8, error('iGloNorm-iGMsca(8) mismatch for PropSca'), end
	gSF    = GM./g;
	g      = GM*ones(nScan,1);
case {1,2,3,4,5,6,7}
	%-Grand mean scaling according to iGMsca
	gSF    = GM./spm_meanby(g,eval(CCforms{iGMsca}));
	g      = g.*gSF;
case 9
	%-No grand mean scaling
	gSF    = ones(nScan,1);
otherwise
	error('illegal iGMsca')
end


%-Apply gSF to memory-mapped scalefactors to implement scaling
%-----------------------------------------------------------------------
for i=1:nScan, VY(i).pinfo(1:2,:)=VY(i).pinfo(1:2,:)*gSF(i); end


%-AnCova: Construct global nuisance covariates partition (if AnCova)
%-----------------------------------------------------------------------
if any(iGloNorm==[1:7])

	%-Centre global covariate as requested
	%---------------------------------------------------------------
	switch iGC, case {1,2,3,4,5,6,7}	%-Standard sCC options
		gc = spm_meanby(g,eval(CCforms{iGC}));
	case 8					%-No centering
		gc = 0;
	case 9					%-User specified centre
		%-gc set above
	case 10					%-As implied by AnCova option
		gc = spm_meanby(g,eval(CCforms{iGloNorm}));
	case 11					%-Around GM
		gc = GM;
	otherwise				%-unknown iGC
		error('unexpected iGC value')
	end


	%-AnCova - add scaled centred global to DesMtx `G' partition
	%---------------------------------------------------------------
	tI        = [eval(CFIforms{iGloNorm,1}),g-gc];
	tConst    = CFIforms{iGloNorm,2};
	tFnames   = [eval(CFIforms{iGloNorm,3}),{'global'}];
	[g,gname] = spm_DesMtx(tI,tConst,tFnames);
	clear tI tConst tFnames

	%-Save GX info in xC struct for reference
	str = {sprintf('%s: global',dstr{2})};
	if any(iGMsca==[1:7]), str=[str;{['(after ',sGMsca,')']}]; end
	if iGC~=8, str=[str;{['used centered ',sCC{iGC}]}]; end
	if iGloNorm>1, str=[str;{['fitted as interaction ',sCFI{iGloNorm}]}]; end
	tmp       = struct(	'rc',rg.*gSF,	'rcname','global',...
				'c',g,		'cname',{gname},...
				'iCC',iGC,	'iCFI',iGloNorm,...
				'type',3,...
				'cols',[1:size(g,2)]+size([H C B G],2),...
				'descrip',{str}				);

	G = [G,g]; Gnames = [Gnames; gname];
	if isempty(xC), xC=tmp; else, xC=[xC,tmp]; end

elseif iGloNorm==8 | iGXcalc>1

	if iGloNorm==8
		str = { 'global values: (used for proportional scaling)';...
			'("raw" unscaled globals shown)'};
	elseif isfinite(M_T) & ~isreal(M_T)
		str = { 'global values: (used to compute analysis threshold)'};
	else
		str = { 'global values: (computed but not used)'};
	end

	tmp       = struct(	'rc',rg,	'rcname','global',...
				'c',{[]},	'cname',{{}},...
				'iCC',0,	'iCFI',0,...
				'type',3,...
				'cols',{[]},...
				'descrip',{str}				);
	if isempty(xC), xC=tmp; else, xC=[xC,tmp]; end

end


%-Save info on global calculation in xGX structure
%-----------------------------------------------------------------------
xGX = struct(...
	'iGXcalc',iGXcalc,	'sGXcalc',sGXcalc,	'rg',rg,...
	'iGMsca',iGMsca,	'sGMsca',sGMsca,	'GM',GM,'gSF',gSF,...
	'iGC',	iGC,		'sGC',	sCC{iGC},	'gc',	gc,...
	'iGloNorm',iGloNorm,	'sGloNorm',sGloNorm);



%-Construct masking information structure and compute actual analysis
% threshold using scaled globals (rg.*gSF)
%-----------------------------------------------------------------------
if isreal(M_T),	M_TH =      M_T  * ones(nScan,1);	%-NB: -Inf is real
else,		M_TH = imag(M_T) * (rg.*gSF); end

if ~isempty(M_P)
	VM = spm_vol(char(M_P));
	xsM.Explicit_masking = [{'Yes: mask images :'};{VM.fname}'];
else
	VM={};
	xsM.Explicit_masking = 'No';
end
xM = struct('T',M_T, 'TH',M_TH, 'I',M_I, 'VM',{VM}, 'xs',xsM);



%-Design matrix orthogonalisation & canonicalisation
%=======================================================================
%-raw design matrix
rX = [H,C,B,G];

%-Design matrix orthogonalisation
%-----------------------------------------------------------------------
sXorth = {};
%-Orthogonalise G wirit B if D.Xorth.GwB
if D.Xorth.GwB & ~isempty(G) & ~isempty(B)
	G = G     - B*pinv(B)    *G;
	sXorth = [sXorth, {'nuisance effects wirit constant'}];
end
%-Orthogonalise H wirit [B,G] if D.Xorth.HwBG
if D.Xorth.HwBG & ~isempty(H) & ~isempty([B,G])
	H = H - [B,G]*pinv([B,G])*H;
	sXorth = [sXorth, {'condition effects wirit nuisance effects & const'}];
end
%-Orthogonalise C wirit [B,G] if D.Xorth.CwBG
if D.Xorth.CwBG & ~isempty(C) & ~isempty([B,G])
	C = C - [B,G]*pinv([B,G])*C;
	sXorth = [sXorth, {'covariates of interest wirit nuisance effects & const'}];
end
%-Orthogonalise C wirit H if D.Xorth.CwH
if D.Xorth.CwH & ~isempty(C) & ~isempty(H)
	C = C     - H*pinv(H)    *C;
	sXorth = [sXorth, {'covariates of interest wirit condition effects'}];
end
%-Fine tune sXorth description cellstr
if isempty(sXorth), sXorth={'none'};
else, sXorth=[{'Yes: in order...'};sXorth]; end


%-Construct full design matrix (X), parameter names (Xnames),
% and design information structure (xX)
%-----------------------------------------------------------------------
X      = [H C B G];
Xnames = [Hnames; Cnames; Bnames; Gnames];
tmp    = cumsum([size(H,2), size(C,2), size(B,2), size(G,2)]);
xX     = struct(	'I',		I,...
			'sF',		{D.sF},...
			'rX',		rX,...
			'X',		X,...
			'K',		speye(size(X,1)),...
			'iH',		[1:size(H,2)],...
			'iC',		[1:size(C,2)] + tmp(1),...
			'iB',		[1:size(B,2)] + tmp(2),...
			'iG',		[1:size(G,2)] + tmp(3),...
			'Xnames',	{Xnames},...
			'sXorth',	{sXorth}			);


%-Pre-specified contrast for F-test on effects of interest
%=======================================================================
if ~isempty([H,C])
	%-Have a first guess at a simple F-contrast!
	c = [diff(eye(size(H,2))), zeros(size(H,2)-1,size([C,B,G],2));...
		zeros(size(C,2),size(H,2)), eye(size(C,2)),...
			zeros(size(C,2),size([B,G],2))];

	%-If that's not of the same rank as the hypothesised redundant design
	% subspace, or the rows aren't all contrasts, make a fancy F-contrast
	if ~spm_DesUtil('AllCon',X,c) | size(c,1)~=(rank([H,C,B,G])-rank([B,G]))
		c = spm_DesUtil('FCon',X,size([H,C],2)+[1:size([B,G],2)]);
	end
	
	xCon = struct(	'name',	'no effects of interest',...
			'con',	c,...
			'V',	[]);

%elseif ~isempty(G)
else

	xCon = [];

end


%-Design description (an nx2 cellstr) - for saving and display
%=======================================================================
tmp = {	sprintf('%d condition, +%d covariate, +%d block, +%d nuisance',...
		size(H,2),size(C,2),size(B,2),size(G,2));...
	sprintf('%d total, having %d degrees of freedom',...
		size(X,2),rank(X));...
	sprintf('leaving %d degrees of freedom from %d images',...
		size(X,1)-rank(X),size(X,1))				};
xsDes = struct(	'Design',			{D.DesName},...
		'Global_calculation',		{sGXcalc},...
		'Grand_mean_scaling',		{sGMsca},...
		'Global_normalisation',		{sGloNorm},...
		'Design_Orthogonalisation',	{sXorth},...
		'Parameters',			{tmp}			);

%-Save SPMcfg.mat file
save SPMcfg SPMid D xsDes VY xM xC xX xGX xCon

%-Display Design report
%=======================================================================
spm_spm_ui('DesRepUI',VY,xX,xC,xsDes,xM)

%-Analysis Proper - ****
%=======================================================================
%-This is PET data so sigma = 0 (i.e. independent observations) 
% RT is undefined
%-----------------------------------------------------------------------
%spm_spm(VY,H,C,B,G,CONTRAST,ORIGIN,THRESH*GX,HCBGnames,P,0,[])


%-End: Cleanup GUI
%=======================================================================
spm('FigName','Stats: done',Finter,CmdLine); spm('Pointer','Arrow')
fprintf('\n\n')



case 'files&indices'
%=======================================================================
% - Get files and factor indices
%=======================================================================
% [P,I] = spm_spm_ui('Files&Indices',DsF,Dn,DbaTime)
% DbaTime=D.b.aTime; Dn=D.n; DsF=D.sF;
if nargin<4, DbaTime = 1; else, DbaTime = varargin{4}; end
if nargin<3, Dn  = [Inf,Inf,Inf,Inf]; else, Dn=varargin{3}; end
if nargin<2, DsF = {'Fac1','Fac2','Fac3','Fac4'}; else, DsF=varargin{2}; end

%-Initialise variables
%-----------------------------------------------------------------------
i4 = [];		% factor 4 index (usually study)
i3 = [];		% factor 3 index (usually subject), per f4
i2 = [];		% factor 2 index (usually condition), per f3/f4
i1 = [];		% factor 1 index (usually replication), per f2/f3/f4
P  = {};		% cell array of string filenames

%-Accrue filenames and factor level indicator vectors
%-----------------------------------------------------------------------
if isinf(Dn(4)), n4 = spm_input(['#',DsF{4},'''s'],'+1','n1',[],1);
	else, n4 = Dn(4); end
bL4 = n4>1;

ti2 = '';
GUIpos = spm_input('!NextPos');
for j4  = 1:n4
    spm_input('!SetNextPos',GUIpos);
    sF4P=''; if bL4, sF4P=[DsF{4},' ',int2str(j4),': ']; end
    if isinf(Dn(3)), n3=spm_input([sF4P,'#',DsF{3},'''s'],'+1','n1',[],1);
	    else, n3 = Dn(3); end
    bL3 = n3>1;
    
    if DbaTime & Dn(2)>1
	%disp('NB:selecting in time order - manually specify conditions')
	%-NB: This means f2 levels might not be 1:n2
	GUIpos2 = spm_input('!NextPos');
	for j3 = 1:n3
	    sF3P=''; if bL3, sF3P=[DsF{3},' ',int2str(j3),': ']; end
	    str = [sF4P,sF3P];
	    tP = spm_get(Dn(2)*Dn(1),'.img',{[str,'select images...']});
	    n21 = length(tP);
	    ti2 = spm_input([str,' ',DsF{2},'?'],GUIpos2,'c',ti2',n21,Dn(2));
	    %-Work out i1 & check
	    [tl2,null,j] = unique(ti2);
	    tn1 = zeros(size(tl2)); ti1 = zeros(size(ti2));
	    for i=1:length(tl2)
		    tn1(i)=sum(j==i); ti1(ti2==tl2(i))=1:tn1(i); end
	    if isfinite(Dn(1)) & any(tn1~=Dn(1))
		%-#i1 levels mismatches specification in Dn(1)
		error(sprintf('#%s not %d as pre-specified',DsF{1},Dn(1)))
	    end
	    P   = [P;tP];
	    i4 = [i4; j4*ones(n21,1)];
	    i3 = [i3; j3*ones(n21,1)];
	    i2 = [i2; ti2];
	    i1 = [i1; ti1];
	end

    else

	if isinf(Dn(2))
	    n2 = spm_input([sF4P,'#',DsF{2},'''s'],'+1','n1',[],1);
	else
	    n2 = Dn(2);
	end
	bL2 = n2>1;

	if n2==1 & Dn(1)==1 %-single scan per f3 (subj)
	    %disp('NB:single scan per f3')
	    str = [sF4P,'select images, ',DsF{3},' 1-',int2str(n3)];
	    P     = [P;spm_get(n3,'.img',{str})];
	    i4 = [i4; j4*ones(n3,1)];
	    i3 = [i3; [1:n3]'];
	    i2 = [i2; ones(n3,1)];
	    i1 = [i1; ones(n3,1)];
	else
	    %-multi scan per f3 (subj) case
	    %disp('NB:multi scan per f3')
	    for j3 = 1:n3
		sF3P=''; if bL3, sF3P=[DsF{3},' ',int2str(j3),': ']; end
		if Dn(1)==1
			%-No f1 (repl) within f2 (cond)
			%disp('NB:no f1 within f2')
			str = [sF4P,sF3P,'select images: ',DsF{2},...
				 ' 1-',int2str(n2)];
			P = [P;spm_get(n2,'.img',{str})];
			i4 = [i4; j4*ones(n2,1)];
			i3 = [i3; j3*ones(n2,1)];
			i2 = [i2; [1:n2]'];
			i1 = [i1; ones(n2,1)];
		else
		    %-multi f1 (repl) within f2 (cond)
		    %disp('NB:f1 within f2')
		    for j2 = 1:n2
			sF2P='';
			if bL2, sF2P=[DsF{2},' ',int2str(j2),': ']; end
			str = [sF4P,sF3P,sF2P,' select images...'];
			tP  = spm_get(Dn(1),'.img',{str});
			n1 = size(tP,1);
			P   = [P;tP];
			i4 = [i4; j4*ones(n1,1)];
			i3 = [i3; j3*ones(n1,1)];
			i2 = [i2; j2*ones(n1,1)];
			i1 = [i1; [1:n1]'];
		    end                         % (for j2)
		end                             % (if Dn(1)==1)
	    end                                 % (for j3)
	end                                     % (if  n2==1 &...)
    end                                         % (if DbaTime & Dn(2)>1)
end                                             % (for j4)
varargout = {P,[i1,i2,i3,i4]};


case {'desrep','desrepui'}
%=======================================================================
% - D I S P L A Y   A N A L Y S I S   P A R A M E T E R S
%=======================================================================
% spm_spm_ui('DesRepUI',SPMcfg)
% spm_spm_ui('DesRepUI',VY,xX,xC,xsDes,xM)
if nargin==1
	SPMcfg = load(spm_get(1,'SPMcfg.mat','Select SPMcfg.mat'));
elseif nargin==2
	if isstruct(varargin{2})
		SPMcfg = varargin{2};
	elseif isstr(varargin{2}) & exist(varargin{2},'file')==2
		SPMcfg = load(varargin{2});
	else
		error('Mis-specified SPMcfg')
	end
elseif nargin==6
	SPMcfg = struct(	'xM',		varargin{6},...
				'xsDes',	varargin{5},...
				'xC',		varargin{4},...
				'xX',		varargin{3},...
				'VY',		varargin{2}	);
else
	error('insufficient arguments')
end


%-Action 'DesRep' - print all reports
%-----------------------------------------------------------------------
if strcmp(lower(Action),'desrep')
	%-Files, indices, covariates... & masking options
	spm_DesRep('Files&Factors',{SPMcfg.VY.fname}',SPMcfg.xX.I,SPMcfg.xC,...
		SPMcfg.xX.sF,SPMcfg.xM.xs)
	spm_print
	
	%-Design matrix & design descriptions
	spm_DesRep('DesMtx',SPMcfg.xX.X,SPMcfg.xX.Xnames,...
		{SPMcfg.VY.fname}',SPMcfg.xsDes)
	spm_print
	
	%-Covariates
	spm_DesRep('Covs',SPMcfg.xC,SPMcfg.xX.X,SPMcfg.xX.Xnames)
	spm_print
end


%-GUI setup
%-----------------------------------------------------------------------
Labels = {'Files & factors';'Design matrix';'Covariates'};
cb     = {	['spm_DesRep(''Files&Factors'',',...
			'{UD.VY.fname}'',UD.xX.I,UD.xC,UD.xX.sF,UD.xM.xs)'],...
		['spm_DesRep(''DesMtx'',',...
			'UD.xX.X,UD.xX.Xnames,{UD.VY.fname}'',UD.xsDes)'],...
		['spm_DesRep(''Covs'',',...
			'UD.xC,UD.xX.X,UD.xX.Xnames)']	};
if ~length(SPMcfg.xC), Labels(3)=[]; cb(3)=[]; end
h = spm_input('Select design summary','!_','p',Labels,cb,SPMcfg,1);



case 'desdefs_stats'
%=======================================================================
% - Basic Stats Design definitions...
%=======================================================================
% D = spm_spm_ui('DesDefs_Stats');
% These are the SPM98 basic Stats design definitions...

%-Note: struct expands cell array values to give multiple records:
%       => must embed cell arrays within another cell array!
%-Negative indices indicate defaults (first used)

D = struct(...
	'DesName','One sample t-test',...
	'n',	[Inf 1 1 1],	'sF',{{'obs','','',''}},...
	'Hform',		'I(:,2),''-'',''mean''',...
	'Bform',		'[]',...
	'nC',[0,0],'iCC',{{8,8}},'iCFI',{{1,1}},...
	'Xorth',struct('GwB',0,'HwBG',0,'CwBG',0,'CwH',0),...
	'iGXcalc',[-1,2,3],'iGMsca',[1,-9],'GM',[],...
	'iGloNorm',9,'iGC',12,...
	'M_',struct('T',-Inf,'I',Inf,'X',Inf),...
	'b',struct('aTime',0));

D = [D, struct(...
	'DesName','Two sample t-test',...
	'n',	[Inf 2 1 1],	'sF',{{'obs','group','',''}},...
	'Hform',		'I(:,2),''-'',''group''',...
	'Bform',		'I(:,3),''-'',''\mu''',...
	'nC',[0,0],'iCC',{{8,8}},'iCFI',{{1,1}},...
	'Xorth',struct('GwB',0,'HwBG',0,'CwBG',0,'CwH',0),...
	'iGXcalc',[-1,2,3],'iGMsca',[1,-9],'GM',[],...
	'iGloNorm',9,'iGC',12,...
	'M_',struct('T',[-Inf,0,0.8*sqrt(-1)],'I',Inf,'X',Inf),...
	'b',struct('aTime',1))];

D = [D, struct(...
	'DesName','Paired t-test',...
	'n',	[1 2 Inf 1],	'sF',{{'','cond','pair',''}},...
	'Hform',		'I(:,2),''-'',''condition''',...
	'Bform',		'I(:,3),''-'',''\gamma''',...
	'nC',[0,0],'iCC',{{8,8}},'iCFI',{{1,1}},...
	'Xorth',struct('GwB',0,'HwBG',0,'CwBG',0,'CwH',0),...
	'iGXcalc',[-1,2,3],'iGMsca',[1,-9],'GM',[],...
	'iGloNorm',9,'iGC',12,...
	'M_',struct('T',[-Inf,0,0.8*sqrt(-1)],'I',Inf,'X',Inf),...
	'b',struct('aTime',0))];

D = [D, struct(...
	'DesName','One way Anova',...
	'n',	[Inf Inf 1 1],	'sF',{{'repl','group','',''}},...
	'Hform',		'I(:,2),''-'',''group''',...
	'Bform',		'I(:,3),''-'',''\mu''',...
	'nC',[0,0],'iCC',{{8,8}},'iCFI',{{1,1}},...
	'Xorth',struct('GwB',0,'HwBG',0,'CwBG',0,'CwH',0),...
	'iGXcalc',[-1,2,3],'iGMsca',[1,-9],'GM',[],...
	'iGloNorm',9,'iGC',12,...
	'M_',struct('T',[-Inf,0,0.8*sqrt(-1)],'I',Inf,'X',Inf),...
	'b',struct('aTime',0))];

D = [D, struct(...
	'DesName','Simple regression (correlation)',...
	'n',	[Inf 1 1 1],	'sF',{{'repl','','',''}},...
	'Hform',		'[]',...
	'Bform',		'I(:,2),''-'',''\mu''',...
	'nC',[1,0],'iCC',{{1,8}},'iCFI',{{1,1}},...
	'Xorth',struct('GwB',0,'HwBG',0,'CwBG',0,'CwH',0),...
	'iGXcalc',[-1,2,3],'iGMsca',[1,-9],'GM',[],...
	'iGloNorm',9,'iGC',12,...
	'M_',struct('T',[-Inf,0,0.8*sqrt(-1)],'I',Inf,'X',Inf),...
	'b',struct('aTime',0))];

D = [D, struct(...
	'DesName','AnCova',...
	'n',	[Inf Inf 1 1],	'sF',{{'repl','group','',''}},...
	'Hform',		'I(:,2),''-'',''group''',...
	'Bform',		'I(:,3),''-'',''\mu''',...
	'nC',[0,1],'iCC',{{8,1}},'iCFI',{{1,1}},...
	'Xorth',struct('GwB',0,'HwBG',0,'CwBG',0,'CwH',0),...
	'iGXcalc',[-1,2,3],'iGMsca',[1,-9],'GM',[],...
	'iGloNorm',9,'iGC',12,...
	'M_',struct('T',[-Inf,0,0.8*sqrt(-1)],'I',Inf,'X',Inf),...
	'b',struct('aTime',0))];

varargout = {D};


case 'desdefs_spm98pet'
%=======================================================================
% - Standard SPM98 PET/SPECT Design definitions...
%=======================================================================
% D = spm_spm_ui('DesDefs_SPM98PET');
% These are the PET SPM98 design definitions...

%-Single subject
%-----------------------------------------------------------------------
D = struct(...
	'DesName','Single-subject: conditions & covariates',...
	'n',	[Inf Inf 1 1],	'sF',{{'repl','condition','',''}},...
	'Hform',		'I(:,2),''-'',''cond''',...
	'Bform',		'I(:,3),''-'',''\mu''',...
	'nC',[Inf,Inf],'iCC',{{[-1,3,8],[-1,8]}},'iCFI',{{[1,3],1}},...
	'Xorth',struct('GwB',0,'HwBG',0,'CwBG',0,'CwH',0),...
	'iGXcalc',[1,2,-3],'iGMsca',[-1,9],'GM',50,...
	'iGloNorm',[1,8,9],'iGC',10,...
	'M_',struct('T',[-Inf,0.8*sqrt(-1)],'I',Inf,'X',Inf),...
	'b',struct('aTime',1));

D = [D, struct(...
	'DesName','Single-subject: covariates only',...
	'n',	[Inf 1 1 1],	'sF',{{'repl','','',''}},...
	'Hform',		'[]',...
	'Bform',		'I(:,3),''-'',''\mu''',...
	'nC',[Inf,Inf],'iCC',{{[-1,8],[-1,8]}},'iCFI',{{1,1}},...
	'Xorth',struct('GwB',0,'HwBG',0,'CwBG',0,'CwH',0),...
	'iGXcalc',[1,2,-3],'iGMsca',[-1,9],'GM',50,...
	'iGloNorm',[1,8,9],'iGC',10,...
	'M_',struct('T',[-Inf,0,0.8*sqrt(-1)],'I',Inf,'X',Inf),...
	'b',struct('aTime',1))];

%-Multi-subject
%-----------------------------------------------------------------------
D = [D, struct(...
	'DesName','Multi-subj: conditions & covariates',...
	'n',[Inf Inf Inf 1],	'sF',{{'repl','condition','subject',''}},...
	'Hform',		'I(:,2),''-'',''cond''',...
	'Bform',		'I(:,3),''-'',''subj''',...
	'nC',[Inf,Inf],'iCC',{{[1,3,4,8],[1,4,8]}},'iCFI',{{[1,3,4],[1,4]}},...
	'Xorth',struct('GwB',0,'HwBG',0,'CwBG',0,'CwH',0),...
	'iGXcalc',[1,2,-3],'iGMsca',[-4,9],'GM',50,...
	'iGloNorm',[4,8,9],'iGC',10,...
	'M_',struct('T',[-Inf,0,0.8*sqrt(-1)],'I',Inf,'X',Inf),...
	'b',struct('aTime',1))];

D = [D, struct(...
	'DesName',...
	'Multi-subj: cond x subj  interaction & covariates (FixedFX)',...
	'n',[Inf Inf Inf 1],	'sF',{{'repl','condition','subject',''}},...
	'Hform',		'I(:,[3,2]),''-'',{''subj'',''cond''}',...
	'Bform',		'I(:,3),''-'',''subj''',...
	'nC',[Inf,Inf],'iCC',{{[1,3,4,8],[1,4,8]}},'iCFI',{{[1,3,4],[1,4]}},...
	'Xorth',struct('GwB',0,'HwBG',0,'CwBG',0,'CwH',0),...
	'iGXcalc',[1,2,-3],'iGMsca',[-4,9],'GM',50,...
	'iGloNorm',[4,8,9],'iGC',10,...
	'M_',struct('T',[-Inf,0,0.8*sqrt(-1)],'I',Inf,'X',Inf),...
	'b',struct('aTime',1))];

D = [D, struct(...
	'DesName','Multi-subj: covariates only',...
	'n',[Inf 1 Inf 1],	'sF',{{'repl','','subject',''}},...
	'Hform',		'[]',...
	'Bform',		'I(:,3),''-'',''subj''',...
	'nC',[Inf,Inf],'iCC',{{[1,4,8],[1,4,8]}},'iCFI',{{[1,4],[1,4]}},...
	'Xorth',struct('GwB',0,'HwBG',0,'CwBG',0,'CwH',0),...
	'iGXcalc',[1,2,-3],'iGMsca',[-4,9],'GM',50,...
	'iGloNorm',[4,8:9],'iGC',10,...
	'M_',struct('T',[-Inf,0,0.8*sqrt(-1)],'I',Inf,'X',Inf),...
	'b',struct('aTime',0))];

%-Multi-study
%-----------------------------------------------------------------------
D = [D, struct(...
	'DesName','Multi-study: conditions & covariates',...
	'n',[Inf Inf Inf Inf],	'sF',{{'repl','condition','subject','study'}},...
	'Hform',		'I(:,[4,2]),''-'',{''stud'',''cond''}',...
	'Bform',		'I(:,[4,3]),''-'',{''stud'',''subj''}',...
	'nC',[Inf,Inf],'iCC',{{[5:8],[5,7,8]}},'iCFI',{{[1,5,6,7],[1,5,7]}},...
	'Xorth',struct('GwB',0,'HwBG',0,'CwBG',0,'CwH',0),...
	'iGXcalc',[1,2,-3],'iGMsca',[-7,9],'GM',50,...
	'iGloNorm',[7,8,9],'iGC',10,...
	'M_',struct('T',[-Inf,0,0.8*sqrt(-1)],'I',Inf,'X',Inf),...
	'b',struct('aTime',1))];

D = [D, struct(...
	'DesName','Multi-study: covariates only',...
	'n',[Inf 1 Inf Inf],	'sF',{{'repl','','subject','study'}},...
	'Hform',		'[]',...
	'Bform',		'I(:,[4,3]),''-'',{''stud'',''subj''}',...
	'nC',[Inf,Inf],'iCC',{{[5,7,8],[5,7,8]}},'iCFI',{{[1,5,7],[1,5,7]}},...
	'Xorth',struct('GwB',0,'HwBG',0,'CwBG',0,'CwH',0),...
	'iGXcalc',[1,2,-3],'iGMsca',[-7,9],'GM',50,...
	'iGloNorm',[7,8,9],'iGC',10,...
	'M_',struct('T',[-Inf,0,0.8*sqrt(-1)],'I',Inf,'X',Inf),...
	'b',struct('aTime',0))];

%-Population comparisons
%-----------------------------------------------------------------------
D = [D, struct(...
	'DesName',...
	'Population main effect: 2 cond''s, 1 scan/cond (paired t-test)',...
	'n',[1 2 Inf 1],	'sF',{{'','condition','subject',''}},...
	'Hform',		'I(:,2),''-'',''cond''',...
	'Bform',		'I(:,3),''-'',''\mu''',...
	'nC',[0,0],'iCC',{{8,8}},'iCFI',{{1,1}},...
	'Xorth',struct('GwB',0,'HwBG',0,'CwBG',0,'CwH',0),...
	'iGXcalc',[1,2,-3],'iGMsca',[-1,9],'GM',50,...
	'iGloNorm',[8,9],'iGC',10,...
	'M_',struct('T',[-Inf,0,0.8*sqrt(-1)],'I',Inf,'X',Inf),...
	'b',struct('aTime',0))];

D = [D, struct(...
	'DesName',...
	'Dodgy population main effect: >2 cond''s, 1 scan/cond',...
	'n',[1 Inf Inf 1],	'sF',{{'','condition','subject',''}},...
	'Hform',		'I(:,2),''-'',''cond''',...
	'Bform',		'I(:,3),''-'',''\mu''',...
	'nC',[0,0],'iCC',{{8,8}},'iCFI',{{1,1}},...
	'Xorth',struct('GwB',0,'HwBG',0,'CwBG',0,'CwH',0),...
	'iGXcalc',[1,2,-3],'iGMsca',[-1,9],'GM',50,...
	'iGloNorm',[8,9],'iGC',10,...
	'M_',struct('T',[-Inf,0,0.8*sqrt(-1)],'I',Inf,'X',Inf),...
	'b',struct('aTime',0))];

D = [D, struct(...
	'DesName','Compare-populations: 1 scan/subject (two sample t-test)',...
	'n',[Inf 2 1 1],	'sF',{{'subject','group','',''}},...
	'Hform',		'I(:,2),''-'',''group''',...
	'Bform',		'I(:,3),''-'',''\mu''',...
	'nC',[0,0],'iCC',{{8,8}},'iCFI',{{1,1}},...
	'Xorth',struct('GwB',0,'HwBG',0,'CwBG',0,'CwH',0),...
	'iGXcalc',[1,2,-3],'iGMsca',[-1,9],'GM',50,...
	'iGloNorm',[8,9],'iGC',10,...
	'M_',struct('T',[-Inf,0,0.8*sqrt(-1)],'I',Inf,'X',Inf),...
	'b',struct('aTime',0))];

%-The Full Monty!
%-----------------------------------------------------------------------
D = [D, struct(...
	'DesName','The Full Monty...',...
	'n',[Inf Inf Inf Inf],	'sF',{{'repl','cond','subj','group'}},...
	'Hform',		'I(:,[4,2]),''-'',{''stud'',''cond''}',...
	'Bform',		'I(:,[4,3]),''-'',{''stud'',''subj''}',...
	'nC',[Inf,Inf],'iCC',{{[1:8],[1:8]}},'iCFI',{{[1:7],[1:7]}},...
	'Xorth',struct('GwB',0,'HwBG',0,'CwBG',0,'CwH',0),...
	'iGXcalc',[1,2,3],'iGMsca',[1:7],'GM',50,...
	'iGloNorm',[1:9],'iGC',[1:11],...
	'M_',struct('T',[-Inf,0,0.8*sqrt(-1)],'I',Inf,'X',Inf),...
	'b',struct('aTime',1))];


varargout = {D};

case 'desdefs_spm96pet'
%=======================================================================
% - SPM96 PET/SPECT Design definitions...
%=======================================================================
% D = spm_spm_ui('DesDefs_SPM96PET');

%-Single subject
%-----------------------------------------------------------------------
D = struct(...
	'DesName','SPM96:Single-subject: replicated conditions',...
	'n',	[Inf Inf 1 1],	'sF',{{'repl','condition','',''}},...
	'Hform',		'I(:,2),''-'',''cond''',...
	'Bform',		'I(:,3),''-'',''\mu''',...
	'nC',[0,0],'iCC',{{8,8}},'iCFI',{{1,1}},...
	'Xorth',struct('GwB',0,'HwBG',0,'CwBG',0,'CwH',0),...
	'iGXcalc',3,'iGMsca',[1,9],'GM',50,...
	'iGloNorm',[1,8,9],'iGC',10,...
	'M_',struct('T',[0.8*sqrt(-1)],'I',0,'X',0),...
	'b',struct('aTime',0));

D = [D, struct(...
	'DesName','SPM96:Single-subject: replicated conditions & covariates',...
	'n',	[Inf Inf 1 1],	'sF',{{'repl','condition','',''}},...
	'Hform',		'I(:,2),''-'',''cond''',...
	'Bform',		'I(:,3),''-'',''\mu''',...
	'nC',[Inf,Inf],'iCC',{{1,1}},'iCFI',{{1,1}},...
	'Xorth',struct('GwB',0,'HwBG',0,'CwBG',0,'CwH',0),...
	'iGXcalc',3,'iGMsca',[1,9],'GM',50,...
	'iGloNorm',[1,8,9],'iGC',10,...
	'M_',struct('T',[0.8*sqrt(-1)],'I',0,'X',0),...
	'b',struct('aTime',0))];

D = [D, struct(...
	'DesName','SPM96:Single-subject: covariates only',...
	'n',	[Inf 1 1 1],	'sF',{{'repl','','',''}},...
	'Hform',		'[]',...
	'Bform',		'I(:,3),''-'',''\mu''',...
	'nC',[Inf,Inf],'iCC',{{1,1}},'iCFI',{{1,1}},...
	'Xorth',struct('GwB',0,'HwBG',0,'CwBG',0,'CwH',0),...
	'iGXcalc',3,'iGMsca',[1,9],'GM',50,...
	'iGloNorm',[1,8,9],'iGC',10,...
	'M_',struct('T',[0.8*sqrt(-1)],'I',0,'X',0),...
	'b',struct('aTime',0))];

%-Multi-subject
%-----------------------------------------------------------------------
D = [D, struct(...
	'DesName','SPM96:Multi-subject: different conditions',...
	'n',	[1 Inf Inf 1],	'sF',{{'','condition','subject',''}},...
	'Hform',		'I(:,2),''-'',''scancond''',...
	'Bform',		'I(:,3),''-'',''subj''',...
	'nC',[0,0],'iCC',{{8,8}},'iCFI',{{1,1}},...
	'Xorth',struct('GwB',0,'HwBG',0,'CwBG',0,'CwH',0),...
	'iGXcalc',3,'iGMsca',[1,9],'GM',50,...
	'iGloNorm',[1,4,8,9],'iGC',10,...
	'M_',struct('T',[0.8*sqrt(-1)],'I',0,'X',0),...
	'b',struct('aTime',0))];

D = [D, struct(...
	'DesName','SPM96:Multi-subject: replicated conditions',...
	'n',[Inf Inf Inf 1],	'sF',{{'repl','condition','subject',''}},...
	'Hform',		'I(:,2),''-'',''cond''',...
	'Bform',		'I(:,3),''-'',''subj''',...
	'nC',[0,0],'iCC',{{8,8}},'iCFI',{{1,1}},...
	'Xorth',struct('GwB',0,'HwBG',0,'CwBG',0,'CwH',0),...
	'iGXcalc',3,'iGMsca',[1,9],'GM',50,...
	'iGloNorm',[1,4,8,9],'iGC',10,...
	'M_',struct('T',[0.8*sqrt(-1)],'I',0,'X',0),...
	'b',struct('aTime',0))];

D = [D, struct(...
	'DesName','SPM96:Multi-subject: different conditions, & covariates',...
	'n',	[1 Inf Inf 1],	'sF',{{'','condition','subject',''}},...
	'Hform',		'I(:,2),''-'',''cond''',...
	'Bform',		'I(:,3),''-'',''subj''',...
	'nC',[Inf,Inf],'iCC',{{1,1}},'iCFI',{{[1,4],[1,4]}},...
	'Xorth',struct('GwB',0,'HwBG',0,'CwBG',0,'CwH',0),...
	'iGXcalc',3,'iGMsca',[1,9],'GM',50,...
	'iGloNorm',[1,4,8,9],'iGC',10,...
	'M_',struct('T',[0.8*sqrt(-1)],'I',0,'X',0),...
	'b',struct('aTime',0))];

D = [D, struct(...
	'DesName','SPM96:Multi-subject: replicated conditions & covariates',...
	'n',[Inf Inf Inf 1],	'sF',{{'repl','condition','subject',''}},...
	'Hform',		'I(:,2),''-'',''condition''',...
	'Bform',		'I(:,3),''-'',''subj''',...
	'nC',[Inf,Inf],'iCC',{{1,1}},'iCFI',{{[1,3,4],[1,4]}},...
	'Xorth',struct('GwB',0,'HwBG',0,'CwBG',0,'CwH',0),...
	'iGXcalc',3,'iGMsca',[1,9],'GM',50,...
	'iGloNorm',[1,4,8,9],'iGC',10,...
	'M_',struct('T',[0.8*sqrt(-1)],'I',0,'X',0),...
	'b',struct('aTime',0))];

D = [D, struct(...
	'DesName','SPM96:Multi-subject: covariates only',...
	'n',[Inf 1 Inf 1],	'sF',{{'repl','','subject',''}},...
	'Hform',		'[]',...
	'Bform',		'I(:,3),''-'',''subj''',...
	'nC',[Inf,Inf],'iCC',{{[1,4,8],[1,4,8]}},'iCFI',{{[1,4],[1,4]}},...
	'Xorth',struct('GwB',0,'HwBG',0,'CwBG',0,'CwH',0),...
	'iGXcalc',3,'iGMsca',[1,9],'GM',50,...
	'iGloNorm',[1,4,8,9],'iGC',10,...
	'M_',struct('T',[0.8*sqrt(-1)],'I',0,'X',0),...
	'b',struct('aTime',0))];

%-Multi-study
%-----------------------------------------------------------------------
D = [D, struct(...
	'DesName','SPM98:Multi-study: different conditions',...
	'n',[1 Inf Inf Inf],	'sF',{{'','cond','subj','study'}},...
	'Hform',		'I(:,[4,2]),''-'',{''study'',''cond''}',...
	'Bform',		'I(:,[4,3]),''-'',{''study'',''subj''}',...
	'nC',[0,0],'iCC',{{8,8}},'iCFI',{{1,1}},...
	'Xorth',struct('GwB',0,'HwBG',0,'CwBG',0,'CwH',0),...
	'iGXcalc',3,'iGMsca',[1,5,9],'GM',50,...
	'iGloNorm',[1,5,7,8,9],'iGC',10,...
	'M_',struct('T',[0.8*sqrt(-1)],'I',0,'X',0),...
	'b',struct('aTime',0))];

D = [D, struct(...
	'DesName','SPM98:Multi-study: replicated conditions',...
	'n',[Inf Inf Inf Inf],	'sF',{{'repl','cond','subj','study'}},...
	'Hform',		'I(:,[4,2]),''-'',{''study'',''condition''}',...
	'Bform',		'I(:,[4,3]),''-'',{''study'',''subj''}',...
	'nC',[0,0],'iCC',{{8,8}},'iCFI',{{1,1}},...
	'Xorth',struct('GwB',0,'HwBG',0,'CwBG',0,'CwH',0),...
	'iGXcalc',3,'iGMsca',[1,5,9],'GM',50,...
	'iGloNorm',[1,5,7,8,9],'iGC',10,...
	'M_',struct('T',[0.8*sqrt(-1)],'I',0,'X',0),...
	'b',struct('aTime',0))];

D = [D, struct(...
	'DesName','SPM98:Multi-study: different conditions & covariates',...
	'n',[1 Inf Inf Inf],	'sF',{{'','cond','subj','study'}},...
	'Hform',		'I(:,[4,2]),''-'',{''study'',''cond''}',...
	'Bform',		'I(:,[4,3]),''-'',{''study'',''subj''}',...
	'nC',[Inf,Inf],'iCC',{{1,1}},'iCFI',{{[1,5,6,7],[1,5,7]}},...
	'Xorth',struct('GwB',0,'HwBG',0,'CwBG',0,'CwH',0),...
	'iGXcalc',3,'iGMsca',[1,5,9],'GM',50,...
	'iGloNorm',[1,5,7,8,9],'iGC',10,...
	'M_',struct('T',[0.8*sqrt(-1)],'I',0,'X',0),...
	'b',struct('aTime',0))];

D = [D, struct(...
	'DesName','SPM98:Multi-study: replicated conditions & covariates',...
	'n',[Inf Inf Inf Inf],	'sF',{{'','cond','subj','study'}},...
	'Hform',		'I(:,[4,2]),''-'',{''study'',''condition''}',...
	'Bform',		'I(:,[4,3]),''-'',{''study'',''subj''}',...
	'nC',[Inf,Inf],'iCC',{{1,1}},'iCFI',{{[1,5,6,7],[1,5,7]}},...
	'Xorth',struct('GwB',0,'HwBG',0,'CwBG',0,'CwH',0),...
	'iGXcalc',3,'iGMsca',[1,5,9],'GM',50,...
	'iGloNorm',[1,5,7,8,9],'iGC',10,...
	'M_',struct('T',[0.8*sqrt(-1)],'I',0,'X',0),...
	'b',struct('aTime',0))];

D = [D, struct(...
	'DesName','SPM98:Multi-study: covariates only',...
	'n',[Inf 1 Inf Inf],	'sF',{{'repl','','subj','study'}},...
	'Hform',		'[]',...
	'Bform',		'I(:,[4,3]),''-'',{''study'',''subj''}',...
	'nC',[Inf,Inf],'iCC',{{1,1}},'iCFI',{{[1,5,7],[1,5,7]}},...
	'Xorth',struct('GwB',0,'HwBG',0,'CwBG',0,'CwH',0),...
	'iGXcalc',3,'iGMsca',[1,5,9],'GM',50,...
	'iGloNorm',[1,5,7,8,9],'iGC',10,...
	'M_',struct('T',[0.8*sqrt(-1)],'I',0,'X',0),...
	'b',struct('aTime',0))];

%-Group comparisons
%-----------------------------------------------------------------------
D = [D, struct(...
	'DesName','Compare-groups: 1 scan per subject',...
	'n',[Inf Inf 1 1],	'sF',{{'subject','group','',''}},...
	'Hform',		'I(:,2),''-'',''group''',...
	'Bform',		'[]',...
	'nC',[0,0],'iCC',{{8,8}},'iCFI',{{1,1}},...
	'Xorth',struct('GwB',0,'HwBG',0,'CwBG',0,'CwH',0),...
	'iGXcalc',3,'iGMsca',[1,9],'GM',50,...
	'iGloNorm',[1,8,9],'iGC',10,...
	'M_',struct('T',[0.8*sqrt(-1)],'I',0,'X',0),...
	'b',struct('aTime',0))];

varargout = {D};


otherwise
%=======================================================================
% - U N K N O W N   A C T I O N
%=======================================================================
warning(['Illegal Action string: ',Action])


%=======================================================================
% - E N D
%=======================================================================
end




%=======================================================================
%- S U B - F U N C T I O N S
%=======================================================================

function str = sf_estrrep(str,srstr)
%=======================================================================
for i = 1:size(srstr,1)
	str = strrep(str,srstr{i,1},srstr{i,2});
end

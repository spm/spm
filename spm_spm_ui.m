function spm_spm_ui
%_______________________________________________________________________
% %W% Andrew Holmes, Karl Friston %E%
SCCSid  = '%I%';

%=======================================================================
% - S E T U P
%=======================================================================
SPMid = spm('FnBanner',mfilename,SCCSid);
[Finter,Fgraph,CmdLine] = spm('FnUIsetup','Setup statistical analysis',0);
spm_help('!ContextHelp',[mfilename,'.m'])

%-Design parameters
%=======================================================================
%-Note: struct expands cell array values to give multiple records:
%       => must embed cell arrays within another cell array!
%-Negative indices for iCC indicates default (first taken)
D = struct(...
	'DesName',	'Single-subject: replication of conditions',...
	'n',	[Inf 2 1 1],	'sF',{{'repl','condition','subject','study'}},...
	'Hform',		'i2,''+0m'',''condition''',...
	'Bform',		'i3,''-'',''constant''',...
	'nC',[Inf,Inf],'iCC',{{[1:8],1}},'iCFI',{{[1,-2,3:7],1}},...
	'Xorth',struct('GwB',0,'HwBG',0,'CwBG',0,'CwH',0),...
	'iGXcalc',[1:2,-3],'iGMsca',[1:7,8,-9],'GM',[],...
	'iGloNorm',[1:9],'iGC',[1:9,-10,-11],...
	'b',struct('aTime',1));
iDD = length(D); %-Above design is the default design definition


%-Option definitions
%-----------------------------------------------------------------------
%-Generic factor names
sF = {'sF1','sF2','sF3','sF4'};

%-Covariate by factor interaction options
sCFI = {'<none>';...							%-1
	'with sF1';'with sF2';'with sF3';'with sF4';...			%-2:5
	'with sF2 (by sF4)';'with sF3 (by sF4)'};			%-6,7

%-DesMtx argument components for covariate by factor interaction options
% (Used for CFI's Covariate Centering (CC), GMscale & Global normalisation)
CFIforms = {	'[]',		'C',	'{}';...			%-1
		'i1',		'FxC',	'{D.sF{1}}';...			%-2
		'i2',		'FxC',	'{D.sF{2}}';...			%-3
		'i3',		'FxC',	'{D.sF{3}}';...			%-4
		'i4',		'FxC',	'{D.sF{4}}';...			%-5
		'[i4,i2]',	'FxC',	'{D.sF{4},D.sF{2}}';...		%-6
		'[i4,i3]',	'FxC',	'{D.sF{4},D.sF{3}}'	};	%-7

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
CCforms = {'ones(size(i1))',CFIforms{2:end,1},''}';


%-Global normalization options (options 1-7 match CFIforms)       (GloNorm)
sGloNorm = {	'AnCova';...						%-1
		'AnCova by sF1';...					%-2
		'AnCova by sF2';...					%-3
		'AnCova by sF3';...					%-4
		'AnCova by sF4';...					%-5
		'AnCova by sF2 (by sF4)';...				%-6
		'AnCova by sF3 (by sF4)';...				%-7
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
dGM = 50;			%-Default GM (if isempty(D.GM))


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
D = D(spm_input('Select design type...',1,'m',{D.DesName}',[],iDD));

%-Set factor names for this design
%-----------------------------------------------------------------------
sCC       = sf_estrrep(sCC,[sF',D.sF']);
sCFI      = sf_estrrep(sCFI,[sF',D.sF']);
sGloNorm  = sf_estrrep(sGloNorm,[sF',D.sF']);
sGMsca    = sf_estrrep(sGMsca,[sF',D.sF']);

%-Get filenames, accrue study, subject, condition & replication indicies
%-----------------------------------------------------------------------
i4 = [];		% factor 4 index (usually study)
i3 = [];		% factor 3 index (usually subject), per f4
i2 = [];		% factor 2 index (usually condition), per f3/f4
i1 = [];		% factor 1 index (usually replication), per f2/f3/f4
P  = {};		% cell array of string filenames

if isinf(D.n(4)), n4 = spm_input(['#',D.sF{4},'.s'],'+1','w',[],1);
	else, n4 = D.n(4); end
bL4 = n4>1;

ti2 = '';
GUIpos = spm_input('!NextPos');
for j4  = 1:n4
    spm_input('!SetNextPos',GUIpos);
    sF4P=''; if bL4, sF4P=[D.sF{4},'.',int2str(j4),': ']; end
    if isinf(D.n(3)), n3=spm_input([sF4P,'#',D.sF{3},'.s'],'+1','w',[],1);
	    else, n3 = D.n(3); end
    bL3 = n3>1;
    
    if D.b.aTime & D.n(2)>1
	disp('NB:selecting in time order - manually specify conditions')%-**
	%-NB: This means f2 levels might not be 1:n2
	GUIpos2 = spm_input('!NextPos');
	for j3 = 1:n3
	    sF3P=''; if bL3, sF3P=[D.sF{3},'.',int2str(j3),': ']; end
	    str = [sF4P,sF3P];
	    tP = spm_get(D.n(2)*D.n(1),'.img',{[str,'select scans...']});
	    n21 = length(tP);
	    ti2 = spm_input([str,'iCond'],GUIpos2,'c',ti2',n21,D.n(2));
	    %-Work out i1 & check
	    [tl2,null,j] = unique(ti2);
	    tn1 = zeros(size(tl2)); ti1 = zeros(size(ti2));
	    for i=1:length(tl2)
		    tn1(i)=sum(j==i); ti1(ti2==tl2(i))=1:tn1(i); end
	    if isfinite(D.n(1)) & any(tn1~=D.n(1))
		%-#i1 levels mismatches specification in D.n(1)
		error(sprintf('#%s not %d as pre-specified',D.sF{1},D.n(1)))
	    end
	    P   = [P;tP];
	    i4 = [i4, j4*ones(n21,1)];
	    i3 = [i3, j3*ones(n21,1)];
	    i2 = [i2, ti2];
	    i1 = [i1, ti1];
	end

    else

	if isinf(D.n(2))
	    n2 = spm_input([sF4P,'#',D.sF{2},'.s'],'+1','w',[],1);
	else
	    n2 = D.n(2);
	end
	bL2 = n2>1;

	if n2==1 & D.n(1)==1 %-single scan per f3 (subj)
	    disp('NB:single scan per f3')%-**
	    str = [sF4P,'select scans, ',D.sF{3},' 1-',int2str(n3)];
	    P     = [P;spm_get(n3,'.img',{str})];
	    i4 = [i4, j4*ones(n3,1)];
	    i3 = [i3, [1:n3]'];
	    i2 = [i2, ones(n3,1)];
	    i1 = [i1, ones(n3,1)];
	else
	    %-multi scan per f3 (subj) case
	    disp('NB:multi scan per f3')%-**
	    for j3 = 1:n3
		sF3P=''; if bL3, sF3P=[D.sF{3},'.',int2str(j3),': ']; end
		if D.n(1)==1
			%-No f1 (repl) within f2 (cond)
			disp('NB:no f1 within f2')%-**
			str = [sF4P,sF3P,'select scans: ',D.sF{2},...
				 ' 1-',int2str(n2)];
			P = [P;spm_get(n2,'.img',{str})];
			i4 = [i4, j4*ones(n2,1)];
			i3 = [i3, j3*ones(n2,1)];
			i2 = [i2, [1:n2]'];
			i1 = [i1, ones(n2,1)];
		else
		    %-multi f1 (repl) within f2 (cond)
		    disp('NB:f1 within f2')%-**
		    for j2 = 1:n2
			sF2P='';
			if bL2, sF2P=[D.sF{2},'.',int2str(j2),': ']; end
			str = [sF4P,sF3P,sF2P,' select scans...'];
			tP  = spm_get(D.n(1),'.img',{str});
			n1 = size(tP,1);
			P   = [P;tP];
			i4 = [i4, j4*ones(n1,1)];
			i3 = [i3, j3*ones(n1,1)];
			i2 = [i2, j2*ones(n1,1)];
			i1 = [i1, [1:n1]'];
		    end                         % (for j2)
		end                             % (if D.n(1)==1)
	    end                                 % (for j3)
	end                                     % (if  n2==1 &...)
    end                                         % (if D.b.aTime & D.n(2)>1)
end                                             % (for j4)
clear n1 n2 n3 n4 bL1 bL2 bL3 bL4
nScan = length(P);			%-#observations
I     = [i1,i2,i3,i4];			%-Factor indices

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
       	cname  = spm_input([str,' name?'],'+1','s',...
       		[dcname{i},'^{',num2str(nc(i)),'}']);
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
        %---------------------------------------------------------------
        if length(D.iCC{i})>1
        	%-User has a choice of centering options
		%-Only offer factor specific for appropriate factor combinations
		iCC = intersect(abs(D.iCC{i}),find([1,bFI,1]) );
        	%-Default is max -ve option in D, overridden by iCFI if CFI
		if iCFI==1, dCC=-D.iCC{i}(D.iCC{i}<0); else, dCC=iCFI; end
		dCC = max([1,intersect(iCC,dCC)]);
		iCC = spm_input([str,': centre?'],'+1','m',...
			sCC(iCC),iCC,find(iCC==dCC));
        else
        	iCC = abs(D.iCC{i});		%-AutoSelect default option
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
else                                    %-Set value of GM
    if ~isempty(D.GM) & D.GM~=0 %-D.GM contains GM value to use
        GM = D.GM;
    else                                %-Ask user for GM
        %-Default response is dGM
        if iGloNorm==8
            str='PropSca global mean to';
        else
            str=[strrep(sGMsca{iGMsca},'scaling of','scale'),' to'];
        end
        GM = spm_input(str,'+1','r',dGM,1);
    end
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
	gc  = 0;
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
	%-Omit '(as implied by AnCova)' if same as GM
	if iGloNorm==iGMsca, iGC = setdiff(iGC,10); end

	%-Set defaults (if any), & get answer
	%---------------------------------------------------------------
	dGC = max([0,intersect(iGC,-D.iGC(D.iGC<0))]);
	str = 'Centre global covariate';
	if iGMsca<8, str = [str,' (after grand mean scaling)']; end
	iGC = spm_input(str,'+1','m',sCC(iGC),iGC,find(iGC==dGC));

	%-If "as implied by AnCova" then set iGC accordingly
	%---------------------------------------------------------------
	if iGC==10, iGC = iGloNorm; end

	%-If 'user specified' then get value
	%---------------------------------------------------------------
	if iGC==9
		gc     = spm_input('Centre globals around','+0','r',dGM,1)
		sCC{9} = sprintf('%s of %g',sCC{iGC},gc);
	else
		gc  = 0;
	end
end


%-Global calculation                                            (GXcalc)
%-----------------------------------------------------------------------
iGXcalc = abs(D.iGXcalc);
%-Only offer "omit" option if not doing any GloNorm or GMsca
if ~(iGloNorm==9 & iGMsca==9)
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


%-Thresholds & masks defining voxels to analyse                  (THRESH)
%-----------------------------------------------------------------------
%-Threshold mask as proportion of global
if iGXcalc>1, str = ' (times global)'; else, str = ''; end
M_THp = spm_input(['Analysis threshold',str],'+1','br1','None',-Inf,0.8);
if isfinite(M_THp)
	xsM.Analysis_threshold = sprintf('images thresholded at %6f%s',...
		M_THp,str);
else
	xsM.Analysis_threshold = 'None (-Inf)';
end

%-Implicit masking: Ignore zero voxels in low data-types?
% (Implicit mask is NaN in higher data-types.)
type = getfield(spm_vol(P{1}),'dim')*[0,0,0,1]';
if ~spm_type(type,'nanrep')
	M_I = spm_input('Implicit mask (ignore zero''s)?','+1','y/n',[1,0],1);
	if M_I
		xsM.Implicit_masking = 'Yes: zero''s treated as missing';
	else
		xsm.Implicit_masking = 'No';
	end
else
	M_I = 1;
	xsM.Implicit_masking = 'Yes: NaN''s treated as missing';
end

%-Explicit mask images (map them later...)
M_X = spm_input('explicit mask images?','+1','y/n',[1,0],2);
if M_X, M_P = spm_get(Inf,'*.img',{'select mask images'}); else, M_P = {}; end


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
	elseif isfinite(M_THp)
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
if iGXcalc>1, M_TH = (rg.*gSF)*M_THp; else, M_TH = M_THp; end
if ~isempty(M_P)
	VM = spm_vol(char(M_P));
	xsM.Explicit_masking = [{'Yes: mask images :'};{VM.fname}'];
else
	VM={};
	xsM.Explicit_masking = 'No';
end
xM = struct('THp',M_THp, 'TH',M_TH, 'I',M_I, 'VM',{VM}, 'xs',xsM);



%-Design matrix canonicalisation
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
save


%=======================================================================
% - D I S P L A Y   A N A L Y S I S   P A R A M E T E R S
%=======================================================================
Fgraph = spm_figure('GetWin','Graphics');
spm_figure('Clear',Fgraph), figure(Fgraph)
FS = spm_figure('FontSizes');


%-Files, indices & covariates... & masking options
%-----------------------------------------------------------------------
spm_DesRep('Files&Vars',{VY.fname}',xX.I,xC,D.sF,xM.xs)
spm_print


%-Design matrix & design descriptions
%-----------------------------------------------------------------------
spm_DesRep('DesMtx',xX.X,xX.Xnames,{VY.fname}',xsDes)
spm_print


%-Covariates
%-----------------------------------------------------------------------
spm_DesRep('Covs',xC,xX.X,xX.Xnames)
spm_print


%=======================================================================
% - A N A L Y S I S   P R O P E R
%=======================================================================
%-This is PET data so sigma = 0 (i.e. independent observations) 
% RT is undefined
%-----------------------------------------------------------------------
%spm_spm(VY,H,C,B,G,CONTRAST,ORIGIN,THRESH*GX,HCBGnames,P,0,[])



%=======================================================================
% - E N D
%=======================================================================
spm('FigName','Stats: done',Finter,CmdLine); spm('Pointer','Arrow')
fprintf('\n\n')







%=======================================================================
%- S U B - F U N C T I O N S
%=======================================================================

function str = sf_estrrep(str,srstr)
%=======================================================================
for i = 1:size(srstr,1)
	str = strrep(str,srstr{i,1},srstr{i,2});
end

%****function spm_spm_ui
%_______________________________________________________________________
% %E% Andrew Holmes, Karl Friston %W%

%=======================================================================
% - S E T U P
%=======================================================================
SCCSid = '%I%';
spm('FnBanner',mfilename,SCCSid);
[Finter,Fgraph,CmdLine] = spm('FnUIsetup','Statistical analysis',1);
spm_help('!ContextHelp',[mfilename,'.m'])

%-Design parameters
%=======================================================================
%-Note: struct expands cell array values to give multiple records:
%       => must embed cell arrays within another cell array!
%-Negative indices for iCC indicates default (first taken)
D = struct(...
	'DesName',	'Single-subject: replication of conditions',...
	'n',	[Inf 2 1 1],...
	'sF',	{{'repl','condition','subject','study'}},...
	'Hform',		'i2,''+0m'',''cond''',...
	'Bform',		'i3,''-'',''const''',...
	'nC',[Inf,Inf],...
	'iCC',{{[1:8],1}},...
	'iCFI',{{[1,-2,3:7],1}},...
	'iGloNorm',[1:9],	'iGMsca',[1:7,8,-9],	'GM',[],...
	'b',struct('aTime',1));
iDD = length(D); %-Above design is the default design definition


%-Option definitions
%-----------------------------------------------------------------------
%-Generic factor names
sF = {'sF1','sF2','sF3','sF4'}

%-Covariate by factor interaction options
sCFI = {'<none>',...							%-1
	'with sF1','with sF2','with sF3','with sF4',...			%-2:5
	'with sF2 (by sF4)','with sF3 (by sF4)'};			%-6,7

%-DesMtx argument components for covariate by factor interaction options
% (Used for CFI's Covariate Centering (CC), GMscale & Global normalisation)
CFIforms = {	'[]',		'C',	'{}';...			%-1
		'i1',		'FxC',	'{D.sF{1}}';...			%-2
		'i2',		'FxC',	'{D.sF{2}}';...			%-3
		'i3',		'FxC',	'{D.sF{3}}';...			%-4
		'i4',		'FxC',	'{D.sF{4}}';...			%-5
		'[i4,i2]',	'FxC',	'{D.sF{4},D.sF{2}}';...		%-6
		'[i4,i3]',	'FxC',	'{D.sF{4},D.sF{3}}'	};	%-7

%-Centre (mean correction) options for covariates                      (CC)
sCC = {		'around overall mean',...				%-1
		'around sF1 means',...					%-2
		'around sF2 means',...					%-3
		'around sF3 means',...					%-4
		'around sF4 means',...					%-5
		'around sF2 (within sF4) means',...			%-6
		'around sF3 (within sF4) means',...			%-7
		'<no centering>'};					%-8
%-DesMtx I forms for centering options
CCforms = {'ones(size(i1))',CFIforms{2:end,1},''};

%-Global normalization options                                    (GloNorm)
sGloNorm = {	'AnCova',...						%-1
		'AnCova by sF1',...					%-2
		'AnCova by sF2',...					%-3
		'AnCova by sF3',...					%-4
		'AnCova by sF4',...					%-5
		'AnCova by sF2 (by sF4)',...				%-6
		'AnCova by sF3 (by sF4)',...				%-7
		'proportional scaling',...				%-8
		'<no global normalisation>'};				%-9

%-Grand mean scaling options                                        (GMsca)
sGMsca = {	'scaling of overall grand mean',...			%-1
		'scaling of sF1 grand means',...			%-2
		'scaling of sF2 grand means',...			%-3
		'scaling of sF3 grand means',...			%-4
		'scaling of sF4 grand means',...			%-5
		'scaling of sF2 (within sF4) grand means',...		%-6
		'scaling of sF3 (within sF4) grand means',...		%-7
		'(implicit in PropSca global normalisation)',...	%-8
		'<no grand Mean scaling>'	};			%-9
%-NB: Grand mean scaling by subject is redundent for proportional scaling
dGM = 50;			%-Default GM (if isempty(D.GM))


%-Centering options for AnCova GloNorm                                 (GC)
sGC  = {sCC{1:7},...							%-1:7
	'around specified value','<no centering>',...			%-8,9
	'(as implied by AnCova)','GM',...				%-10,11
	'(redundant: not doing AnCova)'};				%-12
DiGC = [1:9,-10,-11];

%-Global calculation options                                       (GXcalc)
sGXcalc  = {	'mean voxel value (within per image fullmean/8 mask)',...%-1
		'user specified',...					%-2
		'omit'};						%-3
DiGXcalc = [-1,2:3];


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
sGC       = sf_estrrep(sGC,[sF',D.sF']);

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

%-Additional design parameters
%-----------------------------------------------------------------------
nScan = length(P);			%-#observations
I   = [i1,i2,i3,i4];			%-Factor indices
bL  = any(diff([i1,i2,i3,i4],1),1); 	%-Multiple factor levels?
	% NB: bL(2) might be thrown by user specified f1 levels
	%     (D.b.aTime & D.n(2)>1) - assumme user is consistent?
bFI = [bL(1),bL(2:3)&~bL(4),bL(4),bL([2,3])&bL(4)];
	%-Allowable interactions for covariates
	%-Only offer interactions with multi-level factors, and
	% don't offer by F2/F3 if bL(4)!

%-Build Condition (H) and Block (B) partitions
%=======================================================================
eval(['[H,Hnames] = spm_DesMtx(',D.Hform,');'])
if rank(H)==nScan, error('unestimable condition effects'), end
eval(['[B,Bnames] = spm_DesMtx(',D.Bform,');'])
if rank(B)==nScan, error('unestimable block effects'), end


%-Covariate partition(s): interest & nuisance (excluding global)
%=======================================================================
nC = D.nC;			%-Default #covariates
C  = {[],[]}; Cnames = {{},{}};	%-Covariate DesMtx partitions & names
rC = [];			%-Struct array to hold raw covariates
				% Fields:{'c','cname','iCC','iCFI','type','cols'}

dcname = {'CovInt','NusCov'};	%-Root names for covariates
dstr   = {'covariate','nuisance variable'};

GUIpos = spm_input('!NextPos');
nc     = [0,0];
for i=1:2			% 1:covariates of interest, 2:nuisance variables
    if isinf(nC(i)), nC(i)=spm_input(['# ',dstr{i},'s'],GUIpos,'w1'); end

    while nc(i) < nC(i)

	if nC(i)==1, str=dstr{i}; else, str=sprintf('%s %d',dstr{i},nc(i)+1); end
        c = spm_input(str,GUIpos,'r',[],[nScan,Inf]);
        if any(isnan(c(:))), break, end		%-NaN is dummy value to exit
	nc(i)  = nc(i)+1;			%-#Covariates (so far)
       	cname  = [dcname{i},'^{',num2str(nc(i)),'}'];
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
        	iCFI = spm_input(['Covariate',int2str(nc(i)),...
			': interaction'],GUIpos+1,'m',...
			sCFI(iCFI),iCFI,find(iCFI==dCFI));
	    else
		iCFI = abs(D.iCFI{i});
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
		iCC = spm_input('Centre covariate',GUIpos+1,'m',...
			sCC(iCC),iCC,find(iCC==dCC));
        else
        	iCC = abs(D.iCC{i});
        end
	%-Centre within factor levels as appropriate
        if any(iCC==[1:7]), c = c - spm_meanby(c,eval(CCforms{iCC})); end

        %-Do Interaction (if single covariate vector entered)?
        %---------------------------------------------------------------
        if iCFI>1				%-(NB:iCFI=1 if size(c,2)>1)
       		tI        = [eval(CFIforms{iCFI,1}),c];
		tConst    = CFIforms{iCFI,2};
		tFnames   = [eval(CFIforms{iCFI,3}),{cname}];
		[c,cname] = spm_DesMtx(tI,tConst,tFnames);
	elseif size(c,2)>1			%-Design matrix block
		[null,cname] = spm_DesMtx(c,'X',cname);
	end

	%-Store raw covariate details in rC struct for reference
	%-Pack c into appropriate DesMtx partition
        %---------------------------------------------------------------
	tmp       = struct(	'c',rc,		'cname',rcname,...
				'iCC',iCC,	'iCFI',iCFI,...
				'type',i,...
				'cols',[1:size(c,2)]+size(C{i},2)	);
	if isempty(rC), rC=tmp; else, rC=[rC,tmp]; end
	C{i}      = [C{i},c];
	Cnames{i} = [Cnames{i}; cname];

    end	% (while)

end % (for)
clear i %-****
clear tI tConst tFnames

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
sGloNorm = sGloNorm{iGloNorm};


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
            str='GM: PropSca global mean to';
        else
            str=['GM: ',strrep(sGMsca{iGMsca},'scaling of','scale'),' to'];
        end
        GM = spm_input(str,'+0','e',dGM);
    end
    if GM==0, iGMsca=9; end             %-If GM is zero then don't GMsca!
end
sGMsca = sGMsca{iGMsca};


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
%-Setting DiGC=[-10,-11] gives the standard choices above

%-If not doing AnCova then GC is irrelevant
if ~any(iGloNorm==[1:7])
	iGC = 12;
else
	%-Annotate options 10 & 11 with specific details
	%---------------------------------------------------------------
	%-Tag '(as implied by AnCova)' with actual AnCova situation
	sGC{10} = [sCC{iGloNorm},' (as implied by ',sGloNorm{iGloNorm},')'];
	%-Tag 'GM' case with actual GM & GMsca case
	sGC{11} = sprintf('around GM=%g (i.e. %s after grand mean scaling)',...
		GM,strrep(sCC{iGMsca},'around ',''));

	%-Constuct vector of allowable iGC
	%---------------------------------------------------------------
	%-Weed out redundent factor combinations from pre-set allowable options
	iGC = intersect(abs(DiGC),find([1,bFI,1,1,1,1]));
	%-Omit 'GM' option if didn't GMsca (iGMsca~=8 'cos doing AnCova)
	if any(iGMsca==[8,9]), iGC = setdiff(iGC,11); end
	%-Omit '(as implied by AnCova)' if same as GM
	if iGloNorm==iGMsca, iGC = setdiff(iGC,10); end

	%-Set defaults (if any), & get answer
	%---------------------------------------------------------------
	dGC = max([0,intersect(iGC,-DiGC(DiGC<0))]);
	str = 'GC: Centre global covariate';
	if iGMsca<8, str = [str,' (after grand mean scaling)']; end
	iGC = spm_input(str,'+1','m',sGC(iGC),iGC,find(iGC==dGC));

	%-If 'user specified' then get value
	%---------------------------------------------------------------
	if iGC=8, GC = spm_input('Centre globals around','+0','r',dGM,[1,1]);
		else, GC = 0; end
end
sGC = sGC{iGC};


%-Global calculation                                            (GXcalc)
%-----------------------------------------------------------------------
%-Only offer "omit" option if not doing any GloNorm or GMsca
iGXcalc = intersect(abs(DiGXcalc),find([1,1,iGloNorm==9 & iGMsca==9]));
if isempty(iGXcalc)
	error('no GXcalc options')
elseif length(iGXcalc)>1
	%-User choice of global calculation options, default is negative
	dGXcalc = max([1,intersect(iGXcalc,-DiGXcalc(DiGXcalc<0))]);
	iGXcalc = spm_input('Global calculation','+1','m',...
	    	sGXcalc(iGXcalc),iGXcalc,find(iGXcalc==dGXcalc));
else
	iGXcalc = abs(DiGXcalc);
end
sGXcalc = sGXcalc{iGXcalc};

if iGXcalc==2				%-Get user specified globals
	GX = spm_input('globals','+0','r',[],[nScan,1]);
end
sGXcalc = sGXcalc{iGXcalc};


%-Threshold defining voxels to analyse                          (THRESH)
%-----------------------------------------------------------------------
vMask = spm_input('Analysis (grey&white) threshold','+1','e',0.8);
%-**** Allow -Inf, default, or mask image? Default -Inf for -ve data?



%=======================================================================
% - C O N F I G U R E   D E S I G N
%=======================================================================
spm('FigName','Stats: configuring',Finter,CmdLine); fprintf('\tconfiguring: ')
spm('Pointer','Watch');


%-Images & image info
%=======================================================================

%-File handles
%-----------------------------------------------------------------------
V = spm_vol(char(P));

%-Check for consistency of image dimensions and orientation / voxel size
%-----------------------------------------------------------------------
if any(any(diff(cat(1,V.dim),1,1),1)&[1,1,1,0])	%NB: Bombs for single image
	error('images do not all have the same dimensions'), end
if any(any(any(diff(cat(3,V.mat),1,3),3)))
	error('images do not all have same orientation & voxel size'), end

%-Work out required Analyze header info from handles
%-----------------------------------------------------------------------
DIM    = V(1).dim(1:3);
VOX    = sqrt(sum(V(1).mat(1:3,1:3).^2));
ORIGIN = (V(1).mat\[0 0 0 1]')';
ORIGIN = round(ORIGIN(1:3));


%-Global values, scaling and global normalisation
%=======================================================================
%-Compute global values
%-----------------------------------------------------------------------
fprintf('(globals)')
switch iGXcalc, case 1
	%-Compute as mean voxel value (within per image fullmean/8 mask)
	GX = zeros(nScan,1);
	for i = 1:nScan, GX(i) = spm_global(V(i)); end
case 2
	%-User specified globals
case 3
	%-Don't compute => no GMsca (iGMsca==9) or GloNorm (iGloNorm==9)
	GX = [];
otherwise
	error('illegal iGXcalc')
end
fprintf('\b - done)\n')


%-Scaling: compute global scaling factors gSF required to implement proportional
% scaling global normalisation (PropSca) or grand mean scaling (GMsca),
% as specified by iGMsca (& iGloNorm)
%-----------------------------------------------------------------------
rGX = GX;
switch (iGMsca), case 8
	%-Proportional scaling global normalisation
	if iGloNorm~=8, error('iGloNorm-iGMsca(8) mismatch for PropSca'), end
	gSF    = GM./GX;
	GX     = GM*ones(nScan,1);
case {1,2,3,4,5,6,7}
	%-Grand mean scaling according to iGMsca
	gSF    = GM./spm_meanby(GX,eval(CCforms{iGMsca}));
	GX     = GX.*gSF;
	sGMsca = sprintf('%s to %g',sGMsca,GM);
case 9
	%-No grand mean scaling
	gSF    = ones(nScan,1);
otherwise
	error('illegal iGMsca')
end


%-AnCova: Construct global nuisance covariates partition (if AnCova)
%-----------------------------------------------------------------------
%-Save info in rC

%-===============================-*@*-=================================-


%-Construct Global part of covariates of no interest partition.
%-Centre global means if included in AnCova models, by mean correction.
%-----------------------------------------------------------------------
%-Save scaled globals for printing later on
Gc      = [Gc,GX];
Gcnames = strvcat(Gcnames,'Global');

if iGloNorm == 1				%-No global adjustment
%-----------------------------------------------------------------------

elseif iGloNorm == 2				%-Proportional scaling
%-----------------------------------------------------------------------
    if (GM ~= 0)
	V(7,:) = GM*V(7,:)./GX'; GX = ones(size(GX))*GM;
    else
	V(7,:) = V(7,:)./GX'; GX = ones(size(GX));
    end

elseif iGloNorm == 3				%-AnCova
%-----------------------------------------------------------------------
    G      = [G,(GX - mean(GX))];
    Gnames = strvcat(Gnames,'Global');

elseif iGloNorm == 4				%-AnCova by subject
%-----------------------------------------------------------------------
    [GL,GLnames] = spm_DesMtx([iSUBJ',GX-mean(GX)],'FxC',['SUBJ  ';'Global']);
    G      = [G,GL];
    Gnames = strvcat(Gnames,GLnames);

elseif iGloNorm == 5				%-AnCova by study
%-----------------------------------------------------------------------
    [GL,GLnames] = spm_DesMtx([i4',GX-mean(GX)],'FxC',['Stud  ';'Global']);
    G      = [G,GL];
    Gnames = strvcat(Gnames,GLnames);
else

    fprintf('%cError: invalid iGloNorm option\n',7)

end % (if iGloNorm)

%-Construct full design matrix and name matrices for display
%-----------------------------------------------------------------------
[nHCBG,HCBGnames] = spm_DesMtx('sca',H,Hnames,C,Cnames,B,Bnames,G,Gnames);


%-Ensure validity of contrast of condition effects, zero pad
%-----------------------------------------------------------------------
if ~isempty(CONTRAST)
	if ~isempty(H)
		%-Ensure contrasts of cond effects sum to zero within study
		if ~bMStud | bBetGrp
			tmp = ones(size(H,2),1);
		else
			tmp = zeros(size(H,2),1);
			tmp(cumsum([1,n2(1:nStud-1)])) = ones(nStud,1);
			tmp = cumsum(tmp);
		end
		d        = 1:size(H,2);
		CONTRAST(:,d) = CONTRAST(:,d)-spm_meanby(CONTRAST(:,d)',tmp)';
	end
	%-Remove zero contrasts
	CONTRAST(find(all(CONTRAST'==0)),:)=[];
	
	%-zero pad for B & G partitions
	CONTRAST = [CONTRAST, zeros(size(CONTRAST,1),size([B G],2))];
end

%-Display analysis parameters
%=======================================================================

%-Compute common path components - all paths will begin with '/'
%-----------------------------------------------------------------------
d     = max(find(P(1,1:min(find(~all(P == ones(q,1)*P(1,:))))-1)=='/')) - 1;
CPath = P(1,1:d);
Q     = P(:,[(d + 1):size(P,2)]);


%-Display
%-----------------------------------------------------------------------
figure(Fgraph); spm_clf; axis off
text(0.30,1.02,'Statistical analysis','Fontsize',16,'Fontweight','Bold');
text(-0.10,0.85,'Scan Index','Rotation',90)
if bMStud, text(-0.05,0.85,sStud,        'Rotation',90); end
if bMSubj, text(+0.00,0.85,'Subject',    'Rotation',90); end
if bMCond, text(+0.05,0.85,'Condition',  'Rotation',90); end
if bMRepl, text(+0.10,0.85,'Replication','Rotation',90); end
x0    = 0.15; y0 = 0.83;
dx    = 0.10; dy = 0.02;
x     = x0;
for i = 1:size(Cc,2)
	text(x + 0.02,0.85,Ccnames(i,:),'Rotation',90);
	x = x + dx; end
for i = 1:size(Gc,2)
	text(x + 0.02,0.85,Gcnames(i,:),'Rotation',90);
	x = x + dx; end
text(x,0.92,'Base directory:','FontSize',10,'Fontweight','Bold');
text(x,0.90,CPath,'FontSize',10);
text(x,0.87,'Filename Tails');
y     = y0;
for i = 1:q
	text(-0.12,y,sprintf('%02d :',i));
	if bMStud, text(-0.06,y,sprintf('%2d',i4(i))); end
	if bMSubj, text(-0.01,y,sprintf('%2d',i3(i))); end
	if bMCond, text(+0.04,y,sprintf('%2d',i2(i))); end
	if bMRepl, text(+0.09,y,sprintf('%2d',i1(i))); end
	x     = x0;
	for j = 1:size(Cc,2)
		text(x,y,sprintf('%-8.6g',Cc(i,j)),'FontSize',10)
		x = x + dx; end
	for j = 1:size(Gc,2)
		text(x,y,sprintf('%-8.6g',Gc(i,j)),'FontSize',10)
		x = x + dx; end
	text(x,y,Q(i,:),'FontSize',10);
	y     = y - dy;
	if y < 0;
		spm_print
		spm_clf(Fgraph); axis off
		y = y0;
		text(0.16,1.02,['Statistical analysis (continued)'],...
		    'Fontsize',16,'Fontweight','Bold');
	end
end

y      = y - dy;
dy     = dy*1.2;
if iGMsca==2
	text(0,y,sprintf(['Images scaled to an overall grand mean of %g'],GM))
	y = y - dy;
elseif iGMsca==3
	text(0,y,sprintf(...
		['Images scaled to ',lower(sStud),' grand means of %g'],GM))
	y = y - dy;
end
text(0,y,sprintf(...
    'Gray matter threshold is %6.0f%% of the whole brain mean',THRESH*100))
y = y - dy;
tmp = [strvcat(' ',Cnames,Gnames)];
if any(tmp(:,1)=='r'), str=' (except those tagged ''r'')'; else str=''; end
text(0,y,['Covariates',str,' are centered before inclusion in design matrix'])

spm_print


%-Depict and label design matrix, show numbers of parameters.
%=======================================================================
spm_clf(Fgraph); axis off
text(0.30,1.02,'Design Matrix','Fontsize',16,'Fontweight','Bold');

%-Label the effects
%-----------------------------------------------------------------------
hDesMtx = axes('Position',[0.2 0.3 0.6 0.5]);
image((nHCBG + 1)*32);
ylabel('Observations')
xlabel('effects')
hEfLabs = axes('Position',[0.2 0.82 0.6 0.1],'Visible','off');
y     = 0.1;
dx    = 1/size(nHCBG,2);
for i = 1:size(nHCBG,2)
	text((i - 0.5)*dx,y,deblank(HCBGnames(i,:)),...
		'Fontsize',8,'Rotation',90)
end


%-Display parameter summary
%-----------------------------------------------------------------------
hPramAxes = axes('Position',[0.1 0.1 0.8 0.15],'Visible','off');
text(0,1,['Design: ',DesName]);
text(0,.8,['Global normalisation: ',deblank(sGloNorm(iGloNorm,:))]);

text(0,.6,'Parameters:','Fontsize',12,'Fontweight','Bold');
text(0,.4,sprintf(['%d Condition + %d Covariate ',...
	'+ %d Block + %d Confound'],...
	size(H,2),size(C,2),size(B,2),size(G,2)),...
	'Fontsize',10);
text(0,.25,sprintf(['= %d parameters, having %d degrees of freedom, ',...
	'giving %d residual df (%d scans).'],...
	size([H C B G],2),rank([H C B G]),q - rank([H C B G]),q),...
	'Fontsize',10);

spm_print

%-Implement analysis proper
%=======================================================================

%-This is PET data so sigma = 0 (i.e. independent observations) 
% RT is undefined
%-----------------------------------------------------------------------
spm_spm(V,H,C,B,G,CONTRAST,ORIGIN,THRESH*GX,HCBGnames,P,0,[])

%-Clear figure
%-----------------------------------------------------------------------
spm_clf(Finter)
set(Finter,'Name',' ')
spm('Pointer','Arrow')

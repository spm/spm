%****function spm_spm_ui
% Setting up the general linear model
% FORMAT spm_spm_ui
%_______________________________________________________________________
%
% spm_spm_ui configures the design matrix, data specification and
% thresholds that specify the ensuing statistical analysis. These
% arguments are passed to spm_spm that then performs the actual
% analysis. The design matrix defines the experimental design and the
% nature of hypothesis testing to be implemented.  The design matrix
% has one row for each scan and one column for each effect or
% parameter.  These parameters (e.g. mean condition effect, subject or
% block effect, regression slope of rCBF on global CBF etc) are
% estimated in a least squares sense using the general linear model.
% Specific profiles within these parameters are tested using a linear
% compound or CONTRAST with the t statistic.  The resulting map of t
% values constitutes the SPM{t}.  The SPM{t} is then characterized in
% terms of focal or regional differences by assuming that (under the
% null hypothesis) the SPM{t} behaves as a smooth stationary Gaussian
% field.
%
% The effects are designated as (i) of interest or of no interest and (ii)
% levels of a treatment (indicator type variables) or parameters
% (covariates). E.g:
% level of interest        = condition effect
% covariate of interest    = reaction time. symptom severity, perfomance, etc
% level of no interest     = subject effect
% covariate of no interest = global activity, time, age, dose, etc
%
% From the user's perspective it is important to specify the design
% matrix and contrasts correctly.  The design matrix is built when you
% specifiy the number of studies, subjects and conditions.  An experiment
% consists of one or more STUDIES. A STUDY is the repetition of one or
% more CONDITIONS in one or more SUBJECTS.  The number of SUBJECTS and 
% CONDITIONS can very from STUDY to STUDY.  CONDITIONS from one study are
% treated mathematically as different from CONDITIONS of another STUDY (even
% if the same task was used). Each CONDITION can be replicated within a
% subject.
%
% The CONTRAST is simply a list or vector of coefficients that are used
% to test for a pattern of effects.  The number of coefficients (length
% of the CONTRAST) should be the same as the number of effects of interest
% (the number of conditions for each study plus the number of covariates)
% By specifying different contrasts one can effect a wide variety of analyses,
% including subtractive, parametric and factorial tests.
%
% The data should all have the same image and voxel size and these are
% taken from the first image specified.
%
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
%       => embed cell arrays within another cell array
D = struct(...
	'DesName',	'Single-subject: replication of conditions',...
	'n',	[Inf 2 1 1],	'sL',	{{'repl','cond','subj','stud'}},...
	'nCcov',Inf,		'nGcov',Inf,...
	'Hform',		'i2,''+0m'',''cond''',...
	'Bform',		'i3,''-'',''const''',...
	'iGMsca',1,		'GM',[],...
	'b',struct('aTime',1,'aNCCs',0),...
	'iCFI',1,		'iGloNorm',[1 2 3]);
iDD = length(D); %-Above design is the default design definition

%D = [D, struct(...
%	'DesName',	'Single-subject: replicate conditions & covariates',...
%	'n',[Inf Inf 1 1],	'sL',	{{'repl','scan','',''}},...
%	'Hform',	'i2,''+0m'',''scan''',...
%	'iGMsca',1,	'GM',[],...
%	'b',struct('aTime',1,'aCov',1,'aFDCs',0,'aNCCs',0),...
%	'iGloNorm',[1 2 3])];


%-Parameters for all design specs
%-----------------------------------------------------------------------
% Naming contrasts?
% Writing contrasts?
% Selecting scans in time order & specifying conditions as input

%-Option definitions
%-----------------------------------------------------------------------
%-Generic level names
sL = {'sL1','sL2','sL3','sL4'}

%-Covariate by factor interaction options
sCFI = {'no','sL2','sL3','sL4'};

%-Global normalization options
sGloNorm = {	'no global normalisation',...				%-1
		'proportional scaling',...				%-2
		'AnCova {sL2-specific}',...				%-3
		'AnCova {sL3-specific}',...				%-4
		'AnCova {sL4-specific}',...				%-5
		'AnCova'};						%-6

%-Grand mean scaling options
sGMsca = {	'no Grand Mean scaling',...				%-1
		'scaling of sL1. Grand means',...			%-2
		'scaling of sL2. Grand means',...			%-2
		'scaling of sL3. Grand means',...			%-3
		'scaling of sL4. Grand means',...			%-4
		'scaling of overall Grand mean',...			%-5
		'(implicit in PropSca global normalisation)'};		%-6
%-NB: Grand mean scaling by subject is redundent for proportional scaling

%-Adjustment options for AnCova designs (for centering of globals)
%-If Grand mean scaling, then would usually AnCova adjust in a similar
% fashion, i.e. to GM.
sAdjTo = {	'specify...',...					%-1
		'Grand mean (mean of all globals)',...			%-2
		'sL3 grand mean (mean of sL3. globals)',...		%-3
		'sL4 grand mean (mean of sL4. globals)',...		%-4
		'(redundant: not doing AnCova)'};			%-5

%-Variable "decoder"
%-----------------------------------------------------------------------
%-****

%-Get design type
%-----------------------------------------------------------------------
D = D(spm_input('Select design type...',1,'m',{D.DesName}',[],iDD));


%-Get filenames, accrue study, subject, condition & replication indicies
%-----------------------------------------------------------------------
i4 = [];		% Level 4 index (usually study)
i3 = [];		% Level 3 index (usually subject), per l4
i2 = [];		% Level 2 index (usually condition), per l3/l4
i1 = [];		% Level 1 index (usually replication), per l2/l3/l4
P  = {};		% cell array of string filenames

if isinf(D.n(4)), n4 = spm_input(['#',D.sL{4},'.s ?'],'+1');
	else, n4 = D.n(4); end
bL4 = n4>1;

ti2 = '121212...';
GUIpos = spm_input('!NextPos');
for j4  = 1:n4
	spm_input('!SetNextPos',GUIpos);
	sL4P=''; if bL4, sL4P = [D.sL{4},'.',int2str(j4),': ']; end
	if isinf(D.n(3)), n3 = spm_input([sL4P,'#',D.sL{3},'.s ?'],'+1');
		else n3 = D.n(3); end
	bL3 = n3>1;
	
	if D.b.aTime & D.n(2)>1
	    disp('selecting in time order - manually specify conditions')%-**
	    %-NB: This means l2 levels might not be 1:n2
	    GUIpos2 = spm_input('!NextPos');
	    for j3 = 1:n3
		sL3P=[D.sL{3},'.',int2str(j3),': '];
		str = [sL4P,sL3P,'select scans...'];
	 	tP = spm_get(D.n(2)*D.n(1),'.img',{str});
	 	n21 = length(tP);
		str = [sL4P,sL3P,'iCond?'];
		ti2 = spm_input(str,GUIpos2,'c',ti2,n21,D.n(2));
		[tl2,null,j] = unique(ti2);
		tn1 = zeros(size(tl2)); ti1 = zeros(size(ti2));
		for i=1:length(tl2)
			tn1(i)=sum(j==i); ti1(ti2==tl2(i))=1:tn1(i); end
		if isfinite(D.n(1)) & any(tn1~=D.n(1))
			%-#i1 levels mismatches specification in D.n(1)
			error(sprintf('#%s not as pre-specified',D.sL{1})), end
	 	P   = [P;tP];
		i4 = [i4, j4*ones(1,n21)];
		i3 = [i3, j3*ones(1,n21)];
		i2 = [i2, ti2];
		i1 = [i1, ti1];
	    end

	else

	    if isinf(D.n(2)), n2 = spm_input([sL4P,'#',D.sL{2},'.s ?'],'+1');
		else n2 = D.n(2); end
	    bL2 = n2>1;

	    if n2==1 & D.n(1)==1 %-single scan per l3 (subj)
		disp('single scan per l3')%-**
		str = [sL4P,'select scans, ',D.sL{3},' 1-',int2str(n3)];
 		P     = [P;spm_get(n3,'.img',{str})];
		i4 = [i4, j4*ones(1,n3)];
		i3 = [i3, [1:n3]];
		i2 = [i2, ones(1,n3)];
		i1 = [i1, ones(1,n3)];
	    else
		%-multi scan per l3 (subj) case
		disp('multi scan per l3')%-**
		for j3 = 1:n3
			sL3P=''; if bL3, sL3P=[D.sL{3},'.',int2str(j3),': ']; end
			if D.n(1)==1
				%-No l1 (repl) within l2 (cond)
				disp('no l1 within l2')%-**
				str = [sL4P,sL3P,'select scans: ',D.sL{2},...
					 ' 1-',int2str(n2)];
 				P = [P;spm_get(n2,'.img',{str})];
				i4 = [i4, j4*ones(1,n2)];
				i3 = [i3, j3*ones(1,n2)];
				i2 = [i2, 1:n2];
				i1 = [i1, ones(1,n2)];
			else
			    disp('l1 within l2')%-**
			    for j2 = 1:n2
				sL2P='';
				if bL2, sL2P=[D.sL{2},'.',int2str(j2),': ']; end
				str = [sL4P,sL3P,sL2P,' select scans...'];
				tP  = spm_get(D.n(1),'.img',{str});
				n1 = size(tP,1);
				P   = [P;tP];
				i4 = [i4, j4*ones(1,n1)];
				i3 = [i3, j3*ones(1,n1)];
				i2 = [i2, j2*ones(1,n1)];
				i1 = [i1, 1:n1];

			    end 		% (for j2)
			end 			% (if D.n(1)==1)
		end 				% (for j3)
	    end					% (if  n2==1 &...)
	end					% (if D.b.aTime & D.n(2)>1)
end 						% (for j4)

%-Total #observations
%-----------------------------------------------------------------------
nScan = length(P);
clear n1 n2 n3 n4


%-Build Condition (H) and Block (B) partitions
%=======================================================================
eval(['[H,Hnames] = spm_DesMtx(',D.Hform,');'])
if rank(H)==size(H,1), error('Unestimable condition effects'), end
eval(['[B,Bnames] = spm_DesMtx(',D.Bform,');'])
if rank(B)==size(B,1), error('Unestimable block effects'), end

return


%-Covariate partition
%=======================================================================

%-Generate options for factor by covariate interactions
%-----------------------------------------------------------------------
sCFI = sf_estrrep(sCFI,[sL',D.sL']);


% sCovEffInt = strvcat('no','cond','stud','subj');
% iCovEffInt = 1;
% if any(nCOND>1), iCovEffInt = [iCovEffInt,2]; end
% if nStud  > 1, iCovEffInt = [iCovEffInt,3]; end
% if (nSUBJ > 1) & (nSUBJ < q), iCovEffInt = [iCovEffInt,4]; end
% bAskFDCs   = bAskFDCs & (length(iCovEffInt) > 1);


%-Get covariates of interest
%-----------------------------------------------------------------------
C = []; Cnames = ''; Cc = []; Ccnames = '';
if bAskCov
    c = spm_input('# of covariates (of interest)','+1','0|1|2|3|4|5|>',0:6,1);
    if (c == 6), c = spm_input('# of covariates (of interest)','+0'); end
    GUIpos = spm_input('!NextPos');
    while size(Cc,2) < c
        nCcs = size(Cc,2);
        d    = spm_input(sprintf('[%d] - Covariate %d',[q,nCcs+1]),GUIpos);
        if size(d,1) == 1, d = d'; end
        if size(d,1) == q
            %-Save raw covariates for printing later on
            Cc = [Cc,d];
            %-Centre the covariate?
            if bAskNCCs
                tmp=spm_input('Centre covariate(s) ?',GUIpos,'yes|no',[1,0],1);
            else, tmp = 1; end
            if tmp, d  = d - ones(q,1)*mean(d); str=''; else, str='r'; end
            dnames = [str,'CovInt#',int2str(nCcs + 1)];
            if size(d,2) == 1
                %-Single covariate entered - ask about interactions?
                Ccnames = strvcat(Ccnames,dnames);
                if bAskFDCs
                    i = spm_input(['Covariate',int2str(nCcs+1),...
                            ': specific fits'],GUIpos,'b',...
                            sCovEffInt(iCovEffInt,:),iCovEffInt,1);
                    if (i==2) %-Interaction with condition
                        [d,dnames] = spm_DesMtx([iCOND',d],...
                        'FxC',strvcat('Cond',dnames));
                    elseif (i==3) %-Interaction with study
                        [d,dnames] = spm_DesMtx([i4',d],...
                        'FxC',strvcat('Stud',dnames));
                    elseif (i==4) %-Interaction with subject
                        [d,dnames]=spm_DesMtx([iSUBJ',d],...
                        'FxC',strvcat('SUBJ',dnames));
                    end
                end
                C = [C, d];
                Cnames = strvcat(Cnames,dnames);
            else
                %-Block of covariates entered - add to design matrix
                for i = nCcs+1:nCcs+size(d,1)
                    dnames = strvcat(dnames,['rCovInt#',int2str(i)]); end
                Ccnames = strvcat(Ccnames,dnames);
                C = [C, d];
                Cnames = strvcat(Cnames,dnames);
            end % (if)
        end % (if)
    end % (while)
end % (if bAskCov)


%-Get confounding covariates
%-----------------------------------------------------------------------
G = []; Gnames = ''; Gc = []; Gcnames = '';
if bAskCov
    g = spm_input('# of confounding covariates','+1','0|1|2|3|4|5|>',0:6,1);
    if (g == 6), g = spm_input('# of confounding covariates','+0'); end
    GUIpos = spm_input('!NextPos');
    while size(Gc,2) < g
        nGcs = size(Gc,2);
        d = spm_input(sprintf('[%d] - Covariate %d',[q,nGcs + 1]),GUIpos);
        if (size(d,1) == 1), d = d'; end
        if size(d,1) == q
            %-Save raw covariates for printing later on
            Gc = [Gc,d];
            %-Centre the covariate?
            if bAskNCCs
                tmp=spm_input('Centre covariate(s) ?',GUIpos,'yes|no',[1,0],1);
            else, tmp = 1; end
            if tmp, d  = d - ones(q,1)*mean(d); str=''; else, str='r'; end
            dnames = [str,'ConfCov#',int2str(nGcs+1)];
            if size(d,2) == 1
                %-Single covariate entered - ask about interactions
                Gcnames = strvcat(Gcnames,dnames);
                if bAskFDCs
                    iCovEffInt(iCovEffInt==2)=[];
                    i = spm_input(['Confound',int2str(nGcs + 1),...
                            ': specific fits'],GUIpos,'b',...
                            sCovEffInt(iCovEffInt,:),iCovEffInt,1);
                     if (i==3) %-Interaction with study      
                         [d,dnames] = spm_DesMtx([i4',d],...
                         'FxC',strvcat('Stud',dnames));
                     elseif (i==4) %-Interaction with subject
                         [d,dnames]=spm_DesMtx([iSUBJ',d],...
                         'FxC',strvcat('SUBJ',dnames));
                     end
                end
            else
                %-Block of covariates entered - add to design matrix
                for i = nGcs+1:nGcs+size(d,1)
                     dnames = strvcat(dnames,['rConfCov#',int2str(i)]); end
                Gcnames = strvcat(Gcnames,dnames);
            end
            G = [G, d];
            Gnames = strvcat(Gnames,dnames);
        end % (if)
    end % (while)
end % (if bAskCov)


%-Global normalization options
%-----------------------------------------------------------------------
if iGloNorm>9
	%-User has a choice from the options in iGloNorm.
	%-iGloNorm contains an integer, each digit specifies an option
	%---------------------------------------------------------------
	str = int2str(iGloNorm);
	tmp = []; for i = 1:length(str), tmp = [tmp, eval(str(i))]; end
	iGloNorm=spm_input...
	    ('Select global normalisation','+1','m',sGloNorm(tmp,:),tmp);
end

%-Grand mean scaling options
%-----------------------------------------------------------------------
if iGMsca>9
	%-User has a choice from the options in iGMsca.
	%-iGMsca contains an integer, each digit specifies an option
	%---------------------------------------------------------------
	str = int2str(iGMsca);
	tmp = []; for i = 1:length(str), tmp = [tmp, eval(str(i))]; end
	%-Scaling by study redundent if proportional scaling,
	% don't offer study specifics if not bMStud
	if (iGloNorm==2 | ~bMStud) & any(tmp==3), tmp(find(tmp==3))=[]; end
	iGMsca=spm_input...
	    ('Grand mean scaling','+1','m',sGMsca(tmp,:),tmp,max(tmp));
end
if iGMsca>1,	GM = spm_input('Value for grand mean ?','+1','e',50);
		if GM==0, iGMsca=1; end
else, GM=0; end

%-Get threshold defining voxels to analyse
%-----------------------------------------------------------------------
THRESH = spm_input('Gray matter threshold ?','+1','e',0.8);


%-Get contrasts or linear compound for parameters of interest [H C]
%-----------------------------------------------------------------------
a     = size([H C],2);
if a>0, t = spm_input('# of contrasts','+1'); else, t=0; end
CONTRAST = [];
while size(CONTRAST,1) < t
	d = spm_input(sprintf('[%d] - contrast %d',a,size(CONTRAST,1)+1),'+0');
 	if (size(d,2) ~= a), d = d'; end
	if (size(d,2) == a), CONTRAST = [CONTRAST; d]; end
end


%-The interactive parts of spm_spm_ui are now finished
%-----------------------------------------------------------------------
set(Finter,'Name','thankyou')
spm('Pointer','Watch')



%-Computation
%=======================================================================

%-Get file identifiers
%-----------------------------------------------------------------------
V     = zeros(12,q);
for i = 1:q; V(:,i) = spm_map(P(i,:));  end


%-Check for consistency of image size and voxel size
%-----------------------------------------------------------------------
if ~(all(all(~diff(V([1:6],:)'))))
	error('data do not have the same image and voxel size'); end


%-Get ORIGIN
%-----------------------------------------------------------------------
[DIM VOX SCALE TYPE OFFSET ORIGIN] = spm_hread(P(1,:));

%-Compute global values
%-----------------------------------------------------------------------
GX     = zeros(q,1);
for i  = 1:q
	GX(i) = spm_global(V(:,i)); end

%-Grand mean scaling (if required): Scale scaling coefficients so that
% Grand mean (mean of global means) (for each study) is GM
%-----------------------------------------------------------------------
if iGMsca==1
	%-No grand mean scaling
elseif iGMsca==2
	%-Grand mean scaling (overall)
	V(7,:) = V(7,:)*GM/mean(GX);
	GX     = GX*GM/mean(GX);
elseif iGMsca==3
	%-Grand mean scaling by block
	for i = 1:max(i4)
		SStu_d = find(i4==i); SStu_m = mean(GX(SStu_d));
		V(7,SStu_d)  = V(7,SStu_d) *GM/SStu_m;
		GX(SStu_d)   = GX(SStu_d)  *GM/SStu_m;
	end
else, error('Invalid iGMsca option'), end

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

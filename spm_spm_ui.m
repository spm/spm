function spm_spm_ui
% Setting up the general linear model
% FORMAT spm_spm_ui
%____________________________________________________________________________
%
% spm_spm_ui configures the design matrix, data specification and
% thresholds that specify the ensuing statistical analysis. These
% arguments are passed to spm_spm that then performs the actual
% analysis The design matrix defines the experimental design and the
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
% level of interest = condition effect
% covariate of interest = reaction time. symptom severity, perfomance, etc
% level of no interest = subject effect
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
%__________________________________________________________________________
% %W% Andrew Holmes, Karl Friston %E%

%-----------------------------------------------------------------------
clear

%-Get filenames and other user specified parameters
%=======================================================================

%-Get filenames and create design matrix
%-----------------------------------------------------------------------
set(2,'Name','Statistical analysis'); drawnow

%-Design parameters
%-----------------------------------------------------------------------

Designs = str2mat(...
	'SPM94 style questioning',...					%-0
	'Single-subject: replication of conditions',...			%-1a
	'Single-subject: conditions & covariates',...			%-1b
	'Single-Subject: covariates only');				%-1c
Designs = str2mat(Designs,...
	'Multi-subject: different conditions',...			%-2a
	'Multi-subject: with replications',...				%-2b
	'Multi-subject: with replications & covariates',...		%-2c
	'Multi-subject: condition & covariates',...			%-2d
	'Multi-subject: covariates only');				%-2e
Designs = str2mat(Designs,...
	'Multi-study: different conditions',...				%-3a
	'Multi-study: with replications',...				%-3b
	'Multi-study: condition & covariates',...			%-3c
	'Multi-study: covariates only');				%-3d


DesPrams = str2mat(...
	'bMStud','bMSubj','bMCond','bMRepl','iHForm',...
	'bAskCov','bAskFDCs','iGloNorm');

DesDefaults = [ ...
1,	1,	1,	0,	4,	1,	1,	12345;...	%-0
0,	0,	1,	1,	1,	0,	0,	123;...		%-1a
0,	0,	1,	1,	1,	1,	0,	123;...		%-1b
0,	0,	0,	1,	3,	1,	0,	123;...		%-1c
0,	1,	1,	0,	1,	0,	0,	1234;...	%-2a
0,	1,	1,	1,	2,	0,	0,	1234;...	%-2b
0,	1,	1,	1,	2,	1,	1,	1234;...	%-2c
0,	1,	1,	0,	1,	1,	1,	1234;...	%-2d
0,	1,	0,	1,	3,	1,	1,	1234;...	%-2e
1,	1,	1,	0,	4,	0,	0,	12345;...	%-3a
1,	1,	1,	1,	5,	0,	0,	12345;...	%-3b
1,	1,	1,	0,	4,	1,	1,	12345;...	%-3b
1,	1,	0,	1,	3,	1,	1,	12345	];	%-3d

sGloNorm=str2mat(...
	'No Global Normalisation',...
	'Proportional scaling',...
	'AnCova',...
	'AnCova {subject-specific}',...
	'AnCova {study-specific}');

HForms = str2mat(...
	'iCond,''-'',''Scan''',...
	'iCond,''-'',''Cond''',...
	'[]',...
	'[iStud'',iCond''],''-'',str2mat(''Stud'',''Scan'')',...
	'[iStud'',iCond''],''-'',str2mat(''Stud'',''Cond'')');


if ( size(Designs,1) ~= size(DesDefaults,1)) | ...
	(size(DesPrams,1) ~= size(DesDefaults,2))
	fprintf('%cSize mismatch in design parameter specification\n',7)
	return
end

%-Variable "decoder"
%-----------------------------------------------------------------------
% bMStud   - Multi-study ?
% bMSubj   - Multi-Subject ?
% bMCond   - Multi-Condition ? (No for correlation studies)
% bMRepl   - Multiple replications per condition ?
% iHForm   - Index to HForm for this design
% bAskCov  - Ask for covariates and confounds?
% iGloNorm - Global normalisation code, or allowable codes
% bAskFDCs - Ask for [confounding] factor dependent covariate modelling
%            (Seperate slope parameter for each subject/study)
% sStud    - study number string (for i/o)
% sSubj    - subject number string (for i/o)
% nStud    - number of studies
% nSubj    - number of subjects (current study)
% bBalCond - balanced conditions (scans) per subject (current study)
% nCond    - number of conditions (for subjects within current study)
% subj     - current subject index (within current study)
% cond     - current condition (within current subject)

%-Initialise indicies, and get study parameters
%-----------------------------------------------------------------------
iStud   = [];		% study index
iSubj   = [];		% subject (per study) index
iCond   = [];		% condition (or scan) (per subject) index
iRepl   = [];		% replication (per condition) index
P       = [];		% string matrix of filenames
J       = 1;		% Position of user interafce input

DesType = spm_input('Select design type...',J,'m',Designs); J = J + 1;
DesName = deblank(Designs(DesType,:));
for p   = 1:size(DesPrams,1)
    eval([deblank(DesPrams(p,:)),' = DesDefaults(DesType,p);']), end
HForm   = HForms(iHForm,:);


%-Get filenames, accrue study, subject, condition & replication indicies
%-----------------------------------------------------------------------
nStud     = 1;
if bMStud
    nStud = spm_input('Number of studies ?',J); J = J + 1; end
J0        = J;
for stud  = 1:nStud
	J     = J0;
	sStud = []; if bMStud sStud = ['Study ',int2str(stud),': ']; end
	nSubj = 1;
	if bMSubj nSubj = spm_input([sStud,'# of subjects ?'],J); J = J + 1; end

	if bMCond
		nCond = spm_input([sStud,'# of conditions ? '],J);
		J = J + 1;
	else
		nCond = 1;
	end

	if nCond == 1 & ~bMRepl %-Single scan per subject case
		if nSubj<=1
		    fprintf('%cOnly 1 subj. 1cond & no repl?\n',7)
		    return
		end

		t_str = [sStud,'Select scans, subjects 1 - ',int2str(nSubj)];
 		P     = str2mat(P,spm_get(nSubj,'.img',t_str));
		iStud = [iStud, stud*ones(1,nSubj)];
		iSubj = [iSubj, [1:nSubj]];
		iCond = [iCond, ones(1,nSubj)];
		iRepl = [iRepl, ones(1,nSubj)];

	else %-Multi scan per subject case

		for subj      = 1:nSubj
			sSubj = [];
			if bMSubj sSubj = ['Subj.',int2str(subj),': ']; end
			if ~bMRepl %-No replications within conditions
				t_str = [sStud,sSubj,'Select scans ',...
					 'conditions 1 -',int2str(nCond)];
 				P = str2mat(P,spm_get(nCond,'.img',t_str));
				iStud = [iStud, stud*ones(1,nCond)];
				iSubj = [iSubj, subj*ones(1,nCond)];
				iCond = [iCond, 1:nCond];
				iRepl = [iRepl, ones(1,nCond)];
			else

			    for cond  = 1:nCond
				t_str = [];
				if nCond > 1
					t_str = sprintf('Condition %d:',cond);
				end
				t_str = [sStud,sSubj,t_str,' Select scans...'];
				tP    = spm_get(Inf,'.img',t_str);
				nRepl = size(tP,1);
				P     = str2mat(P,tP);
				iStud = [iStud, stud*ones(1,nRepl)];
				iSubj = [iSubj, subj*ones(1,nRepl)];
				iCond = [iCond, cond*ones(1,nRepl)];
				iRepl = [iRepl, 1:nRepl];

			    end 		% (for cond)

			end 			% (if ~bMRepl)
		end 				% (for subject)
	end 					% (if  nCond==...)
end 						% (for study)

P(1,:)  = [];

%-clear working variables
%-----------------------------------------------------------------------
clear sStud sSubj nSubj bBalCond nCond subj cond

%-Total #observations
%-----------------------------------------------------------------------
q       = length(iStud);

%-Create iSUBJ & iCOND indicators, indicating subjects and conditions
% uniquely across studies. Watch out for unbalanced designs!
%-Construct nSUBJ & nCOND as total numbers of subjects and conditions
% across studies and studies & subjects respectively.
%-----------------------------------------------------------------------
nSubj	= iSubj([diff(iStud),1]);   			%-#subject per study
nSUBJ	= sum(nSubj);               			%-#subjects in total
temp	= cumsum([0,nSubj]);
iSUBJ	= iSubj+temp(cumsum([1,diff(iStud)])); 		%-Index to subjects
nCond   = []; for stud=1:nStud, nCond=[nCond,max(iCond(iStud==stud))]; end
nCOND	= sum(nCond);               			%-#conditions in total
temp	= cumsum([0,nCond]);
iCOND   = iCond+temp(cumsum([1,diff(iStud)])); 		%-Index to conditions


%-Build Constant, Condition and Block partitions
%=======================================================================

%-Condition partition
%-----------------------------------------------------------------------
eval(['[H,Hnames] = spm_DesMtx(',HForm,');'])
if (nCOND == 1) | (nCOND == q)			%-Unestimable effects
	H = []; Hnames = ''; end

%-Include a constant term if empty H partition - design subpartition K
%-----------------------------------------------------------------------
if isempty(H)
	K = ones(q,1); Knames = 'Constant';
else
	K = []; Knames = '';
end

%-Always model block (subject) effects if possible
%-----------------------------------------------------------------------
if (nSUBJ == 1) | (nSUBJ == q)		%-No | Unestimable effects
	B = []; Bnames = [];
else
	%-Use implicit SumToZero constraints via relative block effects & pinv.
	%-See spm_DesMtx for more information on this.
	[B,Bnames] = spm_DesMtx(iSUBJ,'+0m','SUBJ');
end % (if)



%-Covariate partition
%=======================================================================

%-Generate options for factor by covariate interactions
%-----------------------------------------------------------------------
sCovEffInt = str2mat('no','study','subject');
iCovEffInt = 1;
if nStud  > 1, iCovEffInt = [iCovEffInt,2]; end
if (nSUBJ > 1) & (nSUBJ < q), iCovEffInt = [iCovEffInt,3]; end
bAskFDCs   = bAskFDCs & (length(iCovEffInt) > 1);


%-Get covariates of interest
%-----------------------------------------------------------------------
C = []; Cnames = ''; Cc = []; Ccnames = '';
if bAskCov
    c = spm_input('# of covariates (of interest)',J,'0|1|2|3|4|5|>',0:6);
    if (c == 6), c = spm_input('# of covariates (of interest)',J); end
    J = J + 1;
    while size(Cc,2) < c
        nCcs = size(Cc,2);
        d = spm_input(sprintf('[%d] - Covariate %d',[q,nCcs+1]),J);
        if size(d,1) == 1, d = d'; end
        if size(d,1) == q
            if spm_input('Centre this covariate ?',J,'yes|no',[1 0]);
                d = d - ones(q,1)*mean(d); end
            dnames = ['CovInt#',int2str(nCcs + 1)];
            if size(d,2) == 1
                %-Single covariate entered - ask about interactions
                Cc = [Cc,d];
                Ccnames = str2mat(Ccnames,dnames);
                if bAskFDCs
                    i = spm_input(['Covariate',int2str(nCcs+1), ': specific fits'],J,'b',sCovEffInt(iCovEffInt,:),iCovEffInt);
                    if (i==2) %-Interaction with study
                        [d,dnames] = spm_DesMtx([iStud',d],...
                        'FxC',str2mat('Stud',dnames));
                    elseif (i==3) %-Interaction with subject
                        [d,dnames]=spm_DesMtx([iSUBJ',d],...
                        'FxC',str2mat('SUBJ',dnames));
                    end % (if)
                end % (if bAskFDCs)
                C = [C, d];
                Cnames = str2mat(Cnames,dnames);
            else
                %-Block of covariates entered - add to design matrix
                for i = nCcs+1:nCcs+size(d,1)
                    dnames = str2mat(dnames,['CovInt#',int2str(i)]); end
                Cc = [Cc,d];
                Ccnames = str2mat(Ccnames,dnames);
                C = [C, d];
                Cnames = str2mat(Cnames,dnames);
            end % (if)
        end % (if)
    end % (while)
end % (if bAskCov)

%-Strip off blank line from str2mat concatenations
%-----------------------------------------------------------------------
if size(Cc,2), Cnames(1,:) = []; Ccnames(1,:) = []; end


%-Get confounding covariates
%-----------------------------------------------------------------------
G = []; Gnames = ''; Gc = []; Gcnames = '';
if bAskCov
    g = spm_input('# of confounding covariates',J,'0|1|2|3|4|5|>',0:6);
    if (c == 6), g = spm_input('# of confounding covariates',J); end
    J = J + 1;
    while size(Gc,2) < g
        nGcs = size(Gc,2);
        d = spm_input(sprintf('[%d] - Covariate %d',[q,nGcs + 1]),J);
        if (size(d,1) == 1), d = d'; end
        if size(d,1) == q
            if spm_input('Centre this covariate ?',J,'yes|no',[1 0]);
                d = d - ones(q,1)*mean(d); end
            dnames = ['ConfCov#',int2str(nGcs+1)];
            if size(d,2) == 1
                %-Single covariate entered - ask about interactions
                Gc = [Gc,d];
                Gcnames = str2mat(Gcnames,dnames);
                if bAskFDCs
                    i = spm_input(['Confound',int2str(nGcs + 1) ': specific fits'],J,'b',sCovEffInt(iCovEffInt,:),iCovEffInt);
                     if (i==2) %-Interaction with study      
                         [d,dnames] = spm_DesMtx([iStud',d],...
                         'FxC',str2mat('Stud',dnames));
                     elseif (i==3) %-Interaction with subject
                         [d,dnames]=spm_DesMtx([iSUBJ',d],...
                         'FxC',str2mat('SUBJ',dnames));
                     end % (if)
                end % (if bAskFDCs)
            else
                %-Block of covariates entered - add to design matrix
                for i = nGcs+1:nGcs+size(d,1)
                     dnames = str2mat(dnames,['ConfCov#',int2str(i)]); end
                Gc = [Gc,d];
                Gcnames = str2mat(Gcnames,dnames);
            end % (if)
            G = [G, d];
            Gnames = str2mat(Gnames,dnames);
        end % (if)
    end % (while)
end % (if bAskCov)

%-Strip off blank line from str2mat concatenations
%-----------------------------------------------------------------------
if size(Gc,2), Gnames(1,:)=[]; Gcnames(1,:)=[]; end


%-Global normalization options
%-----------------------------------------------------------------------
if iGloNorm>9
	%-User has a choice from the options in iGloNorm.
	%-iGloNorm contains an integer, each digit specifies an allowable option
	str = int2str(iGloNorm);
	tmp = []; for i = 1:length(str), tmp = [tmp, eval(str(i))]; end
	iGloNorm=spm_input...
	    ('Select global normalisation',J,'m',sGloNorm(tmp,:),tmp);
	J   = J + 1;
end


%-Get orientation of images
%-----------------------------------------------------------------------
FLIP   = spm_input('Image left = subject''s ',J,'right|left',[1,0]);


%-Get threshold defining voxels to analyse
%-----------------------------------------------------------------------
THRESH = spm_input('Gray matter threshold ?',J,'e',0.8);


%-Get value to be assigned to grand mean:
% 50 is usual for rCBF. 0 for no scaling
%-----------------------------------------------------------------------
GM    = spm_input('Value for grand mean ?',J,'e',50);
J     = J + 1;

%-Get contrasts or linear compound for parameters of interest [H C]
%-----------------------------------------------------------------------
t     = spm_input('# of contrasts',J);
a     = size([H C],2);
J     = J + 1;

CONTRAST = []; % row matrix of contrasts
while size(CONTRAST,1) < t
	d = spm_input(sprintf('[%d] - contrast %d',a,size(CONTRAST,1) + 1),J);
 	if (size(d,2) ~= a), d = d'; end
	if (size(d,2) == a), CONTRAST = [CONTRAST; d]; end
end % (while)



%-The interactive parts of spm_spm_ui are now finished
%-----------------------------------------------------------------------
set(2,'Name','Thank You','Pointer','Watch')



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
if DIM(3) == 1; ORIGIN = [0 0 0]; FLIP = 0; end

%-Compute global values
%-----------------------------------------------------------------------
GX     = zeros(q,1);
for i  = 1:q
	GX(i) = spm_global(V(:,i)); end

%-Scale scaling coefficients so that Grand mean (mean of global means)
% for each study is GM (if GM~=0)
%-----------------------------------------------------------------------
if GM ~= 0
   for i = 1:max(iStud)
	SStu_d       = find(iStud==i);
	SStu_m	     = mean(GX(SStu_d));
	V(7,SStu_d)  = V(7,SStu_d) *GM/SStu_m;
	GX(SStu_d)   = GX(SStu_d)  *GM/SStu_m;
   end
end

%-Construct Global part of covariates of no interest partition.
%-Centre global means if included in AnCova models, by mean correction.
%-----------------------------------------------------------------------
Gc    = [Gc,GX];
if isempty(Gcnames), Gcnames = 'Global';
    else Gcnames = str2mat(Gcnames,'Global'); end

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
    G = [G,(GX - mean(GX))];
    if isempty(Gnames), Gnames = 'Global';
        else Gnames = str2mat(Gnames,'Global'); end

elseif iGloNorm == 4				%-AnCova by subject
%-----------------------------------------------------------------------
    [GL,GLnames] = spm_DesMtx([iSUBJ',GX-mean(GX)],'FxC',['SUBJ  ';'Global']);
    G = [G,GL];
    if isempty(Gnames), Gnames = GLnames;
        else Gnames = str2mat(Gnames,GLnames); end

elseif iGloNorm == 5				%-AnCova by study
%-----------------------------------------------------------------------
    [GL,GLnames] = spm_DesMtx([iStud',GX-mean(GX)],'FxC',['Stud  ';'Global']);
    G = [G,GL];
    if isempty(Gnames), Gnames = GLnames;
        else Gnames = str2mat(Gnames,GLnames); end
else

    fprintf('%cError: invalid iGloNorm option\n',7)

end % (if)


%-Construct full design matrix and name matrices for display
%-----------------------------------------------------------------------
[nKHCBG,KHCBGnames]=...
	spm_DesMtxSca(K,Knames,H,Hnames,C,Cnames,B,Bnames,G,Gnames);


%-Ensure validity of contrast of condition effects, zero pad
%-----------------------------------------------------------------------
a        = size(CONTRAST,1);
d        = 1:size(H,2);
if d
	for i = 1:a
		CONTRAST(i,d) = CONTRAST(i,d) - mean(CONTRAST(i,d));
	end % (for)
end % (if)

CONTRAST = [zeros(a,size(K,2)), CONTRAST, zeros(a,size([B G],2))];


%-Display analysis parameters
%=======================================================================

%-Compute common path components - all paths will begin with '/'
%-----------------------------------------------------------------------
d     = max(find(P(1,1:min(find(~all(P == ones(q,1)*P(1,:))))-1)=='/')) - 1;
CPath = P(1,1:d);
Q     = P(:,d+1:size(P,2));

%-Display
%-----------------------------------------------------------------------
figure(3); spm_clf; axis off
text(0.30,1.02,'Statistical analysis','Fontsize',16,'Fontweight','Bold');
text(-0.10,0.85,'Scan Index','Rotation',90)
if bMStud, text(-0.05,0.85,'Study',      'Rotation',90); end
if bMSubj, text(+0.00,0.85,'Subject',    'Rotation',90); end
if bMCond, text(+0.05,0.85,'Condition',  'Rotation',90); end
if bMRepl, text(+0.10,0.85,'Replication','Rotation',90); end
x0    = 0.15; y0 = 0.83;
dx    = 0.15; dy = 0.02;
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
	if bMStud, text(-0.06,y,sprintf('%2d',iStud(i))); end
	if bMSubj, text(-0.01,y,sprintf('%2d',iSubj(i))); end
	if bMCond, text(+0.04,y,sprintf('%2d',iCond(i))); end
	if bMRepl, text(+0.09,y,sprintf('%2d',iRepl(i))); end
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
		spm_clf; axis off
		y = y0;
		text(0.16,1.02,['Statistical analysis (continued)'],...
		    'Fontsize',16,'Fontweight','Bold');
	end % (if y < 0)
end % (for i = 1:q)

y      = y - dy;
dy     = dy*1.2;
if (GM~=0)
	text(0,y,sprintf(['Images scaled to a grand mean of %g'],GM))
	y = y - dy;
end
if (FLIP)
	text(0,y,'The data will be flipped in accord with neurological')
	y = y - dy;
	text(0,y,['convention - left corresponds to the subjects left'])
	y = y - dy;
end
text(0,y,sprintf(...
    'Gray matter threshold is %6.0f%% of the whole brain mean',THRESH*100))

spm_print


%-Depict and label design matrix, show numbers of parameters.
%=======================================================================
figure(3); spm_clf; axis off
text(0.30,1.02,'Design Matrix','Fontsize',16,'Fontweight','Bold');

%-Label the effects
%-----------------------------------------------------------------------
hDesMtx = axes('Position',[0.2 0.3 0.6 0.5]);
image((nKHCBG + 1)*32);
ylabel('Observations')
xlabel('effects')
hEfLabs = axes('Position',[0.2 0.82 0.6 0.1],'Visible','off');
y     = 0.1;
dx    = 1/size(nKHCBG,2);
for i = 1:size(nKHCBG,2)
	text((i - 0.5)*dx,y,deblank(KHCBGnames(i,:)),...
		'Fontsize',8,'Rotation',90)
end % (for)


%-Display parameter summary
%-----------------------------------------------------------------------
hPramAxes=axes('Position',[0.1 0.1 0.8 0.15],'Visible','off');
text(0,1,['Design: ',DesName]);
text(0,.8,['Global normalisation: ',deblank(sGloNorm(iGloNorm,:))]);

text(0,.6,'Parameters:','Fontsize',12,'Fontweight','Bold');
text(0,.4,sprintf(['%d (Constant) + %d Condition + %d Covariate ',...
	'+ %d Block + %d Confound'],...
	size(K,2),size(H,2),size(C,2),size(B,2),size(G,2)),...
	'Fontsize',10);
text(0,.25,sprintf(['= %d parameters, having %d degrees of freedom, ',...
	'giving %d residual df (%d scans).'],...
	size([K H C B G],2),rank([K H C B G]),q-rank([K H C B G]),q),...
	'Fontsize',10);

spm_print

%-Implement analysis proper
%=======================================================================
Ut = spm_invNcdf(1 - 0.01);
spm_spm(V,K,H,C,B,G,CONTRAST,Ut,ORIGIN,THRESH*GX,FLIP,KHCBGnames);


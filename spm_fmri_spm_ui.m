function spm_fmri_spm_ui
% Setting up the general linear model for fMRI time-series
% FORMAT spm_fmri_spm_ui
%____________________________________________________________________________
%
% spm_fmri_spm_ui configures the design matrix, data specification and
% thresholds that specify the ensuing statistical analysis. These
% arguments are passed to spm_spm that then performs the actual analysis.
%
% This routine has two uses:
%
% i)  To implement a unvariate analysis testing for effects of interest
%     at each and every voxel.  Effects of interest are specified as
%     covariates of interest and confounds are specified as covariates
%     of no interest.
%
% ii) To remove any confounds as a preprocessing step for further analysis
%      e.g. SVD.  In this instance only confounds are specified and the
%      adjusted data will be found in XA.mat
%
% The design matrix defines the experimental design and the nature of
% hypothesis testing to be implemented.  The design matrix has one row
% for each scan and one column for each effect or parameter.  These
% parameters (e.g. reference waveform, subject or block effect, 
% regression slope of voxel on a confound etc) are estimated in a
% least squares sense using the general linear model.  Specific profiles
% within these parameters are tested using a linear compound or CONTRAST
% with the t statistic.  The resulting map of t values constitutes the
% SPM{t}.  The SPM{t} is then characterized in terms of focal or regional
% differences by assuming that (under the null hypothesis) the SPM{t}
% behaves as a smooth stationary Gaussian field.
%
%     From the user's perspective it is important to specify the design
% matrix and contrasts correctly.  The design matrix is built when you
% specify the number of sessions/subjects and conditions.  The covariates 
% (that constitute the columns of the design matrix) can be thought of as
% reference vectors and can be specified as such.  Alternatively one
% can specify reference vectors in terms of response functions or
% waveforms: Waveforms are specified for each EPOCH of scans that
% constitute a particular condition.  Note that if there are two
% conditions, then two covariates are specified each expressing the
% same waveform[s] every time the condition occurs.  A waveform is
% simply a transient response to the onset of an epoch, that lasts for
% the duration of that epoch.  The form of this response waveform can
% be fixed (e.g. a box-car or half sine wave) or allowed to vary
% between conditions by using two (exponentially modulated sine)
% waveforms.  In the latter case there are two covariates for each
% condition (and the CONTRAST must be specified with this in mind).
% These two response functions correspond to the early and late
% components of a modeled response and differential adaptation of the
% hemodynamics during condition-specific epochs can be tested with the
% appropriate CONTRAST (see the final reference below) Covariates of no
% interest (called confounds) can also be specfied.  You will be
% prompted for some specific confounds such as low frequency artifacts
% and whole brain activity.
%
% Epochs can vary in length (and order) within and between subjects or runs.
% If multiple subjects or sessions are specified, then subject or run-specific
% waveforms are used.  This means that main effects of conditions and
% interactions between conditions and subjects (or runs) can be evaluated
% with the appropriate contrast.  If you want to treat all your sessions (or
% subjects) as one then specify just one session/subject.
%
% The CONTRAST is simply a list or vector of coefficients that are used
% to test for a pattern of effects.  The number of coefficients (length
% of the CONTRAST) should be the same as the number of covariates of interest
% By specifying different contrasts one can effect a wide variety of analyses.
%
% Refs:
%
% Friston KJ, Holmes A, Poline J-B, Grasby PJ, Williams SCR, Frackowiak
% RSJ & Turner R (1995) Analysis of fMRI time-series revisited. NeuroImage
% 2:45-53
%
% Worsley KJ and Friston KJ (1995) Analysis of fMRI time-series revisited -
% again. NeuroImage 2:178-181
%
% Friston KJ, Frith CD, Frackowiak RSJ, & Turner R (1995) Characterising
% dynamic brain responses with fMRI: A multivariate approach NeuroImage -
% 2:166-172
%
% Frith CD, Turner R & Frackowiak RSJ (1995) Characterising evoked 
% hemodynamics with fMRI Friston KJ, NeuroImage 2:157-165
%
%__________________________________________________________________________
% %W% Karl Friston, Jean-Baptiste Poline %E%


%----------------------------------------------------------------------------


% get filenames and other user specified parameters
%============================================================================

% get filenames and create design matrix
%----------------------------------------------------------------------------
Finter = spm_figure('FindWin','Interactive');
Fgraph = spm_figure('FindWin','Graphics');
set(Finter,'Name','fMRI analysis');

Q     = [];				% matrix of filename strings
X     = [];				% data matrix {q x pixels}
H     = [];				% Factors of interest
C     = [];				% covariates of interest
B     = [];				% Factors of no interest
G     = [];				% covariates of no interest
F     = [];				% Low frequency confounds
CONTRAST = [];				% row matrix of contrasts


% get Repeat time
%----------------------------------------------------------------------------
RT     = spm_input('Interscan interval {secs}',1);


% set value to be assigned to the rand mean; 100 is usual for fMRI %----------------------------------------------------------------------------
GM     = 100;


% get filenames
%----------------------------------------------------------------------------
nsubj  = spm_input(['# of subjects or sessions'],'!+1','e',1);
nscan  = zeros(1,nsubj);
for i  = 1:nsubj
	str      = sprintf('select scans for session %0.0f',i);
	P        = spm_get(Inf,'.img',str);
 	Q        = str2mat(Q,P);
	nscan(i) = size(P,1);
end
Q(1,:) = []; P = Q;
q      = size(P,1);


% design matrix effect names
%---------------------------------------------------------------------------
Hnames = ' ';
Bnames = ' ';
Cnames = ' ';
Gnames = ' ';


% get ORIGIN
%----------------------------------------------------------------------------
[DIM VOX SCALE TYPE OFFSET ORIGIN DESCRIP] = spm_hread(P(1,:));


% covariates of interest - Type
%----------------------------------------------------------------------------
CovF  = str2mat(...
	'User specified',...
	'2 temporal basis functions',...
	'half sine-wave',...
	'delayed box-car');
Cov   = spm_input('Select type of response','!+1','m',CovF,[1:size(CovF,1)]);

% if box-car ask for temporal differences
%----------------------------------------------------------------------------
TD    = 0;
if Cov == 4
	TD = spm_input('add temporal derivative','!+1','b','no|yes',[0 1]);
end


% for each subject
%---------------------------------------------------------------------------
for v = 1:nsubj

    % reset name
    %-----------------------------------------------------------------------
    set(Finter,'Name',sprintf('Session or subject %0.0f',v));

    % conditions for this session
    %-----------------------------------------------------------------------
    k = nscan(v);

    % user specified
    %-----------------------------------------------------------------------
    if Cov == 1

    	% number of covariates of interest
	%-------------------------------------------------------------------
	c     = spm_input('# of waveforms or conditions',4);
	D     = [];
	while size(D,2) < c
		u   = size(D,2) + 1;
		str = sprintf('[%d]-covariate %d, session %d',k,u,v);
		d   = spm_input(str,'!+1');
		if size(d,2) == k
			d    = d';    end
		if size(d,1) == k
			D    = [D d]; end	
	end
	for i = 1:size(D,2)
		Cnames = str2mat(Cnames,sprintf('Sess %0.0f Cond %0.0f',v,i));
	end

    else

	% vector of conditions
	%-------------------------------------------------------------------
	a     = spm_input('order of epochs eg 1 2 1 2....',4);
	a     = a(:); a = a';
	
	% vector of epoch lengths
	%-------------------------------------------------------------------
	e     = [];
	while length(e) ~= length(a) & all(e > 0) & (sum(e) ~= k)
		e = spm_input('scans/epoch eg 10 or 8 6 8.... ','!+1');
		if length(e) == 1; e = e*ones(1,length(a)); end
	end
	e     = e(:); e = e';

	% onsets for all conditions  
	%-------------------------------------------------------------------
	ons   = [1 (cumsum(e) + 1)]; ons = ons(1:(length(ons) - 1));
	

	% 2 modulated sine waves
	%-------------------------------------------------------------------
	if Cov == 2			
		D = zeros(k,max(a)*2);
		for i = 1:length(e)
		  W = sin(pi*[0:(e(i) + 1)]'/(e(i) + 1)).*exp((-[0:(e(i) + 1)]')/e(i)*4);
		  W = [W sin(pi*[0:(e(i) + 1)]'/(e(i) + 1)).*exp(([0:(e(i) + 1)]')/e(i))];
		  W = W./(ones(size(W,1),1)*max(W));
		  D(([1:size(W,1)] + (ons(i) - 1)),((a(i) - 1)*2 + [1 2])) = W;
		end
		D     = D(1:k,:);
		for i = 1:size(D,2)/2
			str    = sprintf('Sess %0.0f Cond %0.0f e',v,i);
			Cnames = str2mat(Cnames,str);
			str    = sprintf('Sess %0.0f Cond %0.0f l',v,i);
			Cnames = str2mat(Cnames,str);
		end

	
	% single sine wave
	%-------------------------------------------------------------------
	elseif Cov == 3
		D = zeros(k,max(a));
		for i = 1:length(e)
		  W = sin(pi*[0:(e(i) + 1)]'/(e(i) + 1));
		  W = W./ ( ones(size(W,1),1)*max(W) );
		  D(([1:size(W,1)] + (ons(i) - 1)),a(i)) = W;
		end
		D = D(1:k,:);
		for i = 1:size(D,2)
			str    = sprintf('Sess %0.0f Cond %0.0f',v,i);
			Cnames = str2mat(Cnames,str);
		end


	% 6 second delayed box-car
	%-------------------------------------------------------------------
	elseif Cov == 4
		D     = zeros(k,max(a)); 
		for i = 1:length(e)
			D([ons(i):(ons(i) + e(i) - 1)],a(i)) = ones(e(i),1);
		end
		delay = round(6/RT);
		D((delay + 1):k,:) = D(1:(k - delay),:); 
		D(1:delay,:) = zeros(delay,max(a));
		for i = 1:size(D,2)
			str    = sprintf('Sess %0.0f Cond %0.0f',v,i);
			Cnames = str2mat(Cnames,str);
		end
		if TD
			for i = 1:size(D,2)
				str    = sprintf('Cond %0.0f - timing',i);
				Cnames = str2mat(Cnames,str);
			end
			D      = [D [diff(D); zeros(1,size(D,2))]];

		end
	end

	

    end % (else)

    % append to C
    %-----------------------------------------------------------------------
    [x y]  = size(C);
    [i j]  = size(D);
    C(([1:k] + sum(nscan(1:(v - 1)))),([1:j] + y)) = D;

    % append to B
    %-----------------------------------------------------------------------
    [x y]  = size(B);
    B(([1:k] + x),(1 + y))   = ones(k,1);
    Bnames = str2mat(Bnames,sprintf('Session %0.0f',v));

    % covariates not of interest
    %---------------------------------------------------------------------------
    g      = spm_input('# of covariates of no interest','!+1','e',0);
    D      = [];
    while size(D,2) < g
	u   = size(D,2) + 1;
	str = sprintf('[%d]-confound %d, session %d',k,u,v);
	d   = spm_input(str,'!+1','e',0);
	if size(d,2) == k
		d = d'; 	end
	if size(d,1) == k
		D = [D d];	end

    end

    % append to G
    %-----------------------------------------------------------------------
    [x y]  = size(G);
    [i j]  = size(D);
    G(([1:k] + sum(nscan(1:(v - 1)))),([1:j] + y)) = D;
    for i = 1:size(D,2)
    	Gnames = str2mat(Gnames,sprintf('Confound %0.0f',i));
    end

end % loop over sessions


% high pass filter using discrete cosine set
%---------------------------------------------------------------------------
if spm_input('high pass filter','!+1','yes|no',[1 0])
	CUT   = spm_input('cut-off period {secs}','!+0','e',120);
	for v = 1:nsubj
		D     = [];
		k     = nscan(v);
		u     = [1:k/2];
		u     = find(u <= 2*(k*RT)/CUT);
		for i = 1:length(u)
			d = cos(pi*[0:(k - 1)]*u(i)/(k - 1));
			D = [D d(:)];
		end

    		% append to G
   		%-----------------------------------------------------------
    		[x y] = size(F);
    		[i j] = size(D);
   		F(([1:i] + x),([1:j] + y))   = D;

	end
end


%-combine confounds and append confound effect names
%---------------------------------------------------------------------------
G      = spm_en([G F]);
for i  = 1:size(F,2);
	Gnames = str2mat(Gnames,sprintf('Low Hz %0.0f',i));
end

% global normalization with ANCOVA as the default
%----------------------------------------------------------------------------
GLOBAL = 'None';
str    = 'global normalization';
if spm_input(str,'!+1','yes|no',[1 0])
	GLOBAL = spm_input(str,'!+0','Scaling|AnCova'); end


% get contrasts or linear compound for parameters of interest [H C]
%---------------------------------------------------------------------------
t         = 0;
if size([H C],2)
	t = spm_input('# of contrasts','!+1'); end

while size(CONTRAST,1) ~= t
	d   = [];
        str = sprintf('[%0.0f] contrast %0.0f',size(C,2),size(CONTRAST,1) + 1);
	while size(d,2) ~= size([H C],2)
		d = spm_input(str,'!+0');
	end
     	CONTRAST = [CONTRAST; d]; end
end

% temporal smoothing
%---------------------------------------------------------------------------
if spm_input('Temporal smoothing','!+1','yes|no',[1 0])
	SIGMA  = sqrt(8)/RT;
else
	SIGMA  = 0;
end



% the interactive parts of spm_spm_ui are now finished
%---------------------------------------------------------------------------
set(Finter,'Name','thankyou','Pointer','Watch')



% get file identifiers and global values
%===========================================================================
V     = zeros(12,q);
for i = 1:q; V(:,i) = spm_map(P(i,:));  end


% check for consistency of image size and voxel size
%---------------------------------------------------------------------------
if ~(all(all(~diff(V([1:3],:)'))))
	error('images are not compatible'); end


% adjust scaling coefficients so that Grand mean - <GM> = GM
%---------------------------------------------------------------------------
GX     = zeros(q,1);
for i  = 1:q
	GX(i) = spm_global(V(:,i)); end


V(7,:) = V(7,:)*GM/mean(GX);
GX     = GX*GM/mean(GX);

if strcmp(GLOBAL,'AnCova'); G = [G GX]; Gnames = str2mat(Gnames,'Global'); end
if strcmp(GLOBAL,'Scaling'); V(7,:) = GM*V(7,:)./GX'; GX = GM*GX./GX;	   end



% mean correct the contrast for condition effects and zero pad
%---------------------------------------------------------------------------
for i = 1:size(CONTRAST,1)
	if length(find(CONTRAST(i,:))) > 1
		CONTRAST(i,:) = CONTRAST(i,:) - mean(CONTRAST(i,:));
	end
end
[i j]    = size(CONTRAST);
CONTRAST = [CONTRAST zeros(i,(size([H C B G],2) - j))];

%-Construct full design matrix and name matrices for display
%---------------------------------------------------------------------------
Hnames(1,:) = [];
Cnames(1,:) = [];
Bnames(1,:) = [];
Gnames(1,:) = [];


%-centre confounds (Not block, which can embody the constant term)
%---------------------------------------------------------------------------
G     = spm_detrend(G,0);

[nHCBG,HCBGnames] = spm_DesMtxSca(H,Hnames,C,Cnames,B,Bnames,G,Gnames);

% display analysis parameters
%===========================================================================
figure(Fgraph); spm_clf; axis off
axes('Position',[0.2 0.5 0.6 0.4])
imagesc(spm_DesMtxSca(nHCBG)')
set(gca,'FontSize',8)
set(gca,'YTick',[1:size(nHCBG)],'YTickLabels',HCBGnames)
xlabel('scan')
title(['Design Matrix'],'Fontsize',16,'Fontweight','Bold')

axes('Visible','off')
text(0.00,0.40,...
	['AnCova for ' spm_str_manip(pwd,'a24')],'Fontsize',16);
text(0.00,0.32,'Session, scan and Filename');
y     = 0.28;
dy    = 0.03;
if size(B,2)
	for i = 1:size(B,2)
		u = min(find(B(:,i) > 0));
		v = max(find(B(:,i) > 0));
		d = sprintf('%d  %d to %d ',i,u,v);
		text(0,y,d,'FontSize',10)
		text(0.2,y, [P(u,:) '...'],'FontSize',8)
		y = y - dy;
	end
else
	d = sprintf('%d  %d to %d ',1,1,size(P,1));
	text(0,y,d,'FontSize',10)
	text(0.2,y, [P(1,:) '...'],'FontSize',10)
	y = y - dy;
end
y = y - dy;
if SIGMA
	text(0,y,['Temporal Smoothing {2.8sec Gaussian Kernel}']); y = y - dy;
end
if size(F,2)
	str = sprintf('Hi-Pass cutoff %0.1f cycles/min',60/CUT);
	text(0,y,str); y = y - dy;
end
text(0,y,sprintf('Grand Mean = %0.2f a.u.',GM)); y = y - dy;
text(0,y,sprintf('Interscan Interval = %0.2f secs',RT)); y = y - dy;
text(0,y,['Response Form: ' CovF(Cov,:)  ]); y = y - dy;
text(0,y,['Global Normalization: ' GLOBAL]); y = y - dy;


spm_print


% implement analysis proper
%--------------------------------------------------------------------------
spm_spm(V,H,C,B,G,CONTRAST,ORIGIN,GX,HCBGnames,P,SIGMA,RT);

%-Clear figure
%-----------------------------------------------------------------------
spm_clf(Finter); set(Finter,'Name',' ','Pointer','Arrow');

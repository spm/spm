function spm_fmri_spm_ui
% Setting up the general linear model for fMRI time-series
% FORMAT spm_fmri_spm_ui
%____________________________________________________________________________
%
% spm_fmri_spm_ui configures the design matrix, data specification and
% thresholds that specify the ensuing statistical analysis. These
% arguments are passed to spm_fmri_spm that then performs the actual analysis.
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
% regression slope of voxel on global activity etc) are estimated in a
% least squares sense using the general linear model.  Specific profiles
% within these parameters are tested using a linear compound or CONTRAST
% with the t statistic.  The resulting map of t values constitutes the
% SPM{t}.  The SPM{t} is then characterized in terms of focal or regional
% differences by assuming that (under the null hypothesis) the SPM{t}
% behaves as a smooth stationary Gaussian field.
%
%     From the user's perspective it is important to specify the design
% matrix and contrasts correctly.  The design matrix is built when you
% specifiy the number of subjects and conditions.  The covariates (that
% constitute the columns of the design matrix) can be thought of as
% reference vectors and can be specified as such.  Alternatively one
% can specify reference vectors in terms of response functions or
% waveforms: Waveforms are specified for each epoch of scans that
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
% The CONTRAST is simply a list or vector of coefficients that are used
% to test for a pattern of effects.  The number of coefficients (length
% of the CONTRAST) should be the same as the number of covariates of interest
% By specifying different contrasts one can effect a wide variety of analyses.
%
% Refs:
%
% Friston KJ, Holmes A, Poline J-B, Grasby PJ, Williams SCR, Frackowiak
% RSJ & Turner R (1995) Analysis of fMRI time-series revisited. NeuroImage
% - in press
%
% Worsley KJ and Friston KJ (1995) Analysis of fMRI time-series revisited -
% again. NeuroImage submitted
%
% Friston KJ, Frith CD, Frackowiak RSJ, & Turner R (1995) Characterising
% dynamic brain responses with fMRI: A multivariate approach NeuroImage -
% in press
%
% Frith CD, Turner R & Frackowiak RSJ (1995) Characterising evoked 
% hemodynamics with fMRI Friston KJ, NeuroImage - in press
%
%__________________________________________________________________________
% %W% Karl Friston, Jean-Baptiste Poline %E%


%----------------------------------------------------------------------------


% get filenames and other user specified parameters
%============================================================================

% get filenames and create design matrix
%----------------------------------------------------------------------------
set(2,'Name','Statistical analysis'); drawnow

Q     = [];				% matrix of filename strings
X     = [];				% data matrix {q x pixels}
H     = [];				% Factors of interest
C     = [];				% covariates of interest
B     = [];				% Factors of no interest
G     = [];				% covariates of no interest
T     = [];				% row matrix of contrasts
K     = [];				% string matrix for graphical output

CONTRAST = [];

% temporal smoothing
%----------------------------------------------------------------------------
SIGMA  = 1;

% get Repeat time
%----------------------------------------------------------------------------
RT     = spm_input('Interscan interval {secs}',1);


% set value to be assigned to the rand mean; 100 is usual for fMRI %----------------------------------------------------------------------------
GM     = 100;


% cycle over studies, subjects and conditions
%----------------------------------------------------------------------------
n      = spm_input(['# of subjects'],1);
k      = spm_input(['# of scans/subject'],2);
q      = n*k;
for j  = 1:n
	d = ['select scans 1 - ' num2str(k) ' {subject ' num2str(j) '}'];
 	Q = str2mat(Q,spm_get(k,'.img',d));
end

Q(1,:) = []; P = Q;

% check that the data have been normalized and get ORIGIN
%----------------------------------------------------------------------------
[DIM VOX SCALE TYPE OFFSET ORIGIN DESCRIP] = spm_hread(P(1,:));
if DIM(3) == 1; ORIGIN = [0 0 0]; end

% design matrix subpartitions - Block or subject
%---------------------------------------------------------------------------
B     = []; for i = 1:k; B = [B; eye(n)]; end
B     = reshape(B,n,q)';
B     = (B' - ones(n,1)*mean(B'))';

% remove any redundant Block partition
%---------------------------------------------------------------------------
if (size(B,2) == q) | (size(B,2) == 1); B = []; end

% constant term
%---------------------------------------------------------------------------
K     = ones(q,1);


% covariates of interest - Type
%----------------------------------------------------------------------------
CovF  = str2mat(...
	'User specified',...
	'2 temporal basis functions',...
	'half sine-wave',...
	'delayed box-car');

Cov   = spm_input('Select type of response',3,'m',CovF,[1:size(CovF,1)]);


% for each subject
%---------------------------------------------------------------------------
for v = 1:n
	
	keyboard;

    % user specified
    %-----------------------------------------------------------------------
    if Cov == 1

    	% number of covariates of interest
	%-------------------------------------------------------------------
	c     = spm_input('# of covariates or conditions',4);
	D     = [];
	while size(D,2) < c
		u   = size(D,2) + 1;
		str = sprintf('[%d]-covariate %d, subject %d',k,u,v);
		d   = spm_input(str,5);
		if size(d,2) == k
			d = d'; 	end
		if size(d,1) == k
			D   = [D d(:)];	end		
	end

    else

	% vector of conditions
	%-------------------------------------------------------------------
	a = spm_input('conditions {1 to n} eg 1 2 1 2....',5);
	a = a(:); a = a';
	
	% vector of epoch lengths
	%-------------------------------------------------------------------
	e	= [];
	while length(e) ~= length(a) & all(e>0)
		e = spm_input(sprintf('Subj. %d: scans/epoch eg 10 ',v),5);
		if length(e) == 1; e = e*ones(1,length(a)); end
	end

	% onsets for all conditions  
	%-------------------------------------------------------------------
	ons  	= [1 cumsum(e)+1]; ons = ons(1:(length(ons)-1));
	

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
		D = D(1:k,:);
	
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

	% 6 second delayed box-car
	%---------------------------------------------------------------
	elseif Cov == 4
		D = zeros(k,max(a)); 
		for i = 1:length(e)
			D(ons(i):ons(i)+e(i)-1,a(i)) = ones(e(i),1);
		end
		delay = round(6/RT);
		D(delay+1:k,:) = D(1:k-delay,:); 
		D(1:delay,:) = zeros(delay,max(a));
	end

	

    end %-- else --%
    C = [C; D];
end


% covariates not of interest
%----------------------------------------------------------------------------
g     = spm_input('# of covariates of no interest',5);
while size(G,2) < g
	d = spm_input(sprintf('[%0.0f]- covariate %0.0f',q,size(G,2) + 1),6);
	if size(d,2) == q
		d = d'; 	end
	if size(d,1) == q
    		G = [G d];	end
end


% high pass filter
%----------------------------------------------------------------------------
F     = [];
if spm_input('high pass filter',6,'yes|no',[1 0])
      d     = spm_input('cut-off period {secs}',6);
      u     = [1:k/2];
      u     = find(u < 2*(k*RT)/d);
      for v = 1:length(u)
              d = sin(pi*[1:k]*u(v)/k);
              F = [F d(:)];
              d = cos(pi*[1:k]*u(v)/k);
              F = [F d(:)];
      end
end



%-replicate this partition for each subject and combine confounds
% Specified confounds, low frequency components and constant
%---------------------------------------------------------------------------
d      = F;   F = [];
for i  = 1:n; F = [F; d]; end
G      = [G F K];


% global normalization with ANCOVA as the default
%----------------------------------------------------------------------------
GLOBAL = 'None';
str    = 'global normalization';
if spm_input(str,7,'yes|no',[1 0])
	GLOBAL = spm_input(str,7,'AnCova|Scaling'); end


% get contrasts or linear compound for parameters of interest [H C]
%----------------------------------------------------------------------------
t       = 0;
if size([H C],2)
	t = spm_input('# of contrasts',8); end

while size(CONTRAST,1) ~= t
        d = sprintf('[%0.0f] contrast %0.0f',size(C,2),size(CONTRAST,1) + 1);
	d = spm_input(d,9);
	d = d(1:size([H C],2));
	d = d(:)';
     	CONTRAST = [CONTRAST; d]; end
end



% the interactive parts of spm_spm_ui are now finished
%----------------------------------------------------------------------------
set(2,'Name','thankyou','Pointer','Watch')



% get file identifiers and global values
%============================================================================
V     = zeros(12,q);
for i = 1:q; V(:,i) = spm_map(P(i,:));  end


% check for consistency of image size and voxel size
%----------------------------------------------------------------------------
if ~(all(all(~diff(V([1:3],:)'))))
	error('images are not compatible'); end


% adjust scaling coefficients so that Grand mean - <GM> = GM
%----------------------------------------------------------------------------
GX     = zeros(q,1);
for i  = 1:q
	GX(i) = spm_global(V(:,i)); end


V(7,:) = V(7,:)*GM/mean(GX);
GX     = GX*GM/mean(GX);

if strcmp(GLOBAL,'AnCova');  G      = [G GX];        			end
if strcmp(GLOBAL,'Scaling'); V(7,:) = GM*V(7,:)./GX'; GX = GM*GX./GX;	end



% mean correct the contrast for condition effects and zero pad
%---------------------------------------------------------------------------
for i = 1:size(CONTRAST,1)
	if length(find(CONTRAST(i,:))) > 1
		CONTRAST(i,:) = CONTRAST(i,:) - mean(CONTRAST(i,:));
	end
end
d     = size(CONTRAST,1);
CONTRAST = [CONTRAST zeros(d,size([B G],2))];


% display analysis parameters
%============================================================================
global CWD
figure(3); spm_clf; axis off
axes('Position',[0.3 0.6 0.4 0.3])
imagesc(spm_DesMtxSca([H C B G]))
xlabel('parameters or effects')
ylabel('scan')
title(['Design Matrix'],'Fontsize',16,'Fontweight','Bold')

axes('Visible','off')
text(0.20,0.36,'Filename');
text(0.00,0.50,['AnCova for ' CWD],'Fontsize',16);
text(0.00,0.36,'Subject',   'Rotation',90);
text(0.06,0.36,'Scan',      'Rotation',90);
y     = 0.32;
dy    = 0.03;
if size(B,2)
	for i = 1:size(B,2)
		u = min(find(B(:,i) > 0));
		v = max(find(B(:,i) > 0));
		d = sprintf('%d  %d to %d ',i,u,v);
		text(0,y,d,'FontSize',10)
		text(0.2,y, [P(u,:) '...'],'FontSize',10)
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
	text(0,y,['High Pass Filtering {0.5 cycles per minute}']); y = y - dy;
end
text(0,y,sprintf('Grand Mean = %0.2f a.u.',GM)); y = y - dy;
text(0,y,sprintf('Interscan Interval = %0.2f secs',RT)); y = y - dy;
text(0,y,['Response Form: ' CovF(Cov,:)  ]); y = y - dy;
text(0,y,['Global Normalization: ' GLOBAL]); y = y - dy;




spm_print


% implement analysis proper
%---------------------------------------------------------------------------
spm_fmri_spm(V,H,C,B,G,CONTRAST,ORIGIN,GX,RT,SIGMA);


function PPI = spm_peb_ppi(SPM)
% computes either physiophysiologic or psychophysiologic interactions
% FORMAT PPI = spm_peb_ppi(SPM);
%
% SPM    - Structure containing generic details about the analysis
%
% PPI.ppi    = (PSY*xn  or xn1*xn2) convolved with the HRF
% PPI.P      = PSY convolved with HRF. Use as covariate of no interest
% PPI.Y      = Original first eigenvariate bold signal. 
% PPI.name   = Name of PPI
% PPI.xY     = Original VOI information
% PPI.xn     = Deconvolved neural signal(s)
% PPI.U.u    = Psychological variable or input function (PPIs only)
% PPI.U.w    = Contrast weights for psychological variable (PPIs only)
% PPI.U.name = Names of psychological conditions (PPIs only)
%---------------------------------------------------------------------
%
% SETTING UP A PPI THAT ACCOUNTS FOR THE HRF
% ==================================================================
% PPI's were initially conceived as a means of identifying regions whose
% reponses can be explained in terms of an interaction between activity in
% a specified source (the physiological factor) and some experimental
% effect (the psychological factor). However, a problem in setting up PPI's
% is that in order to derive a proper estimate of the interaction between
% a psychological variable (P) and measured hemodynamic signal (x), one cannot
% simply convolve the psychological variable with the hrf (HRF) and multiply
% by the signal. Thus:
% 
%                  conv(P,HRF).* x ~= conv((P.*xn),HRF)
%
% P   = psychological variable
% HRF = hemodynamic response function
% xn  = underlying neural signal which in fMRI is convolved with the hrf to
%       give the signal one measures -- x.
% x   = measured fmri signal
%
% It is actually the right hand side of the equation one wants.
% Thus one has to work backwards, in a sense, and deconvolve the hrf
% from x to get xn. This can them be multiplied by P and the resulting
% vector (or matrix) reconvolved with the hrf.
%
% This algorithm uses a least squares strategy to solve for xn.
%
% The source's hemodynamics are x = HRF*xn;
%
% Using the constraint that xn should be continous and smooth 
% we can expand x in terms of some temporal basis set (xb)
%
%      xn  = xb*B
%       B  = parameter estimate
%
% The estimator of x is then
%
%       x  = HRF(k,:)*xn
%       x  = HRF(k,:) * xb * B
%       B  = pinv(HRF(k,:)*Xb) * x
%
% This accounts for different time resolutions between
% our hemodynamic signal and the discrete representation of
% the psychological variable. In this case k is a vector 
% representing the time resolution of the scans.
%
% Version 1.3 and beyond of this routine was renamed spm_ped_ppi
% as it now uses parametric empirical Bayes to estimate the parameters
% producing the deconvolution. See spm_PEB for details of the inputs
% and outputs.
%
% For event related data one then takes the stick functions,
% convolves with signal and resamples like so:
%
%        P  = sf1 + (-sf2);
%        P  = P(k);
%
%	sf1 and sf2 are the stick functions for the event
%	related psycholgocal variable from the original 
%	design matrix
%
%  convolve and resample at each scan
%
%	 y  = conv(full(P*xb*B),HRF);
%	 y  = y(k);
%
%---------------------------------------------------------------------
% %W% Darren Gittleman %E%

% set up the graphical interface
%----------------------------------------------------------------------
Finter = spm_figure('GetWin','Interactive');
Fgraph = spm_figure('GetWin','Graphics');
header = get(Finter,'Name');

% check inputs and set up variables
%----------------------------------------------------------------------
if ~nargin
	swd   = spm_str_manip(spm_get(1,'SPM.mat','Select SPM.mat'),'H');
	load(fullfile(swd,'SPM.mat'))
	cd(swd)
end
RT     = SPM.xY.RT;
dt     = SPM.xBF.dt;
NT     = RT/dt;

% Ask whether to perform physiophysiologic or psychophysiologic interactions
%--------------------------------------------------------------------------
set(Finter,'name','PPI Setup')
ppiflag    = spm_input('Analysis type?',1,'b','PHI-PHI|PSY-PHI',[0 1],2);


if ~ppiflag  % physiophysiologic
    
    % Assume that user will examine interactions between 2 regions
    %----------------------------------------------------------------------
    spm_input('physiological variables:...  ',1,'d');
    P      = spm_get(2,'VOI*.mat',{'select VOIs'});
    for  i = 1:2
 	p     = load(P{i},'xY');
	xY(i) = p.xY;
    end
    Sess   = SPM.Sess(xY(1).Sess);

else        % psychophysiologic
    
    % get hemodynamic response 
    %----------------------------------------------------------------------
    spm_input('physiological variable:...  ',1,'d');
    P      = spm_get(1,'VOI*.mat',{'select VOIs'});
    p      = load(P{1},'xY');
    xY(1)  = p.xY;
    Sess   = SPM.Sess(xY(1).Sess);


    % get 'causes' or inputs U
    %----------------------------------------------------------------------
    spm_input('Psychological variable:...  ',1,'d');
    u      = length(Sess.U);
    U.name = {};
    U.u    = [];
    U.w    = [];
    for  i = 1:u
        for  j = 1:length(Sess.U(i).name)
	    str   = ['include ' Sess.U(i).name{j} '?'];
	    if spm_input(str,2,'y/n',[1 0])
		U.u             = [U.u Sess.U(i).u(33:end,j)];
		U.name{end + 1} = Sess.U(i).name{j};
		str             = 'Contrast weight';
		U.w             = [U.w spm_input(str,3,'e',[],1)];
	    end
        end
    end

end % end setup


% Setup some variables based on BOLD
%-------------------------------------------------------------------------
N     = length(xY(1).u);
k     = 1:NT:N*NT;  	% microtime to scan time conversion index
for i = 1:size(xY,2)
    Y(:,i) = spm_detrend(xY(i).u);
end

% name of PPI file to be saved
%-------------------------------------------------------------------------
PPI.name = spm_input('Name of PPI','!+1','s','PPI');

% create basis functions and hrf in scan time and microtime
%-------------------------------------------------------------------------
spm('Pointer','watch')
hrf    = spm_hrf(dt);

% get neuronal responses to known inputs
%-------------------------------------------------------------------------
xu    = [];
Hxu   = [];
for i = 1:length(Sess.U)
	for j = 1:size(Sess.U(i).u,2)
		u    = full(Sess.U(i).u(33:end,j));
		xu   = [xu u];

	end
end
M     = size(xu,2);

% create convolved explanatory variables.in scan time
%-------------------------------------------------------------------------
xb    = spm_dctmtx(N*NT + 128,N);
Hxb   = zeros(N,N);
Hxu   = zeros(N,M);
for i = 1:N
	Hx       = conv(xb(:,i),hrf);
	Hxb(:,i) = Hx(k + 128);
end
for i = 1:M
	Hx       = conv(xu(:,i),hrf);
	Hxu(:,i) = Hx(k);
end
xb    = xb(129:end,:);

% add filter confounds.(in scan time) and constant term
%-------------------------------------------------------------------------
try
	KH = SPM.xX.K(xY.Sess).KH;
catch
	KH = [];
end
KH     = [ones(N,1) KH];
K      = size(KH,2);

% assume neuronal response is white
%-------------------------------------------------------------------------
C      = {	sparse([1:N],[1:N],1,M + N + K,M + N + K),...
		sparse([1:M] + N,[1:M] + N,1,M + N + K,M + N + K),...
		sparse([1:K] + N + M,[1:K] + N + M,1,M + N + K,M + N + K)};

R      = Y - [Hxu KH]*(pinv(full([Hxu KH]))*Y);
h      = var(Y);

% create structure for spm_PEB
%-------------------------------------------------------------------------
P{1}.X = [Hxb Hxu KH];		% Design matrix for lowest level
P{1}.C = {speye(N,N)};		% i.i.d assumptions
P{2}.X = sparse(N + M + K,1);	% Design matrix for parameters (0's)
P{2}.C = C;

if ~ppiflag

    % COMPUTE PHYSIOPHYSIOLOGIC INTERACTIONS
    %---------------------------------------------------------------------
    C       = spm_PEB(Y(:,1),P);
    xn1     = [xb xu]*C{2}.E(1:(N + M));
    C       = spm_PEB(Y(:,2),P);
    xn2     = [xb xu]*C{2}.E(1:(N + M));
    xn1     = spm_detrend(xn1);
    xn2     = spm_detrend(xn2);
    xnxn    = xn1.*xn2;
    
    % convolve and resample at each scan for bold signal
    %---------------------------------------------------------------------
    ppi     = conv(xnxn,hrf);
    ppi     = ppi(k);
    
    % save (empty) psychological variables
    %---------------------------------------------------------------------
    PPI.psy = [];
    PPI.P   = [];
    PPI.xn  = [xn1 xn2];

else

    % COMPUTE PSYCHOPHYSIOLOGIC INTERACTIONS
    % use basis set in microtime
    %---------------------------------------------------------------------
    % get parameter estimates and neural signal; beta (C) is in scan time
    % This clever trick allows us to compute the betas in scan time which is
    % much quicker than with the large microtime vectors. Then the betas
    % are applied to a microtime basis set generating the correct neural
    % activity to convolve with the psychological variable in mircrotime
    %---------------------------------------------------------------------
    C       = spm_PEB(Y,P);
    xn      = [xb xu]*C{2}.E(1:(N + M));
    xn      = spm_detrend(xn);
    
    % setup psychological variable from inputs and contast weights
    %---------------------------------------------------------------------
    PSY     = zeros(N*NT,1);
    for i = 1:size(U.u,2)
        PSY = PSY + full(U.u(:,i)*U.w(:,i));
    end
    PSY     = spm_detrend(PSY);

    % multiply psychological variable by neural signal
    %----------------------------------------------------------------------
    PSYxn   = PSY.*xn;
        
    % convolve and resample at each scan for bold signal
    %----------------------------------------------------------------------
    ppi	    = conv(PSYxn,hrf);
    ppi     = ppi(k);

    % similarly for psychological effect
    %----------------------------------------------------------------------
    PSYHRF  = conv(PSY,hrf);
    PSYHRF  = PSYHRF(k);
    
    % save psychological variables
    %----------------------------------------------------------------------
    PPI.psy = U;
    PPI.P   = PSYHRF;
    PPI.xn  = xn;
    
end % end computation of physiophysiologic vs. psychophysiologic interactions

% setup other output variables
%--------------------------------------------------------------------------
PPI.ppi = spm_detrend(ppi);
PPI.xY  = xY;
PPI.dt  = dt;
str     = ['PPI_' PPI.name];
save(fullfile(SPM.swd,str),'PPI')


% Plot so the user can see the results
%--------------------------------------------------------------------------
figure(Fgraph);
subplot(2,1,1)
t       = RT*[1:N];
T       = dt*[1:(N*NT)];
plot(t,Y(:,1),T,PPI.xn(:,1))
xlabel('time (secs)')
ylabel('hemodynamic and neuronal responses')
title('Phsyiological variable')
axis tight square
grid on

subplot(2,2,4)
plot(t,PPI.ppi)
title('PPI')
xlabel('time (secs)')
axis tight square
grid on


if ~ppiflag

	subplot(2,2,3)
	plot(t,Y(:,2),T,PPI.xn(:,2))
	ylabel('hemodynamic and neuronal responses')
	title('2nd Phsyiological variable')


else
	subplot(2,2,3)
	plot(t,PPI.P,T,PSY,'--')
	title('psychological variable')
end
xlabel('time (secs)')
axis tight square
grid on


% clean up
%--------------------------------------------------------------------------
spm('Pointer','arrow')
spm('FigName',header);
fprintf('\ndone\n')
return

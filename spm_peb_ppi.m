function PPI = spm_peb_ppi(SPM,VOL,xX,xCon,xSDM,hReg)
% computes either physiophysiologic or psychophysiologic interactions
% FORMAT PPI = spm_peb_ppi(SPM,VOL,xX,xCon,xSDM,hReg);
%
% SPM    - Structure containing SPM, distribution & filtering detals
% VOL    - Structure containing details of volume analysed
% xX     - Design Matrix structure
% xCon   - Contrast definitions structure (see spm_FcUtil.m for details)
% xSDM   - Structure containing contents of SPM.mat file
% hReg   - Handle of results section XYZ registry (see spm_results_ui.m)
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
WS     = spm('WinScale');

% check inputs and set up variables
%----------------------------------------------------------------------
if nargin ~= 6
    error('Must enter input arguments as (SPM,VOL,xX,xCon,xSDM,hReg)')
    return
end
RT     = xX.RT;
dt     = xSDM.Sess{1}.U{1}.dt;
NT     = RT/dt;

% Ask whether to perform physiophysiologic or psychophysiologic interactions
%--------------------------------------------------------------------------
set(Finter,'name','PPI Setup')
ppiflag   = spm_input('Analysis type?',1,'b','PHI-PHI|PSY-PHI',[0 1],2);


if ~ppiflag  % physiophysiologic
    
    % Assume that user will examine interactions between 2 regions
    %----------------------------------------------------------------------
    spm_input('physiological variables:...  ',1,'d');
    P     = spm_get(2,'VOI*.mat',{'select VOIs'});
    for i = 1:2
 	p     = load(P{i},'xY');
	xY(i) = p.xY;
    end
    Sess  = xSDM.Sess{xY(1).Sess};

else        % psychophysiologic
    
    % get hemodynamic response 
    %----------------------------------------------------------------------
    spm_input('physiological variable:...  ',1,'d');
    P     = spm_get(1,'VOI*.mat',{'select VOIs'});
    p     = load(P{1},'xY');
    xY(1) = p.xY;
    Sess  = xSDM.Sess{xY(1).Sess};

    % set name of interactive window
    %----------------------------------------------------------------------
    
    % get 'causes' or inputs U
    %----------------------------------------------------------------------
    spm_input('Psychological variable:...  ',1,'d');
    u      = length(Sess.U);
    U.name = {};
    U.u    = [];
    U.w    = [];
    for  i = 1:u
        for  j = 1:length(Sess.U{i}.Uname)
	    str   = ['include ' Sess.U{i}.Uname{j} '?'];
	    if spm_input(str,2,'y/n',[1 0])
		U.u             = [U.u Sess.U{i}.u(33:end,j)];
		U.name{end + 1} = Sess.U{i}.Uname{j};
		str             = 'Contrast weight'
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
    Y(:,i) = spm_detrend(full(xY(i).u));
end

% name of PPI file to be saved
%-------------------------------------------------------------------------
PPI.name = spm_input('Name of PPI','!+1','s','PPI');

% create basis functions and hrf in scan time
%-------------------------------------------------------------------------
spm('Pointer','watch')
xb       = spm_dctmtx(N,N);
hrf      = spm_hrf(RT);

% create HRF convolution matrix. It is a toeplitz matrix
%-------------------------------------------------------------------------
H        = convmtx(hrf,N);
H        = H(1:N,:);

% create structure for spm_PEB
%-------------------------------------------------------------------------
P{1}.X   = [H*xb];		% Design matrix for lowest level
P{1}.C   = {speye(N,N)};	% Priors on lowest level
P{2}.X   = sparse(N,1);		% Design matrix for parameters (0's)
P{2}.C   = {speye(N,N)};	% Assume neuronal response is white.

if ~ppiflag

    % COMPUTE PHYSIOPHYSIOLOGIC INTERACTIONS
    % use basis set in scantime
    % get parameter estimates and neural signal; beta is in scan time
    %---------------------------------------------------------------------
    for i = 1:2
        C       = spm_PEB(Y(:,i),P);
        xn(:,i) = xb*C{2}.E;
        xn(:,i) = spm_detrend(xn(:,i));
    end
    xnxn    = xn(:,1).*xn(:,2);
    
    % convolve and resample at each scan for bold signal
    %---------------------------------------------------------------------
    ppi     = H*xnxn;
    
    % save (empty) psychological variables
    %---------------------------------------------------------------------
    PPI.psy = [];
    PPI.P   = [];

else

    % COMPUTE PSYCHOPHYSIOLOGIC INTERACTIONS
    % use basis set in microtime
    %---------------------------------------------------------------------
    xbsf    = spm_dctmtx(N*NT,N);
    
    % get parameter estimates and neural signal; beta (C) is in scan time
    % This clever trick allows us to compute the betas in scan time which is
    % much quicker than with the large microtime vectors. Then the betas
    % are applied to a microtime basis set generating the correct neural
    % activity to convolve with the psychological variable in mircrotime
    %---------------------------------------------------------------------
    C       = spm_PEB(Y,P);
    xn      = xbsf*C{2}.E;
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
    
    % generate an hrf in microtime
    %----------------------------------------------------------------------
    hrfmt   = spm_hrf(dt);
    
    % convolve and resample at each scan for bold signal
    %----------------------------------------------------------------------
    ppi	    = conv(PSYxn,hrf);
    ppi     = ppi(k);
    xn      = xn(k);

    % similarly for psychological effect
    %----------------------------------------------------------------------
    PSYHRF  = conv(PSY,hrfmt);
    PSYHRF  = PSYHRF(k);
    
    % save psychological variables
    %----------------------------------------------------------------------
    PPI.psy = U;
    PPI.P   = PSYHRF;
    
end % end computation of physiophysiologic vs. psychophysiologic interactions

% setup other output variables
%--------------------------------------------------------------------------
PPI.ppi = spm_detrend(ppi);
PPI.xY  = xY;
PPI.xn  = xn;
str     = ['PPI_' PPI.name];
save(fullfile(SPM.swd,str),'PPI')


% Plot little graphs so the user can see the results
%--------------------------------------------------------------------------
figure(Fgraph);
subplot(2,2,3)
x       = RT*[1:N];
plot(x,Y(:,1),x,xn(:,1))
xlabel('time (secs)')
ylabel('hemodynamic and neuronal responses')
title('Phsyiological variable')
axis tight square
grid on

subplot(2,2,4)
if ~ppiflag
	plot(x,Y(:,2),x,xn(:,2))
	ylabel('hemodynamic and neuronal responses')
	title('2nd Phsyiological variable')
else
	plot(x,PPI.P)
	title('convolved psychological variable')
end
xlabel('time (secs)')
axis tight square
grid on

% clean up
%--------------------------------------------------------------------------
spm('Pointer','arrow')
spm('FigName',header);
fprintf('\nDONE\n')
return

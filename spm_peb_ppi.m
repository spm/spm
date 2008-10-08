function PPI = spm_peb_ppi(SPM)
% Bold deconvolution to create physio- or psycho-physiologic interactions
% FORMAT PPI = spm_peb_ppi(SPM);
%
% SPM    - Structure containing generic details about the analysis
%
% PPI.ppi    = (PSY*xn  or xn1*xn2) convolved with the HRF
% PPI.Y      = Original BOLD eigenvariate. Use as covariate of no interest.
% PPI.P      = PSY convolved with HRF for psychophysiologic interactions,
%              or in the case of physiophysologic interactions contains
%              the eigenvariate of the second region. 
% PPI.name   = Name of PPI
% PPI.xY     = Original VOI information
% PPI.xn     = Deconvolved neural signal(s)
% PPI.U.u    = Psychological variable or input function (PPIs only)
% PPI.U.w    = Contrast weights for psychological variable (PPIs only)
% PPI.U.name = Names of psychological conditions (PPIs only)
%---------------------------------------------------------------------
%
% This routine is effectively a hemodynamic deconvolution using 
% full priors and EM to deconvolve the HRF from a hemodynamic
% time series to give a neuronal time series [that can be found in
% PPI.xn].  This deconvolution conforms to Weiner filtering 
% The neuronal process is then used to form PPIs.....
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
% Using the constraint that xn should have a uniform spectral density 
% we can expand x in terms of a discrete cosine set (xb)
%
%      xn  = xb*B
%       B  = parameter estimate
%
% The estimator of x is then
%
%       x  = HRF(k,:)*xn
%       x  = HRF(k,:) * xb * B
%
% This accounts for different time resolutions between
% our hemodynamic signal and the discrete representation of
% the psychological variable. In this case k is a vector 
% representing the time resolution of the scans.
%
% Conditional estimates of B allow for priors that ensure
% uniform variance over frequencies.
%
%---------------------------------------------------------------------
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% Darren Gitelman
% $Id: spm_peb_ppi.m 2317 2008-10-08 17:31:06Z Darren $


% set up the graphical interface
%----------------------------------------------------------------------
Finter = spm_figure('GetWin','Interactive');
spm_clf(Finter)
Fgraph = spm_figure;
spm_clf(Fgraph)
header = get(Finter,'Name');


% check inputs and set up variables
%----------------------------------------------------------------------
if ~nargin
    [P, sts] = spm_select(1,'^SPM\.mat$','Select SPM.mat');
    if ~sts, return; end
    swd   = spm_str_manip(P,'H');
    load(fullfile(swd,'SPM.mat'))
    cd(swd)
    clear P
end
RT     = SPM.xY.RT;
dt     = SPM.xBF.dt;
NT     = RT/dt;


% Ask whether to perform physiophysiologic or psychophysiologic interactions
%--------------------------------------------------------------------------
set(Finter,'name','PPI Setup')
ppiflag    = {  'simple deconvolution',...
        'psychophysiologic interaction',...
        'physiophysiologic interaction'};
i          = spm_input('Analysis type?',1,'m',ppiflag);
ppiflag    = ppiflag{i};


switch ppiflag
    
    
case  'simple deconvolution'
    %=====================================================================
    spm_input('physiological variable:...  ',2,'d');
    voi      = spm_select(1,'^VOI.*\.mat$',{'select VOI'});
    p      = load(deblank(voi(:))','xY');
    xY(1)  = p.xY;
    Sess   = SPM.Sess(xY(1).Sess);

    
    
case  'physiophysiologic interaction' % interactions between 2 regions
    %=====================================================================
    spm_input('physiological variables:...  ',2,'d');
    voi      = spm_select(2,'^VOI.*\.mat$',{'select VOIs'});
    for  i = 1:2
        p      = load(deblank(voi(i,:)),'xY');
        xY(i)  = p.xY;
    end
    Sess   = SPM.Sess(xY(1).Sess);

    
    
case  'psychophysiologic interaction'  % get hemodynamic response 
    %=====================================================================
    spm_input('physiological variable:...  ',2,'d');
    voi      = spm_select(1,'^VOI.*\.mat$',{'select VOI'});
    p      = load(deblank(voi(:))','xY');
    xY(1)  = p.xY;
    Sess   = SPM.Sess(xY(1).Sess);
    
    % get 'causes' or inputs U
    %----------------------------------------------------------------------
    spm_input('Psychological variable:...  ',2,'d');
    u      = length(Sess.U);
    U.name = {};
    U.u    = [];
    U.w    = [];
    for  i = 1:u
        for  j = 1:length(Sess.U(i).name)
            str   = ['include ' Sess.U(i).name{j} '?'];
            if spm_input(str,3,'y/n',[1 0])
                U.u             = [U.u Sess.U(i).u(33:end,j)];
                U.name{end + 1} = Sess.U(i).name{j};
                str             = 'Contrast weight';
                U.w             = [U.w spm_input(str,4,'e',[],1)];
            end
        end
    end
    
    
end % (switch setup)


% name of PPI file to be saved
%-------------------------------------------------------------------------
PPI.name    = spm_input('Name of PPI',3,'s','PPI');


% Setup variables
%-------------------------------------------------------------------------
N     = length(xY(1).u);
k     = 1:NT:N*NT;              % microtime to scan time indices


% create basis functions and hrf in scan time and microtime
%-------------------------------------------------------------------------
spm('Pointer','watch')
hrf   = spm_hrf(dt);


% create convolved explanatory {Hxb} variables in scan time
%-------------------------------------------------------------------------
xb    = spm_dctmtx(N*NT + 128,N);
Hxb   = zeros(N,N);
for i = 1:N
    Hx       = conv(xb(:,i),hrf);
    Hxb(:,i) = Hx(k + 128);
end
xb    = xb(129:end,:);


% get confounds (in scan time) and constant term
%-------------------------------------------------------------------------
X0    = xY(1).X0;
M     = size(X0,2);


% get response variable,
%-------------------------------------------------------------------------
for i = 1:size(xY,2)
    Y(:,i) = xY(i).u;
end


% remove confounds and save Y in ouput structure
%-------------------------------------------------------------------------
Yc    = Y - X0*inv(X0'*X0)*X0'*Y;
PPI.Y = Yc(:,1);
if size(Y,2) == 2
    PPI.P  = Yc(:,2);
end


% specify covariance components; assume neuronal response is white
% treating confounds as fixed effects
%-------------------------------------------------------------------------
Q      = speye(N,N)*N/trace(Hxb'*Hxb);
Q      = blkdiag(Q, speye(M,M)*1e6  );

% get whitening matrix (NB: confounds have already been whitened)
%-------------------------------------------------------------------------
W      = SPM.xX.W(Sess.row,Sess.row);

% create structure for spm_PEB
%-------------------------------------------------------------------------
P{1}.X = [W*Hxb X0];        % Design matrix for lowest level
P{1}.C = speye(N,N)/4;      % i.i.d assumptions
P{2}.X = sparse(N + M,1);   % Design matrix for parameters (0's)
P{2}.C = Q;


switch ppiflag
    
    
case  'simple deconvolution'
    %=====================================================================
    C       = spm_PEB(Y,P);
    xn      = xb*C{2}.E(1:N);
    xn      = spm_detrend(xn);
    
    % save variables
    %---------------------------------------------------------------------
    PPI.xn  = xn;
    
    % Plot so the user can see the results
    %---------------------------------------------------------------------
    figure(Fgraph);
    t       = RT*[1:N];
    T       = dt*[1:(N*NT)];
    
    subplot(2,1,1)
    plot(t,Yc,T,PPI.xn)
    title('hemodynamic and neuronal responses')
    xlabel('time (secs)')
    axis tight square
    grid on
    legend('BOLD','neuronal')
    
case  'physiophysiologic interaction' % PHYSIOPHYSIOLOGIC INTERACTIONS
    %=====================================================================
    C       = spm_PEB(Y(:,1),P);
    xn1     = xb*C{2}.E(1:N);
    C       = spm_PEB(Y(:,2),P);
    xn2     = xb*C{2}.E(1:N);
    xn1     = spm_detrend(xn1);
    xn2     = spm_detrend(xn2);
    xnxn    = xn1.*xn2;
    
    % convolve and resample at each scan for bold signal
    %---------------------------------------------------------------------
    ppi     = conv(xnxn,hrf);
    ppi     = ppi(k);
    
    % save variables
    %---------------------------------------------------------------------
    PPI.xn  = [xn1 xn2];
    PPI.ppi = spm_detrend(ppi);
    
    
    % Plot so the user can see the results
    %---------------------------------------------------------------------
    figure(Fgraph);
    t       = RT*[1:N];
    T       = dt*[1:(N*NT)];
    
    subplot(2,1,1)
    plot(t,PPI.ppi)
    title('PPI')
    xlabel('time (secs)')
    axis tight square
    grid on
    
    subplot(2,2,3)
    plot(t,Yc(:,1),T,PPI.xn(:,1))
    title('hemodynamic and neuronal responses (1st)')
    xlabel('time (secs)')
    axis tight square
    grid on
    legend('BOLD','neuronal')

    
    subplot(2,2,4)
    plot(t,Yc(:,2),T,PPI.xn(:,2))
    title('hemodynamic and neuronal responses (2nd)')
    xlabel('time (secs)')
    axis tight square
    grid on
    legend('BOLD','neuronal')
    
    
    
case  'psychophysiologic interaction'  
    %=====================================================================
    
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
    xn      = xb*C{2}.E(1:N);
    xn      = spm_detrend(xn);
    
    % setup psychological variable from inputs and contast weights
    %---------------------------------------------------------------------
    PSY     = zeros(N*NT,1);
    for i = 1:size(U.u,2)
        PSY = PSY + full(U.u(:,i)*U.w(:,i));
    end
    PSY     = spm_detrend(PSY);
    
    % multiply psychological variable by neural signal
    %---------------------------------------------------------------------
    PSYxn   = PSY.*xn;
    
    % convolve and resample at each scan for bold signal
    %---------------------------------------------------------------------
    ppi     = conv(PSYxn,hrf);
    ppi     = ppi(k);
    
    % similarly for psychological effect
    %---------------------------------------------------------------------
    PSYHRF  = conv(PSY,hrf);
    PSYHRF  = PSYHRF(k);
    
    % save psychological variables
    %---------------------------------------------------------------------
    PPI.psy = U;
    PPI.P   = PSYHRF;
    PPI.xn  = xn;
    PPI.ppi = spm_detrend(ppi);
    
    
    % Plot so the user can see the results
    %---------------------------------------------------------------------
    figure(Fgraph);
    t       = RT*[1:N];
    T       = dt*[1:(N*NT)];
    
    subplot(2,1,1)
    plot(t,Yc(:,1),T,PPI.xn(:,1))
    title('hemodynamic and neuronal responses')
    xlabel('time (secs)')
    axis tight square
    grid on
    legend('BOLD','neuronal')
    
    subplot(2,2,3)
    plot(t,PPI.P,T,PSY,'--')
    title('[convolved] psych. variable')
    xlabel('time (secs)')
    axis tight square
    grid on
    
    subplot(2,2,4)
    plot(t,PPI.ppi)
    title('PPI')
    xlabel('time (secs)')
    axis tight square
    grid on
    
    
end % (switch)

% setup other output variables
%-------------------------------------------------------------------------
PPI.xY  = xY;
PPI.dt  = dt;
str     = ['PPI_' PPI.name];

if spm_matlab_version_chk('7') >= 0,
    save(fullfile(SPM.swd,str),'-V6','PPI')
else
    save(fullfile(SPM.swd,str),'PPI')
end

% clean up
%-------------------------------------------------------------------------
spm('Pointer','arrow')
spm('FigName',header);
fprintf('\ndone\n')
return

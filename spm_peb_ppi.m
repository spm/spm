function PPI = spm_peb_ppi(varargin)
% Bold deconvolution to create physio- or psycho-physiologic interactions
% FORMAT PPI = spm_peb_ppi(SPMname,ppiflag,VOI,Uu,ppiname,showGraphics)
%
% SPM          - Structure containing generic details about the analysis or
%                the fully qualified filename of such a structure.
% ppiflag      - Type of analysis. Must be one of:
%                  'simple deconvolution'          or 'sd'
%                  'psychophysiologic interaction' or 'ppi'
%                  'physiophysiologic interaction' or 'phipi'
% VOI          - Structure containing details about a VOI (as produced by
%                spm_regions) or the fully qualified filename of such a
%                structure. If a structure, then VOI should be of size 1x1
%                in the case of simple deconvolution, and psychophysiologic 
%                interactions) or 1x2, in the case of physiophysiologic
%                interactions. If a file name it should be 1xN or 2xN.
% Uu           - Matrix of input variables and contrast weights. This is an
%                [n x 3] matrix. The first column indexes SPM.Sess.U(i). The
%                second column indexes the name of the input or cause, see
%                SPM.Sess.U(i).name{j}. The third column is the contrast
%                weight. Unless there are parametric effects the second
%                column will generally be a 1.
% ppiname      - Basename of the PPI file to save. The saved file will be:
%                <PATH_TO_SPM.MAT>/PPI_<ppiname>.mat
% showGraphics - empty or 1 = yes, 0 = no.
%
%
% PPI.ppi    - (PSY*xn  or xn1*xn2) convolved with the HRF
% PPI.Y      - Original BOLD eigenvariate. Use as covariate of no interest
% PPI.P      - PSY convolved with HRF for psychophysiologic interactions,
%              or in the case of physiophysologic interactions contains
%              the eigenvariate of the second region. 
% PPI.name   - Name of PPI
% PPI.xY     - Original VOI information
% PPI.xn     - Deconvolved neural signal(s)
% PPI.U.u    - Psychological variable or input function (PPIs only)
% PPI.U.w    - Contrast weights for psychological variable (PPIs only)
% PPI.U.name - Names of psychological conditions (PPIs only)
%__________________________________________________________________________
%
% This routine is effectively a hemodynamic deconvolution using full priors
% and EM to deconvolve the HRF from a hemodynamic time series to give a 
% neuronal time series [that can be found in PPI.xn].  This deconvolution 
% conforms to Wiener filtering. The neuronal process is then used to form 
% PPIs. See help text within function for more details.
%__________________________________________________________________________
% Copyright (C) 2002-2011 Wellcome Trust Centre for Neuroimaging

% Darren Gitelman
% $Id: spm_peb_ppi.m 4306 2011-04-13 17:02:00Z guillaume $

% SETTING UP A PPI THAT ACCOUNTS FOR THE HRF
% =========================================================================
% PPI's were initially conceived as a means of identifying regions whose
% reponses can be explained in terms of an interaction between activity in
% a specified source (the physiological factor) and some experimental
% effect (the psychological factor). However, a problem in setting up PPI's
% is that in order to derive a proper estimate of the interaction between
% a psychological variable (P) and measured hemodynamic signal (x), one 
% cannot simply convolve the psychological variable with the hrf (HRF) and 
% multiply by the signal. Thus:
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
% from x to get xn. This can then be multiplied by P and the resulting
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
% This accounts for different time resolutions between our hemodynamic 
% signal and the discrete representation of the psychological variable. In 
% this case k is a vector representing the time resolution of the scans.
%
% Conditional estimates of B allow for priors that ensure uniform variance 
% over frequencies.
%
% PPI STATISTICAL MODEL
% =========================================================================
% Once the PPI.ppi interaction term has been calculated a new GLM must be
% setup to search for the interaction effects across the brain. This is
% done using a standard, first level, fMRI model, which must include 3
% covariates, PPI.ppi (interaction), PPI.Y (main effect: source region bold
% signal) and PPI.P (main effect: "psychological" condition), plus any
% nuisance regressors according to the particular design.
%
% NB: Designs that include only the interaction term without the main
% effects are not proper as inferences on the interaction will include a
% mixture of both main and interaction effects. 
%
% Once the model has been setup and run, a contrast of [1 0 0 ] over the
% PPI.ppi, PPI.Y and PPI.P columns respectively, will show regions with a
% positive relationship to the interaction term, discounting any main
% effects. Negative regressions can be examined with [-1 0 0]. A PPI random
% effects analysis would involve taking the con*.img files from the [1 0 0]
% t-contrast for each subject and forwarding them to a second level
% analysis.


% Set up the graphical interface
%--------------------------------------------------------------------------
Finter = spm_figure('GetWin','Interactive');
header = get(Finter,'Name');
spm_clf(Finter); set(Finter,'name','PPI Setup');

% Check inputs
%--------------------------------------------------------------------------
if nargin && isstruct(varargin{1})
    SPM = varargin{1};
    try
        swd = SPM.pwd;
    catch
        swd = '';
    end
else
    try
        P            = varargin{1};
    catch
        [P, sts]     = spm_select(1,'^SPM\.mat$','Select SPM.mat');
        if ~sts, PPI = struct([]); return; end
    end
    swd = spm_str_manip(P,'H');
    load(fullfile(swd,'SPM.mat'));
end
if isempty(swd) || strcmp(swd,'.'), swd = pwd; end
SPM.swd = swd;
cwd     = pwd;
cd(SPM.swd)

% Setup variables
%--------------------------------------------------------------------------
RT      = SPM.xY.RT;
dt      = SPM.xBF.dt;
NT      = RT/dt;
fMRI_T0 = SPM.xBF.T0;

% Ask whether to perform physiophysiologic or psychophysiologic interactions
%--------------------------------------------------------------------------
try
    ppiflag = varargin{2};
catch
    ppiflag = {'simple deconvolution',...
               'psychophysiologic interaction',...
               'physiophysiologic interaction'};
    i       = spm_input('Analysis type?',1,'m',ppiflag);
    ppiflag = ppiflag{i};
end


switch lower(ppiflag)
    
    case  {'simple deconvolution','sd'}
    %======================================================================
    if nargin>2 && isstruct(varargin{3})
        p.xY = varargin{3};
    else
        try
            VOI = varargin{3};
            p   = load(deblank(VOI(1,:)),'xY');
        catch
            spm_input('physiological variable:...  ',2,'d');
            voi = spm_select(1,'^VOI.*\.mat$',{'select VOI'});
            p   = load(deblank(voi),'xY');
        end
    end
    xY(1) = p.xY;
    Sess  = SPM.Sess(xY(1).Sess);

    case  {'physiophysiologic interaction','phipi'}
    %======================================================================
    if nargin>2 && isstruct(varargin{3})
        xY = varargin{3};
        xY = xY(:)';
        if size(xY) ~= [1 2]
            error('Must include 2 VOI structures for physiophysiologic interactions')
        end
    else
        try
            VOI = varargin{3};
            if size(VOI,1) ~= 2
                error('Must include 2 VOI filenames for physiophygiologic interactions')
            end
            for i = 1:2
                p     = load(deblank(VOI(i,:)),'xY');
                xY(i) = p.xY;
            end
        catch
            spm_input('physiological variables:...  ',2,'d');
            voi      = spm_select(2,'^VOI.*\.mat$',{'select VOIs'});
            for  i = 1:2
                p      = load(deblank(voi(i,:)),'xY');
                xY(i)  = p.xY;
            end
        end
    end
    Sess = SPM.Sess(xY(1).Sess);

    case  {'psychophysiologic interaction','ppi'}
    %======================================================================
    if nargin>2 && isstruct(varargin{3})
        p.xY = varargin{3};
    else
        try
            VOI = varargin{3};
            p   = load(deblank(VOI(1,:)),'xY');
        catch
            spm_input('physiological variable:...  ',2,'d');
            voi = spm_select(1,'^VOI.*\.mat$',{'select VOI'});
            p   = load(deblank(voi(:))','xY');
        end
    end
    xY(1) = p.xY;
    Sess  = SPM.Sess(xY(1).Sess);

    % get 'causes' or inputs U
    %----------------------------------------------------------------------
    U.name = {};
    U.u    = [];
    U.w    = [];
    try
        Uu = varargin{4};
        for i = 1:size(Uu,1)
            U.u           = [U.u Sess.U(Uu(i,1)).u(33:end,Uu(i,2))];
            U.name{end+1} = Sess.U(Uu(i,1)).name{Uu(i,2)};
            U.w           = [U.w Uu(i,3)];
        end
    catch
        spm_input('Psychological variable:...  ',2,'d');
        u      = length(Sess.U);
        for  i = 1:u
            for  j = 1:length(Sess.U(i).name)
                str   = ['include ' Sess.U(i).name{j} '?'];
                if spm_input(str,3,'y/n',[1 0])
                    str             = 'Contrast weight';
                    tmpw            = spm_input(str,4,'e',[],1);
                    % if tmpw==0 then don't include the column in the
                    % design. This takes care of the possibility that
                    % the user would select to include the column but
                    % then give it a 0 weight.
                    %------------------------------------------------
                    if tmpw ~= 0
                        U.w             = [U.w tmpw];
                        U.u             = [U.u Sess.U(i).u(33:end,j)];
                        U.name{end + 1} = Sess.U(i).name{j};
                    end
                end
            end
        end
    end

end % (switch setup)


% Name of PPI file to be saved
%--------------------------------------------------------------------------
try
    PPI.name = varargin{5};
catch
    PPI.name = spm_input('Name of PPI',3,'s','PPI');
end
[tmp ppiFilename] = fileparts(PPI.name);


% Check if Graphical output should be shown
%--------------------------------------------------------------------------
try
    showGraphics = varargin{6};
catch
    showGraphics = 1;
end
if showGraphics
    Fgraph       = spm_figure('GetWin','PPI');
    spm_clf(Fgraph);
    FS           = spm('FontSizes');
end


% Setup more variables
%--------------------------------------------------------------------------
N = length(xY(1).u);
k = 1:NT:N*NT;                             % microtime to scan time indices


% Create basis functions and hrf in scan time and microtime
%--------------------------------------------------------------------------
spm('Pointer','watch')
hrf = spm_hrf(dt);


% Create convolved explanatory {Hxb} variables in scan time
%--------------------------------------------------------------------------
xb  = spm_dctmtx(N*NT + 128,N);
Hxb = zeros(N,N);
for i = 1:N
    Hx       = conv(xb(:,i),hrf);
    Hxb(:,i) = Hx(k + 128);
end
xb = xb(129:end,:);


% Get confounds (in scan time) and constant term
%--------------------------------------------------------------------------
X0 = xY(1).X0;
M  = size(X0,2);


% Get response variable,
%--------------------------------------------------------------------------
for i = 1:size(xY,2)
    Y(:,i) = xY(i).u;
end


% Remove confounds and save Y in ouput structure
%--------------------------------------------------------------------------
Yc    = Y - X0*inv(X0'*X0)*X0'*Y;
PPI.Y = Yc(:,1);
if size(Y,2) == 2
    PPI.P = Yc(:,2);
end


% Specify covariance components; assume neuronal response is white
% treating confounds as fixed effects
%--------------------------------------------------------------------------
Q = speye(N,N)*N/trace(Hxb'*Hxb);
Q = blkdiag(Q, speye(M,M)*1e6  );


% Get whitening matrix (NB: confounds have already been whitened)
%--------------------------------------------------------------------------
W = SPM.xX.W(Sess.row,Sess.row);


% Create structure for spm_PEB
%--------------------------------------------------------------------------
clear P
P{1}.X = [W*Hxb X0];        % Design matrix for lowest level
P{1}.C = speye(N,N)/4;      % i.i.d assumptions
P{2}.X = sparse(N + M,1);   % Design matrix for parameters (0's)
P{2}.C = Q;


switch ppiflag

    case  {'simple deconvolution','sd'}
    %======================================================================
    C  = spm_PEB(Y,P);
    xn = xb*C{2}.E(1:N);
    xn = spm_detrend(xn);

    % Save variables (NOTE: xn is in microtime and does not account for
    % slice timing shifts). To convert to BOLD signal convolve with a hrf.
    % Use a microtime to scan time index to convert to scan time: e.g.,
    % k = 1:NT:N*NT; where NT = number of bins per TR = TR/dt or SPM.xBF.T
    % and N = number of scans in the session. Finally account for slice
    % timing effects by shifting the index accordingly. See approximately
    % lines 420 and 421 below for an example.)
    %----------------------------------------------------------------------
    PPI.xn = xn;

    % Plot results
    %----------------------------------------------------------------------
    if showGraphics
        figure(Fgraph);
        t = RT*[1:N];
        T = dt*[1:(N*NT)];

        ax = subplot(2,1,1);
        plot(t,Yc,T,PPI.xn)
        title('hemodynamic and neuronal responses')
        xlabel('time (secs)')
        axis tight square
        grid on
        legend('BOLD','neuronal')

        str = sprintf('Simple Deconvolution: %s\n',ppiFilename);
        str = [str sprintf('VOI file: %s',xY.name)];
        hAx = axes('Position',[0 0 1 1],...
                'DefaultTextFontSize',FS(8),...
                'DefaultTextInterpreter','Tex',...
                'DefaultTextVerticalAlignment','Baseline',...
                'Parent',Fgraph,...
                'Units','points',...
                'Visible','off');
        AxPos = get(hAx,'Position'); 
        set(hAx,'XLim',[0,AxPos(3)]); set(hAx,'YLim',[0,AxPos(4)])
        dy    = FS(9);
        h     = text(dy,floor(AxPos(4))-2*dy,str);
        set(h,'Parent',hAx);
    end

    case  {'physiophysiologic interaction','phipi'}
    %======================================================================
    C    = spm_PEB(Y(:,1),P);
    xn1  = xb*C{2}.E(1:N);
    C    = spm_PEB(Y(:,2),P);
    xn2  = xb*C{2}.E(1:N);
    xn1  = spm_detrend(xn1);
    xn2  = spm_detrend(xn2);
    xnxn = xn1.*xn2;

    % Convolve, convert to scan time, and account for slice timing shift
    %----------------------------------------------------------------------
    ppi = conv(xnxn,hrf);
    ppi = ppi((k-1) + fMRI_T0);

    % Save variables
    %----------------------------------------------------------------------
    PPI.xn  = [xn1 xn2];
    PPI.ppi = spm_detrend(ppi);

    % Plot results
    %----------------------------------------------------------------------
    if showGraphics
        figure(Fgraph);
        t = RT*[1:N];
        T = dt*[1:(N*NT)];

        ax = subplot(2,1,1);
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
        
        str = sprintf('Physiophysiologic Interaction: %s\n',ppiFilename);
        str = [str, sprintf('VOI File 1: %s\n',xY(1).name)];
        str = [str, sprintf('VOI File 2: %s',xY(2).name)];
        hAx = axes('Position',[0 0 1 1],...
                'DefaultTextFontSize',FS(8),...
                'DefaultTextInterpreter','Tex',...
                'DefaultTextVerticalAlignment','Baseline',...
                'Parent',Fgraph,...
                'Units','points',...
                'Visible','off');
        AxPos = get(hAx,'Position'); 
        set(hAx,'XLim',[0,AxPos(3)]); set(hAx,'YLim',[0,AxPos(4)])
        dy    = FS(9);
        h     = text(dy,floor(AxPos(4))-2*dy,str);
        set(h,'Parent',hAx);
        
    end

    case  {'psychophysiologic interaction','ppi'}
    %======================================================================

    % COMPUTE PSYCHOPHYSIOLOGIC INTERACTIONS
    % use basis set in microtime
    %----------------------------------------------------------------------
    % get parameter estimates and neural signal; beta (C) is in scan time
    % This clever trick allows us to compute the betas in scan time which is
    % much quicker than with the large microtime vectors. Then the betas
    % are applied to a microtime basis set generating the correct neural
    % activity to convolve with the psychological variable in microtime
    %----------------------------------------------------------------------
    C  = spm_PEB(Y,P);
    xn = xb*C{2}.E(1:N);
    xn = spm_detrend(xn);

    % Setup psychological variable from inputs and contast weights
    %----------------------------------------------------------------------
    PSY = zeros(N*NT,1);
    for i = 1:size(U.u,2)
        PSY = PSY + full(U.u(:,i)*U.w(:,i));
    end
    % PSY = spm_detrend(PSY);  <- removed centering of psych variable
    % prior to multiplication with xn. Based on discussion with Karl
    % and Donald McLaren. 

    % Multiply psychological variable by neural signal
    %----------------------------------------------------------------------
    PSYxn = PSY.*xn;

    % Convolve, convert to scan time, and account for slice timing shift
    %----------------------------------------------------------------------
    ppi = conv(PSYxn,hrf);
    ppi = ppi((k-1) + fMRI_T0);

    % Convolve psych effect, convert to scan time, and account for slice
    % timing shift
    %----------------------------------------------------------------------
    PSYHRF = conv(PSY,hrf);
    PSYHRF = PSYHRF((k-1) + fMRI_T0);

    % Save psychological variables
    %----------------------------------------------------------------------
    PPI.psy = U;
    PPI.P   = PSYHRF;
    PPI.xn  = xn;
    PPI.ppi = spm_detrend(ppi);

    % Plot results
    %----------------------------------------------------------------------
    if showGraphics
        figure(Fgraph);
        t = RT*[1:N];
        T = dt*[1:(N*NT)];

        str = sprintf('Psychophyiologic Interaction: %s\n',ppiFilename);
        str = [str, sprintf('VOI File: %s\n',xY(1).name)];
        str = [str, sprintf('Factors: ')];
        for i = 1:numel(U.name)
            str = [str, sprintf('%s [%0.0f]',U.name{i},U.w(i))];
            if i < numel(U.name)
                str = [str, sprintf('; ')];
            end
        end

        ax = subplot(2,1,1);
        plot(t,Yc(:,1),T,PPI.xn(:,1))
        title('hemodynamic and neuronal responses')
        xlabel('time (secs)')
        axis tight square
        grid on
        legend('BOLD','neuronal')

        subplot(2,2,3)
        plot(T,PSY,'LineStyle','--','Color',[0 .65 0]);
        hold on
        plot(t,PPI.P,'LineStyle','-','LineWidth',1,'Color','b');
        hold off
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
        
        hAx = axes('Position',[0 0 1 1],...
                'DefaultTextFontSize',FS(8),...
                'DefaultTextInterpreter','None',...
                'DefaultTextVerticalAlignment','Baseline',...
                'Parent',Fgraph,...
                'Units','points',...
                'Visible','off');
        AxPos = get(hAx,'Position'); 
        set(hAx,'XLim',[0,AxPos(3)]); set(hAx,'YLim',[0,AxPos(4)])
        dy    = FS(9);
        h     = text(dy,floor(AxPos(4))-2*dy,str);
        set(h,'Parent',hAx);
    end
    
end % (switch)

% Setup other output variables and Save
%--------------------------------------------------------------------------
PPI.xY = xY;
PPI.dt = dt;
str    = ['PPI_' PPI.name '.mat'];

if spm_check_version('matlab','7') >= 0,
    save(fullfile(SPM.swd,str),'-V6','PPI')
else
    save(fullfile(SPM.swd,str),'PPI')
end

% Clean up
%--------------------------------------------------------------------------
spm('Pointer','arrow'); set(Finter,'name',header);
fprintf('\nCompleted PPI: %s\n',ppiFilename);
cd(cwd);

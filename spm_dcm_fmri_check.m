function [DCM] = spm_dcm_fmri_check(P, varargin)
% post-hoc diagnostics for DCM (bilinear or nonlinear) of fMRI data
% FORMAT [DCM] = spm_dcm_fmri_check(DCM)
%   DCM     - DCM structure or its filename
%
% FORMAT [GCM] = spm_dcm_fmri_check(GCM)
%   GCM     - Subjects x Models cell array of DCM structures or filenames
%
% FORMAT [DCM] = spm_dcm_fmri_check(DCM, nograph, GCM)
%   DCM     - DCM structure or its filename
%   nograph - (Optional) if true, disables graphical output
%   GCM     - (Optional) full GCM array from which the DCM in P was sourced
%             for use in graphics
%
% This routine provides some diagnostics to ensure model inversion has
% converged. It plots the predicted and observed responses over all regions
% and provides the coefficient of determination - or percent variance
% explained. This should normally be above 10%. An abnormally low
% coefficient of determination is highlighted in red. Quantitatively, one
% would normally expect to see one or more extrinsic (between source)
% connections with the strength of 1/8 Hz or greater. If all the extrinsic
% posterior expectations are below this value, then this suggests a failure
% of convergence or that the data are very noisy (possibly due to using
% very small regions of interest to summarise regional responses). Finally,
% the posterior correlations among all parameters are shown in terms of a
% correlation matrix. The number of effective parameters estimated is
% reported in terms of the (KL) divergence between the posterior and
% prior densities over parameters. This is divided by the log of the
% number of observations, by appealing to the Bayesian information
% criterion. The divergence corresponds to complexity or Bayesian
% surprise. Normally, one would expect the posterior and prior to diverge
% in a non-trivial fashion.
%
% Posterior densities are shown as bars with 90% confidence intervals in
% pink. An informed model inversion would normally provide posterior
% densities with confidence intervals that are, for some connections,
% displaced from prior expectations (at or around zero).
%
% The following diagnostics are stored in the returned DCM:
%
% DCM.diagnostics(1) - Percent variance explained
% DCM.diagnostics(2) - Largest absolute parameter estimate
% DCM.diagnostics(3) - Effective number of parameters estimated
%
% This routine is compatible with DCM8, DCM10 and DCM12 files.
%__________________________________________________________________________
% Copyright (C) 2012-2013 Wellcome Trust Centre for Neuroimaging

% Karl Friston
% $Id: spm_dcm_fmri_check.m 6800 2016-05-25 09:33:15Z peter $


%-Prepare inputs
%--------------------------------------------------------------------------
if isempty(varargin)
    nograph = false;    
else
    nograph = varargin{1};
end

if length(varargin) < 2
    GCM = []; 
else
    GCM = varargin{2};
end

%-Load DCM structure
%--------------------------------------------------------------------------
if ~nargin
    [P, sts] = spm_select(1,'^(D|G)CM.*\.mat$','select DCM or GCM mat');
    if ~sts, DCM = []; return; end
end

if isstruct(P)
    DCM = P;
elseif ischar(P)
    DCM = load(P);
    if isfield(DCM,'DCM')
        DCM = DCM.DCM;
    elseif isfield(DCM,'GCM')
        DCM = DCM.GCM;
    else
        error('Unknown DCM format');
    end
elseif iscell(P)
    DCM = P;
end

%-Handle multiple DCMs (group DCM structure)
%--------------------------------------------------------------------------
if iscell(DCM)
    
    % Prepare figure
    if ~nograph      
        datacursormode off;
        f = spm_figure('GetWin','DCM diagnostics'); spm_clf;  
        colormap(jet);
        add_title('Loading...',...
                  'Position',[0,0.5,1,0.2], ...
                  'ForegroundColor',[0.7 0.7 0.7]);          
        drawnow;
    end
    
    % Call spm_dcm_fmri_check recursively to assemble diagnostics
    [stats,DCM] = get_diagnostics(DCM);

    if nograph, return; end    
    
    % Extract explained variance
    data = cellfun(@(x)x(1),stats);    
    
    % Add title
    spm_clf; subplot(10,2,1:2); axis off;
    p = get(gca,'Position');
    add_title('DCM for fMRI Diagnostics','Position',[0,p(2)+p(4)/2,1,0.06]);
    
    % Plot grid of models    
    subplot(10,2,3:2:13);
    im = plot_model_space(data, f);
    
    % Create diagnostic panel
    subplot(10,2,4:2:14); axis off;
    h = create_diagnostic_panel();           

    % Store handles and DCM in image
    h.DCM = DCM;    
    set(im,'UserData',h);    
    
    % Enable clicks
    datacursormode on;
    
    % Plot performance over empirical bayes iterations
    if isfield(DCM{1},'FEB')
        subplot(10,3,[25 28]), bar(DCM{1}.FEB(:) - DCM{1}.FEB(1))
        xlabel('iteration','FontSize',12);
        title('Free energy','FontSize',16)
        axis square

        subplot(10,3,[26 29]), bar(DCM{1}.EEB(:),'b')
        xlabel('iteration','FontSize',12)
        title('Log precision','FontSize',16)
        axis square

        subplot(10,3,[27 30]), bar(DCM{1}.HEB(:) - DCM{1}.HEB(1),'c')
        xlabel('iteration','FontSize',12)
        title('Posterior uncertainty','FontSize',16)
        axis square
    end
        
    return;
end

% Assemble diagnostics
%==========================================================================

% coefficient of determination (percent variance explained)
%--------------------------------------------------------------------------

% Check if spectral DCM (DCM for cross spectra)
%--------------------------------------------------------------------------
try
    analysis = DCM.options.analysis;
catch
    analysis = '';
end
if strcmp(analysis,'CSD') && isfield(DCM,'Hc')
    PSS   = sum(sum(sum(abs(DCM.Hc).^2)));
    RSS   = sum(sum(sum(abs(DCM.Rc).^2)));
elseif isfield(DCM,'y')
    PSS   = sum(sum(DCM.y.^2));
    RSS   = sum(sum(DCM.R.^2));
else
    PSS   = NaN;
    RSS   = NaN;
end

D(1)  = 100*PSS/(PSS + RSS);

% largest absolute posterior expectation (extrinsic connections)
%--------------------------------------------------------------------------
try
    A = DCM.Ep.A;
catch
    A = DCM.A;
end

if isfield(DCM.options,'two_state') && DCM.options.two_state
    A = exp(A);
end

D(2)  = max(max(abs(A - diag(diag(A)))));

% complexity and effective number of parameters estimated
%--------------------------------------------------------------------------
qE    = spm_vec(DCM.Ep);
pE    = spm_vec(DCM.M.pE);
qC    = DCM.Cp;
pC    = DCM.M.pC;
k     = rank(full(pC));
pC    = pinv(pC);

D(3)  = trace(pC*qC) + (pE - qE)'*pC*(pE - qE) - spm_logdet(qC*pC) - k;
D(3)  = D(3)/log(DCM.v);

D     = full(D);
DCM.diagnostics = D;

if nograph
    return;
end

% Plot summary of inversion
%==========================================================================
spm_figure('GetWin','DCM diagnostics'); clf

% back to GCM link if needed
%--------------------------------------------------------------------------
if ~isempty(GCM)
   p = get(gcf,'Position');
    
   uicontrol('Style','PushButton','String','Return to models',...
        'Enable','Inactive','ButtonDownFcn',@backbutton_clicked, ...
        'Units','Normalized','Position',[0.04 0.96 0.2 0.03],...
        'HorizontalAlignment','Left','UserData',GCM);
end

% plot predicted and observed regional responses
%--------------------------------------------------------------------------
subplot(2,1,1);

% Check if spectral DCM (DCM for cross spectra)
%--------------------------------------------------------------------------
if strcmp(analysis,'CSD')
    
    %-CSD data
    %----------------------------------------------------------------------
    Hz   = DCM.Hz;                 % frequencies
    name = {DCM.xY.name};          % names
    ns   = size(DCM.a,1);          % number of regions
    ns   = min(ns,8);              % bounded number of regions
    
    for i = 1:ns

        Hc(:,i) = abs(DCM.Hc(:,i,i));
        Yc(:,i) = abs(DCM.Hc(:,i,i) + DCM.Rc(:,i,i));
            
    end
        
        plot(Hz,Hc), hold on
        plot(Hz,Yc,':'), hold off
        str = sprintf('variance explained %0.0f%%', D(1));
        str = {'Responses and Predictions',str};
        
        if D(1) > 10
    
            title(str,'FontSize',16);
        else
            title(str,'FontSize',16,'Color','r');
        end
        
        xlabel('frequency (Hz)')
        ylabel('abs(CSD)')
        axis square, spm_axis tight
        legend(name)
else
    t   = (1:DCM.v)*DCM.Y.dt;
    plot(t,DCM.y,t,DCM.y + DCM.R,':');
    str = sprintf('variance explained %0.0f%%', D(1));
    str = {'Responses and Predictions',str};
    
    if D(1) > 10
        title(str,'FontSize',16);
    else
        title(str,'FontSize',16,'Color','r');
    end
    
    xlabel('time {seconds}');
end

    

% posterior densities over A parameters
%--------------------------------------------------------------------------
try
    i = spm_fieldindices(DCM.Ep,'A');
catch
    i = 1 + (1:DCM.n^2);
end
qE  = spm_vec(DCM.Ep);
qC  = DCM.Cp;

if DCM.options.two_state
    qE = exp(qE);
end

subplot(2,2,3)
spm_plot_ci(qE(i),qC(i,i)), hold on
str = sprintf('largest connection strength %0.2f', D(2));
str = {'Intrinsic and Extrinsic connections',str};
if D(2) > 1/8
    title(str,'FontSize',16);
else
    title(str,'FontSize',16,'Color','r');
end
xlabel('parameters');
axis square


% posterior correlations among all parameters
%--------------------------------------------------------------------------
subplot(2,2,4)
imagesc(spm_cov2corr(DCM.Cp))
title('Posterior Correlations','FontSize',16)
str = sprintf('estimable parameters %0.0f', D(3));
str = {'Posterior Correlations',str};
if D(3) > 1
    title(str,'FontSize',16);
else
    title(str,'FontSize',16,'Color','r');
end
axis square

% =========================================================================
function [stats,GCM] = get_diagnostics(GCM)
% Builds a [s x m] cell array of diagnostic stats for an array of DCMs

stats    = cell(size(GCM));
stats(:) = {zeros(1,3)};    
[Ns, Nm] = size(GCM);    
counter = 1;
for s = 1:Ns
    for m = 1:Nm
        if isfield(GCM{s,m},'Ep')
            GCM{s,m}.v  = GCM{s,1}.v;
            GCM{s,m}    = spm_dcm_fmri_check(GCM{s,m}, true);
            stats{s,m}  = GCM{s,m}.diagnostics;
        else
            stats{s,m}  = zeros(3,1);
        end
        counter = counter + 1;
    end
end

% =========================================================================
function add_title(str,varargin)
% Draws a title
uicontrol('Style','Text','String',str,'Fontsize',20,...
          'Units','Normalized','HorizontalAlignment','Center',...
          'BackgroundColor','w',varargin{:});  

% =========================================================================
function im = plot_model_space(data, f)
% Plots the model space. Returns a graphics handle to the image

% Plot explained variances    
im = imagesc(data, [0 100]);

[Ns, Nm] = size(data); 
set(gca,'xtick', linspace(0.5,Nm+0.5,Nm+1), ...
        'ytick', linspace(0.5,Ns+.5,Ns+1));
set(gca,'xgrid', 'on', 'ygrid', 'on', ...
        'gridlinestyle', '-', 'xcolor', 'w', 'ycolor', 'w');
ylabel('Subjects','FontSize',14,'Color','k');
xlabel('Models','FontSize',14,'Color','k');
title('Variance explained (%)','FontSize',16);
axis equal

% Create small colorbar
cbar = colorbar;
p    = get(cbar,'Position');    
t    = p(2) + p(4);
p(2) = t - 0.2;
p(4) = 0.2;
set(cbar,'Position',p);        

% Handle mouseclicks
dcm_obj = datacursormode(f);
set(dcm_obj,'UpdateFcn',@gcm_clicked);      

% =========================================================================
function h = create_diagnostic_panel()
% Create a panel to show details / controls for the selected model

% Get axes top
p = get(gca,'Position');
t = p(2)+p(4);

% Add panel
he = 0.15; % Height
h = struct();
h.panel = uipanel('Units','normalized',...
'BackgroundColor',[1 1 1],'Title','Diagnostics',...
'Position',[0.5 t-he 0.4 he],'FontSize',12);    

% Add explained variance label
h.expvar = uicontrol('Style','text','Units','normalized',...
          'String','Please select a model, left.',...
          'BackgroundColor',[1 1 1],'Position',[0,0.6,1,0.2],...
          'Parent',h.panel,'HorizontalAlignment','center');

% Add button
h.infobutton = uicontrol('Style','pushbutton','Units','normalized',...
          'String','Diagnostics',...
          'Position',[0.1,0.2,0.8,0.25],...
          'Parent',h.panel,'Callback',@info_clicked,'Enable','Off',...
          'Visible','Off');    
    
% =========================================================================
function txt = gcm_clicked(varargin)
% Handle mouseclick on DCM array

% Unpack handles and explained variance matrix
im     = get(varargin{2},'Target');
h      = get(im,'UserData');
expvar = get(im,'CData');

% Get DCM for this click position
xy  = get(varargin{2},'Position');
m   = xy(1);
s   = xy(2);
GCM = h.DCM;
DCM = GCM{s,m};

% Update DCM info panel text
if isnan(expvar(s,m))
    str = 'Variance explained: unavailable';
    set(h.infobutton,'Enable','off');
else
    str = sprintf('Variance explained: %2.2f%%', expvar(s,m));
    set(h.infobutton,'Enable','on');
end
set(h.expvar, 'String', str);
set(h.panel,  'Title', sprintf('Subject %d Model %d',s,m));

% Update DCM info panel button
button_data = struct();
button_data.DCM = DCM;
button_data.GCM = GCM;
set(h.infobutton,'UserData', button_data);
set(h.infobutton,'Visible', 'on');

% Return tooltip
txt = sprintf('Subject %d model %d',s,m);

% =========================================================================
function info_clicked(varargin)
% Callback for the Diagnostics button being clicked

button_data = get(varargin{1},'UserData');
DCM = button_data.DCM;
GCM = button_data.GCM;

if ~isempty(button_data.DCM)
    datacursormode off;
    spm_dcm_fmri_check(DCM, false, GCM);
end

% =========================================================================
function backbutton_clicked(varargin)
% Callback for the back button, to return to the model space
GCM = get(varargin{1},'UserData');
clf; drawnow();
spm_dcm_fmri_check(GCM);
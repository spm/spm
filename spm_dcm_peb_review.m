function spm_dcm_peb_review(PEB, DCM)
% Review tool for DCM PEB models
% FORMAT spm_dcm_peb_review(PEB,DCM)
%
% PEB - PEB model to review
% DCM - (Optional) A single DCM or cell array of DCMs. Data is used to 
%       enhance the GUI.
%
%__________________________________________________________________________
% Copyright (C) 2016 Wellcome Trust Centre for Neuroimaging

% Peter Zeidman
% $Id: spm_dcm_peb_review.m 6946 2016-11-23 15:26:29Z peter $

% Prepare input
% -------------------------------------------------------------------------
if nargin < 1 || isempty(PEB)
    [PEB,sts] = spm_select(1,'mat','Select a PEB model',{},pwd,'^(PEB|BMA)_.*mat$');
    if ~sts, return; end
end

if nargin < 2 || isempty(DCM)
    [DCM,sts] = spm_select(1,'mat','(Optional) Select DCM or GCM file',{},pwd,'^(DCM|GCM)_.*mat$');
end

% Load / validate PEB
if ischar(PEB)
    PEB = load(PEB);
    if isfield(PEB,'PEB')
        PEB = PEB.PEB;
    elseif isfield(PEB,'BMA')
        PEB = PEB.BMA;
    else
        error('Not a valid PEB or BMA file');
    end
end

if length(PEB) > 1, PEB = PEB(1); end

if ~isfield(PEB,'Ep')
    error('Please estimate this PEB model before reviewing');
end

% Load / validate DCM
if isempty(DCM) || (iscell(DCM) && isempty(DCM{1})), DCM = {}; end
if ~iscell(DCM), DCM = {DCM}; end

if ~isempty(DCM) && ischar(DCM{1})
    DCM{1} = load(DCM{1});
    if isfield(DCM{1},'DCM')
        DCM{1} = DCM{1}.DCM;
    elseif isfield(DCM{1},'GCM')
        DCM = DCM{1}.GCM;
    else
        DCM = [];
    end
end

% Correct BMA matrix size (Ep)
np = length(PEB.Pnames); % Parameters
nc = size(PEB.M.X,2);    % Covariates
if size(PEB.Ep,2) ~= nc
    PEB.Ep = reshape(PEB.Ep,np,nc);
end

% Correct BMA matrix size (Cp)
if isvector(PEB.Cp)
    PEB.Cp = diag(PEB.Cp);
end

% Create state variable
xPEB = struct();
xPEB.PEB           = PEB;
xPEB.DCM           = DCM;
xPEB.view          = 1;     % Selected view
xPEB.sel_field_idx = 1;     % Selected first-level DCM field
xPEB.sel_input     = 1;     % Selected first-level DCM input (U)
xPEB.region_names  = {};    % First level region names
xPEB.input_names   = {};    % First level input names
xPEB.mtx_fig       = [];    % Figure handle for connectivity matrix
xPEB.threshold_idx = 1;     % P-threshold following model comparison
xPEB.threshold_method_idx=1;% Method for thresholding parameters

% Get first-level DCM metadata
if ~isempty(DCM) 
    if isfield(DCM{1},'Sname')
        % MEG / EEG
        xPEB.region_names = DCM{1}.Sname;
        xPEB.input_names  = DCM{1}.xU.name;
    elseif isfield(DCM{1},'xY')
        % fMRI
        xPEB.region_names = {DCM{1}.xY.name};    
        xPEB.input_names  = strcat({'Input '}, DCM{1}.U.name);       
    else
        % Don't use the DCM
        warning('Unknown modality');
        xPEB.DCM = [];
    end
end

% Add in posterior correlations
correlations = spm_cov2corr(PEB.Cp);
xPEB.corr    = correlations;    

% Store
assignin('base','xPEB',xPEB);

% Display
update_view();

% =========================================================================
function update_view()
% Updates the view after loading or after user requests a different view

% Read GUI state
xPEB = evalin('base','xPEB');
PEB           = xPEB.PEB;
DCM           = xPEB.DCM;
view          = xPEB.view;
threshold_idx = xPEB.threshold_idx;
threshold_method_idx = xPEB.threshold_method_idx;

% Unpack PEB metadata
np = length(PEB.Pnames); % Parameters
ns = length(PEB.Snames); % Subjects
nc = size(PEB.Ep,2);     % Covariates

% Get covariate names
if isfield(PEB,'Xnames')
    Xnames = PEB.Xnames;
else
    Xnames = cell(1,nc);
    for i = 1:nc
        Xnames{i} = sprintf('Covariate %d',i);
    end
end    

% Prepare probability threshold (Penny et al., NeuroImage, 2004)
% -------------------------------------------------------------------------
effect = view - 1;

has_Pp = isfield(PEB,'Pp') || isfield(PEB,'Pw') || isfield(PEB,'Px');
    
display_threshold = (effect > 0 && effect <= nc);

if xPEB.threshold_method_idx == 1
    % Thresholding based on posterior probability
    thresholds_str = {'No threshold (Pp >= 0)';
                  'Pp > 0.5';
                  'Pp > 0.75';
                  'Pp > 0.95';
                  'Pp > .99'};    
else
    % Thresholding based on free energy
    thresholds_str = {'No threshold (Pp >= 0)';
                  'Weak evidence (Pp > 0.5)';
                  'Positive evidence (Pp > 0.75)';
                  'Strong evidence (Pp > 0.95)';
                  'Very strong evidence (Pp > .99)'};
end
              
threshold_methods_str = {'Probability that parameter > 0';
                         'Free energy (with vs without)'};

% If there's no BMA posterior probability, disable free energy thresholding
if ~has_Pp
    xPEB.threshold_method_idx = 1;
    threshold_methods_str = threshold_methods_str(1);
end

thresholds = [0 0.5 0.75 0.95 0.99]; 
threshold  = thresholds(threshold_idx);

% Get parameters / variance for the selected covariate
% -------------------------------------------------------------------------
if effect > 0 && effect <= nc
    % Identify relevant parameters
    effect_idx         = 1:np:(np*nc);
    peb_param_idx      = effect_idx(effect) : (effect_idx(effect) + np - 1);
    xPEB.peb_param_idx = peb_param_idx;

    % Posterior means / covariance
    Ep = PEB.Ep(:,effect);
    Cp = diag(PEB.Cp);
    Cp = Cp(peb_param_idx);                
end

% Apply threshold
% -------------------------------------------------------------------------
if display_threshold
    
    % BMA from spm_dcm_peb_bmc only has probabilities for parameters which
    % differed across models. Get the original indices of these parameters.
    if isfield(PEB,'Pw') || isfield(PEB,'Px') && isfield(PEB.K)
        k = find(any(~PEB.K));        
    end
    
    if threshold_method_idx == 1
        % Threshold on posterior probability
        T  = 0;
        Pp = 1 - spm_Ncdf(T,abs(Ep),Cp);
    elseif isfield(PEB,'Pw') && effect == 1
        Pp    = zeros(size(Ep,1),1);
        Pp(k) = PEB.Pw(:);
    elseif isfield(PEB,'Px') && effect == 2
        Pp    = zeros(size(Ep,1),1);
        Pp(k) = PEB.Px;
    elseif isfield(PEB,'Pp')
        Pp = PEB.Pp(peb_param_idx);
    else
        Pp = [];
    end
    
    xPEB.Pp = Pp;
    
    if ~isempty(Pp) && threshold > 0
        Ep = Ep .* (Pp(:) > threshold);
        Cp = Cp .* (Pp(:) > threshold);
    end
else
    xPEB.Pp = [];
end

% Posterior random effects variance
Ce = PEB.Ce;

% Posterior weights on covariance components
Eh = PEB.Eh;
Ch = PEB.Ch;

% If first level DCMs are provided, unpack field names
% -------------------------------------------------------------------------

% Get names of DCM fields included in the PEB
fields = {};
parts  = {};
for p = 1:np
    [name,parts{p}] = pname_to_string(PEB.Pnames{p}, ...
                                      xPEB.region_names, ...
                                      xPEB.input_names);

    if isnan(parts{p}.input)
        parts{p}.input = 1;
    end

    if ~any(strcmp(parts{p}.field, fields))
        fields{end+1} = parts{p}.field;
    end        
end

xPEB.fields = fields;


% Reshape PEB parameters according to DCM
% -------------------------------------------------------------------------
    
sel_field_idx = xPEB.sel_field_idx - 1;

display_connectivity_selector = (effect > 0 && effect <= nc && ~isempty(DCM));
display_connectivity          = display_connectivity_selector & (sel_field_idx > 0);

sel_input = xPEB.sel_input;    

if display_connectivity            
    
    sel_field = fields{sel_field_idx};
    
    % Get the size of this field's matrix in the DCM
    [i,j,k] = size(eval(['DCM{1}.Ep.' sel_field]));
    
    % Reshape PEB parameters
    Eq = zeros(i,j,k);    
    for p = 1:np
        if strcmp(parts{p}.field, sel_field)
            Eq(parts{p}.row, parts{p}.col, parts{p}.input) = Ep(p);
        end
    end
    
    % Limit to a specific input (fMRI)
    nu = size(Eq,3);
    Eq = Eq(:,:,sel_input);
        
    xPEB.Eq     = Eq;  
end

% Set GUI constants
% -------------------------------------------------------------------------

% Drop-down menu indices
VIEW_NONE        = 1;
VIEW_COMPONENTS  = nc+2;
VIEW_DIAGNOSTICS = nc+3;

% Drop-down menu labels
views = cell(1,nc+3);
views{VIEW_NONE} = 'Please select...';
for i = 1:nc
    views{i+1} = ['Second-level effect - ' Xnames{i}];
end
views{VIEW_COMPONENTS}  = 'Precision components';
views{VIEW_DIAGNOSTICS} = 'Diagnostics';
    
% Create GUI
% -------------------------------------------------------------------------
f = spm_figure('GetWin','PEB - Review Parameters');
datacursormode off;
spm_clf;    

% Panel for design matrix
h = create_panel('Position',[0 0.70 1 0.29]);
xPEB.panels(1) = h;

% Drop-down menu for selecting effect
create_menu('Position',[0.1 0.58 0.85 0.1],'Tag','peb_select',...
    'Callback',@selected_view_changed, 'String',views,'Value',xPEB.view);

% Drop-down menu for selecting threshold
if display_threshold
    create_text('Threshold (optional):',...
           'Position',[0.1 0.585 0.65 0.05],'HorizontalAlignment','left');

    create_menu('Position',[0.28 0.54 0.33 0.1],'Tag','peb_threshold_method',...
        'Callback',@selected_threshold_method_changed, ...
        'String',threshold_methods_str,'Value',xPEB.threshold_method_idx);
    
    create_menu('Position',[0.62 0.54 0.33 0.1],'Tag','peb_threshold',...
        'Callback',@selected_threshold_changed, ...
        'String',thresholds_str,'Value',xPEB.threshold_idx);
end

% Panel for estimated parameter plot
h = create_panel('Position',[0 0.26 1 0.30]);
xPEB.panels(2) = h;

% Panel for reshaped parameters plot
h = create_panel('Position',[0.05 0 0.55 0.22]);
xPEB.panels(3) = h;

if display_connectivity_selector
    
    create_text('Display as matrix','Position',[0.13 0.8 0.6 0.10],...
                'FontSize',16,'Parent',xPEB.panels(3) );
    
    % Drop-down menu for controlling reshaped parameter plot
    create_menu('Position',[0.13 0.6 0.6 0.1],...
                'Tag','peb_select_field',...
                'Callback',@selected_field_changed,...
                'String',['Please select a field' xPEB.fields],...
                'Value',xPEB.sel_field_idx,...
                'Parent',xPEB.panels(3));
    
    % Drop-down menu for selecting first level input
    if display_connectivity && nu > 1        
        create_menu('Position',[0.13 0.45 0.6 0.1],...
                    'Tag','peb_select_input',...
                    'Callback',@selected_input_changed,...
                    'String',xPEB.input_names,...
                    'Value',xPEB.sel_input,...
                    'Parent',xPEB.panels(3));        
    end
end


% Add plots (upper panel)
% -------------------------------------------------------------------------

% Basic stats
subplot(3,7,1,'Parent',xPEB.panels(1));
create_tile(nc,{'COVARIATES'});

subplot(3,7,2,'Parent',xPEB.panels(1));
create_tile(np,{'DCM'; 'PARAMS.'});

subplot(3,7,3,'Parent',xPEB.panels(1));
create_tile(ns,'SUBJECTS');

% 2nd level design matrix
subplot(3,7,[8:10 15:17],'Parent',xPEB.panels(1));

imagesc(PEB.M.X);

set(gca,'XTick',1:nc);
xlabel('Covariate','FontSize',12); ylabel('Subject','FontSize',12);
axis square;

% Random effects variance
subplot(3,7,[4:7 11:14 18:21],'Parent',xPEB.panels(1));
imagesc(Ce);
set(gca,'Tag','rfx');
xlabel('First-Level Parameter','FontSize',12);
text(np,1,'Random effects variance','FontSize',16,...
        'Color','white','HorizontalAlignment','right',...
        'VerticalAlignment','top');
axis square;

% Add plots (lower panel) and render pop-up plots
% -------------------------------------------------------------------------
if view == VIEW_NONE
    % Welcome text    
    axes('Parent',xPEB.panels(2));
    text(0.5,0.75,'Please select which data to display, above',...
        'HorizontalAlignment','Center','FontSize',16,'Color',[0.6 0.6 0.6]);
    axis off;    
elseif view <= (nc+1)
    % Parameter plots
    
    % Add hint
    uicontrol('Style','text','Units','normalized',...
              'Position',[0.1 0.58 0.85 0.02],'String',...
              ['Tip: Click a plot below for details. (Select Tools->Data '...
              'Cursor if disabled.)'],'FontAngle','italic',...
              'BackgroundColor',[1 1 1],'FontSize',11);
    
    % Plot parameters
    axes('Parent',xPEB.panels(2));
    spm_plot_ci(Ep,Cp);
    set(gca,'Tag','parameters');
    xlabel('Parameter','FontSize',12); ylabel('Posterior','FontSize',12);
    title('Estimated Parameters','FontSize',16);
    
    % Plot connectivity matrix
    if display_connectivity
        
        f = gcf;
        
        % Get / create popup window
        if isempty(xPEB.mtx_fig) || ~ishandle(xPEB.mtx_fig)  
            
            s = get(0,'screensize');
            
            w = round(s(3) * 0.25);
            h = round(s(4) * 0.38);
            
            p = get(gcf,'OuterPosition');            
            b = p(4) - h;
                        
            l = p(1) - w;                 % Attempt left of figure window                        
            if l < 1, l = p(1)+p(3); end  % Attempt right of figure window                        
            if l > s(3), l = p(1); end    % Align with figure window
            
            xPEB.mtx_fig = figure('Name','Connectivity','NumberTitle','off',...
                                  'OuterPosition',[l b w h],'Color','white',...
                                  'ToolBar','none');
                        
        else
            figure(xPEB.mtx_fig);
        end
        
        % Set custom datacursor
        dcm_obj = datacursormode(xPEB.mtx_fig);
        set(dcm_obj,'UpdateFcn',@plot_clicked);
        datacursormode on;        
        
        % Prepare axes
        ax  = get(xPEB.mtx_fig,'CurrentAxes');
        if isempty(ax)
            ax = axes('Parent',xPEB.mtx_fig,'Units','normalized');
            p = get(ax,'Position');
            p(1) = 0.2;
            p(3) = 0.7;
            p(4) = 0.75;
            set(ax,'Position',p);
        end
        
        % Colormap 1 (hot)
        T = [1 0 0   % Red   
             1 1 0]; % Yellow                    
        x = [0 1];
        x = x(end:-1:1);
        c1 = interp1(x,T,linspace(0,1,64));
        
        % Colormap 2 (cold)
        T = [0 1 1             % Turqoise
             0 0.0745 0.6078]; % Dark blue         
        x = [0 1];
        x = x(end:-1:1);
        c2 = interp1(x,T,linspace(0,1,64));        
         
        % Combine and add white
        c = [c2;c1];
        c(64:65,:) = 1;
        
        % Plot
        imagesc(xPEB.Eq,[-1 1]);
        colorbar;
        colormap(c);
        
        % Style
        axis square;     
        xlabel('From','FontSize',12);ylabel('To','FontSize',12);
        set(gca,'XAxisLocation','top','Tag','connectivity');
        
        title_str = ['Connectivity: ' xPEB.fields{sel_field_idx}];
        if nu > 1
            title_str = [title_str ' - ' xPEB.input_names{xPEB.sel_input}];
        end
        
        title(title_str,'FontSize',16);        
        
        if size(xPEB.Eq,1) == size(xPEB.Eq,2) && ...
                size(xPEB.Eq,1) == length(xPEB.region_names)
            set(gca,'YTickLabel',xPEB.region_names,...
                'YTick',1:length(xPEB.region_names),...
                'XTickLabel',{''});            
        end
        
        set(0, 'currentfigure', f);
    end
    
elseif view == VIEW_COMPONENTS
    % Precision components    
    if ~isempty(Eh)
        subplot(1,2,1,'Parent',xPEB.panels(2));
        bar(Eh);
        xlim([0 length(Eh)+1]);
        xlabel('Precision component','FontSize',12);
        title('Precision component weights','FontSize',16);
    
        subplot(1,2,2,'Parent',xPEB.panels(2));
        imagesc(Ch);
        xlabel('Precision component','FontSize',12);
        set(gca,'XTick',1:length(Eh),'YTick',1:length(Eh));
        title('Covariance','FontSize',16);
        colorbar;
    end
elseif view == VIEW_DIAGNOSTICS
    % Plot correlations
    axes('Parent',xPEB.panels(2));
    imagesc(xPEB.corr); colorbar;
    set(gca,'Tag','correlations');
    xlabel('Parameter','FontSize',12); ylabel('Parameter','FontSize',12);
    title('Parameter Correlation','FontSize',16); axis square;    
end

% Set custom datacursor
dcm_obj = datacursormode(f);
set(dcm_obj,'UpdateFcn',@plot_clicked);
datacursormode on;

% Store view data in the workspace
assignin('base','xPEB',xPEB);

% =========================================================================
function h = create_panel(varargin)
% Creates a blank panel
h = uipanel('Units','normalized',...
    'BackgroundColor',[1 1 1],'BorderWidth',0,varargin{:});

% =========================================================================
function h = create_menu(varargin)
% Creates a dropdown menu
h = uicontrol('Style','Popupmenu','Units','normalized', ...
        'ToolTipString','','Enable','on',varargin{:});

% =========================================================================
function h = create_text(str,varargin)
% Creates a text label
h = uicontrol('Style','text','Units','normalized','String',str,...
              'BackgroundColor',[1 1 1],varargin{:});

% =========================================================================
function create_tile(stat,label)
% Creates a square display with a statistic and a label
stat = num2str(stat);

text(0.5,0.7,stat,'FontSize',18,'HorizontalAlignment','Center',...
    'Color',[1 1 1],'FontWeight','Bold');

text(0.5,0.45,label,'FontSize',8,'HorizontalAlignment','Center',...
    'VerticalAlignment','top','Color',[1 1 1],'FontWeight','Bold');

set(gca,'Color',[0.2 0.2 0.2],'XTick',[],'YTick',[]);

% =========================================================================
function [out,parts] = pname_to_string(pname, region_names, input_names)
% Translates a PEB parameter string e.g. B(1,2,3) to a friendly descriptor
% e.g. B(From region 1 to region 3 input 3)
%
% pname        - parameter name string from PEB e.g. A{2}(1,2,3)
% region_names - cell array of region names
% input_names  - cell array of input (condition) names
%
% out          - friendly name for the parameter
% parts        - cell array for row, col and (fMRI) task input

str = ['(?<field>[A-Za-z0-9\{\},]+)\('... % Match field and open bracket
       '(?<row>\d+)(,|\))'...     % Match row and open bracket or comma
       '(?<col>\d+)?(,|\))?'...   % Match column and open bracket or comma
       '(?<input>\d+)?(,|\))?'];  % Match input and open bracket or comma

parts = regexp(pname, str, 'names');

parts.row   = str2double(parts.row);
parts.col   = str2double(parts.col);
parts.input = str2double(parts.input);

if isempty(region_names)
    out = pname;
    return;
end

% Skip G-matrix for CMC model
if strcmp(parts.field,'G')
    out = pname;
    return;
end

out = [parts.field '-matrix '];
if isnan(parts.col)
    % Row only
    out = sprintf('%s %s', out, region_names{parts.row});
else
    % Row and col
    out = sprintf('%s from %s to %s', ...
        out, ...
        region_names{parts.col}, ...
        region_names{parts.row});
end

if ~isnan(parts.input)
    out = sprintf('%s (%s)', out, input_names{parts.input});
end

% =========================================================================
function selected_view_changed(varargin)
% Callback for change of view click

xPEB = evalin('base','xPEB');
xPEB.view = get(varargin{1},'Value');

assignin('base','xPEB',xPEB);

update_view();

% =========================================================================
function selected_threshold_method_changed(varargin)
% Callback for change of threshold method click
xPEB = evalin('base','xPEB');
xPEB.threshold_method_idx = get(varargin{1},'Value');

assignin('base','xPEB',xPEB);

update_view();
% =========================================================================
function selected_threshold_changed(varargin)
% Callback for change of threshold click
xPEB = evalin('base','xPEB');
xPEB.threshold_idx = get(varargin{1},'Value');

assignin('base','xPEB',xPEB);

update_view();

% =========================================================================
function selected_field_changed(varargin)
% Callback for change of DCM field click

xPEB = evalin('base','xPEB');

xPEB.sel_field_idx = get(varargin{1},'Value');
xPEB.sel_input     = 1;

assignin('base','xPEB',xPEB);

update_view();

% =========================================================================
function selected_input_changed(varargin)
% Callback for change of DCM field click

xPEB = evalin('base','xPEB');
xPEB.sel_input = get(varargin{1},'Value');

assignin('base','xPEB',xPEB);

update_view();
% =========================================================================
function txt = plot_clicked(varargin)
% Returns customised data tip for parameters in each plot

xPEB = evalin('base','xPEB');

% Selected axes
ax = get(varargin{2},'Target');
ax = get(ax,'Parent'); 
tag = get(ax,'Tag');

% Selected DCM parameter
pos = get(varargin{2},'Position');
idx1 = pos(1);
idx2 = pos(2);

% Selected 2nd level effect
effect = xPEB.view - 1;

Pnames        = xPEB.PEB.Pnames;
region_names  = xPEB.region_names;
input_names   = xPEB.input_names;

try
    Px_name = xPEB.PEB.Xnames{2};
catch
    Px_name = 'Group difference';
end

try
    Pp = xPEB.Pp;
catch
    Pp = [];
end

switch tag
    case 'parameters'
        peb_param_idx = xPEB.peb_param_idx;
        
        pname1 = pname_to_string(Pnames{idx1}, region_names, input_names);
        
        Ep = full(xPEB.PEB.Ep(idx1, effect));
        
        txt = {pname1; 
               sprintf('%3.3f', Ep);
               ' ';
               sprintf('DCM parameter %d',idx1); 
               sprintf('PEB parameter %d',peb_param_idx(idx1))};
           
        if ~isempty(Pp)
            txt = vertcat(txt,' ',sprintf('Probability: %2.2f',Pp(idx1)));
        end
        
    case 'rfx'        
        pname1 = pname_to_string(Pnames{idx1}, region_names, input_names);
        pname2 = pname_to_string(Pnames{idx2}, region_names, input_names);        
        
        Ce = full(xPEB.PEB.Ce(idx1,idx2));
        
        if idx1==idx2            
            txt = {sprintf('%s',pname1);
                   sprintf('DCM parameter %d',idx1); 
                   sprintf('Variance: %2.2f',Ce)};
        else
                        
            txt = {pname1;
                   'and';
                   pname2;
                   sprintf('DCM parameters %d and %d',idx1,idx2); 
                   sprintf('Covariance: %2.2f',Ce)};            
        end
    case 'correlations'        
        % Get index within covariate from the index across all parameters
        idx1_allparams = mod(idx1-1,length(Pnames)) + 1;
        idx2_allparams = mod(idx2-1,length(Pnames)) + 1;
        
        % Get the covariate of each parameter
        cov1 = fix((idx1-1)/length(Pnames)) + 1;
        cov2 = fix((idx2-1)/length(Pnames)) + 1;
        
        if cov1 == cov2
            cov_str = sprintf('Covariate %d',cov1);
        else
            cov_str = sprintf('Covariates %d and %d',cov1,cov2);
        end
        
        corr = full(xPEB.corr);
        
        txt = {sprintf('%s and %s',Pnames{idx1_allparams},Pnames{idx2_allparams});
               cov_str;
               sprintf('DCM parameters %d and %d',idx1_allparams,idx2_allparams); 
               sprintf('PEB parameters %d and %d',...
                idx1,...
                idx2);
                sprintf('Correlation: %2.2f',corr(idx2,idx1)); };
    case 'connectivity'   
        
        sel_field_idx = xPEB.sel_field_idx - 1;
        
        % Exclude for CMC
        if strcmp(xPEB.fields{sel_field_idx},'G')
            txt = '';
            return;
        end
                
        idx_from = idx1;
        idx_to   = idx2;
        
        r_from = region_names{idx_from};
        r_to   = region_names{idx_to};
        Eq     = xPEB.Eq;
        
        txt = {sprintf('From %s to %s',r_from,r_to);
               sprintf('%2.2f',Eq(idx_to,idx_from))};
    otherwise
        txt = '';
end 
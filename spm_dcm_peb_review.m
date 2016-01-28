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
% $Id: spm_dcm_peb_review.m 6702 2016-01-28 15:10:49Z guillaume $

% Prepare input
% -------------------------------------------------------------------------
if nargin < 1 || isempty(PEB)
    [PEB,sts] = spm_select(1,'mat','Select a PEB model',{},pwd,'^PEB_.*mat$');
    if ~sts, return; end
end

% Load / validate PEB
if ischar(PEB)
    PEB = load(PEB);
    PEB = PEB.PEB;
end

if length(PEB) > 1, PEB = PEB(1); end

if ~isfield(PEB,'Ep')
    error('Please estimate this PEB model before reviewing');
end

% Load / validate DCM
if nargin < 2 || isempty(DCM) || isempty(DCM{1}), DCM = {}; end
if ~iscell(DCM), DCM = {DCM}; end

if ~isempty(DCM) && ischar(DCM{1})
    DCM{1} = load(DCM{1});
    DCM{1} = DCM{1}.DCM;
end

if ~isempty(DCM) && iscell(DCM{1}) && numel(DCM) == 1
    DCM = DCM{1};
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

% Get parameters / variance for the selected covariate
% -------------------------------------------------------------------------
effect = view - 1;

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

% Posterior random effects variance
Ce = PEB.Ce;

% Posterior weights on covariance components
Eh = PEB.Eh;
Ch = PEB.Ch;

% If first level DCMs are provided, unpack
% -------------------------------------------------------------------------
display_connectivity = (effect > 0 && effect <= nc && ~isempty(DCM));

if display_connectivity   
    
    sel_field_idx = xPEB.sel_field_idx;
    sel_input     = xPEB.sel_input;    
    
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
    
    sel_field = fields{sel_field_idx};
    
    % Get the size of this field's matrix in the DCM
    [i,j,k] = size(eval(['DCM{1}.Ep.' sel_field]));
    
    % Reshape PEB parameters
    Eq = zeros(i,j,k);    
    for p = 1:np
        if strcmp(parts{p}.field, sel_field)
            Eq(parts{p}.row, parts{p}.col, parts{p}.input) = PEB.Ep(p,effect);
        end
    end
    
    % Limit to a specific input (fMRI)
    nu = size(Eq,3);
    Eq = Eq(:,:,xPEB.sel_input);
    
    xPEB.fields = fields;
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
f = spm_figure('GetWin','PEB - Review');
spm_figure('ColorMap','gray-jet');
datacursormode off;
spm_clf;    

% Panel for design matrix
h = create_panel('Position',[0 0.70 1 0.29]);
xPEB.panels(1) = h;

% Drop-down menu for selecting effect
uicontrol('Style','Popupmenu','Units','normalized', ...
    'Position',[0.1 0.58 0.85 0.1],'Tag','peb_select',...
    'Callback',@selected_view_changed, ... 
    'ToolTipString','','String',views,'Value',xPEB.view,'Enable','on');

% Panel for estimated parameter plot
h = create_panel('Position',[0 0.3 1 0.30]);
xPEB.panels(2) = h;

% Panel for reshaped parameters plot
h = create_panel('Position',[0.05 0 0.55 0.29]);
xPEB.panels(3) = h;

if display_connectivity
    % Drop-down menu for controlling reshaped parameter plot
    create_menu('Position',[0.55 0.13 0.35 0.1],...
                'Tag','peb_select_field',...
                'Callback',@selected_field_changed,...
                'String',xPEB.fields,...
                'Value',xPEB.sel_field_idx);
    
    % Drop-down menu for selecting first level input
    if nu > 1        
        create_menu('Position',[0.55 0.1 0.35 0.1],...
                    'Tag','peb_select_input',...
                    'Callback',@selected_input_changed,...
                    'String',xPEB.input_names,...
                    'Value',xPEB.sel_input);        
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

image(rescale(PEB.M.X, 1, 64));

set(gca,'XTick',1:nc);
xlabel('Covariate','FontSize',12); ylabel('Subject','FontSize',12);
axis square;

% Random effects variance
subplot(3,7,[4:7 11:14 18:21],'Parent',xPEB.panels(1));
image(rescale(Ce,1,16));
set(gca,'Tag','rfx');
xlabel('First-Level Parameter','FontSize',12);
text(np,1,'Random effects variance','FontSize',16,...
        'Color','white','HorizontalAlignment','right');
axis square;

% Add plots (lower panel)
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
              'Position',[0.1 0.6 0.85 0.03],'String',...
              ['Tip: Click a plot below for details. (Select Tools->Data '...
              'Cursor if disabled.)'],'FontAngle','italic',...
              'BackgroundColor',[1 1 1]);
    
    % Plot parameters
    axes('Parent',xPEB.panels(2));
    spm_plot_ci(Ep,Cp);
    set(gca,'Tag','parameters');
    xlabel('Parameter','FontSize',12); ylabel('Posterior','FontSize',12);
    title('Estimated Parameters','FontSize',16);
    
    % Plot connectivity matrix
    if display_connectivity
        axes('Parent',xPEB.panels(3));
        image(rescale(xPEB.Eq,65,128)); 
        
        axis square;     
        xlabel('From','FontSize',12);ylabel('To','FontSize',12);
        set(gca,'XAxisLocation','top','Tag','connectivity');
        title('Connectivity','FontSize',16);
        
        if size(xPEB.Eq,1) == size(xPEB.Eq,2) && ...
                size(xPEB.Eq,1) == length(xPEB.region_names)
            set(gca,'YTickLabel',xPEB.region_names,'XTickLabel',{''});
        end               
    end
    
elseif view == VIEW_COMPONENTS
    % Precision components    
    subplot(1,2,1,'Parent',xPEB.panels(2));
    bar(Eh);
    xlim([0 length(Eh)+1]);
    xlabel('Precision component','FontSize',12);
    title('Precision component weights','FontSize',16);
    
    subplot(1,2,2,'Parent',xPEB.panels(2));
    image(rescale(Ch, 1, 64));
    xlabel('Precision component','FontSize',12);
    set(gca,'XTick',1:length(Eh),'YTick',1:length(Eh));
    title('Covariance','FontSize',16);
    colorbar;
elseif view == VIEW_DIAGNOSTICS
    % Plot correlations
    axes('Parent',xPEB.panels(2));
    imagesc(xPEB.corr);
    set(gca,'Tag','correlations');
    xlabel('Parameter','FontSize',12); ylabel('Parameter','FontSize',12);
    title('Parameter Correlation','FontSize',16); axis square;    
end

% Set custom datacursor
dcm_obj = datacursormode(f);
set(dcm_obj,'UpdateFcn',@parameter_clicked);
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
function txt = parameter_clicked(varargin)
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

Pnames        = xPEB.PEB.Pnames;
region_names  = xPEB.region_names;
input_names   = xPEB.input_names;

switch tag
    case 'parameters'
        peb_param_idx = xPEB.peb_param_idx;
        
        pname1 = pname_to_string(Pnames{idx1}, region_names, input_names);
        
        txt = {pname1;
               sprintf('DCM parameter %d',idx1); 
               sprintf('PEB parameter %d',peb_param_idx(idx1))};
    case 'rfx'        
        pname1 = pname_to_string(Pnames{idx1}, region_names, input_names);
        pname2 = pname_to_string(Pnames{idx2}, region_names, input_names);
        
        if idx1==idx2            
            txt = {sprintf('%s',pname1);
                   sprintf('DCM parameter %d',idx1); 
                   sprintf('Variance: %2.2f',xPEB.PEB.Ce(idx1,idx1))};
        else
                        
            txt = {pname1;
                   'and';
                   pname2;
                   sprintf('DCM parameters %d and %d',idx1,idx2); 
                   sprintf('Covariance: %2.2f',xPEB.PEB.Ce(idx1,idx2))};            
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
        
        txt = {sprintf('%s and %s',Pnames{idx1_allparams},Pnames{idx2_allparams});
               cov_str;
               sprintf('DCM parameters %d and %d',idx1_allparams,idx2_allparams); 
               sprintf('PEB parameters %d and %d',...
                idx1,...
                idx2);
                sprintf('Correlation: %2.2f',xPEB.corr(idx1_allparams,idx2_allparams)); };
    case 'connectivity'       
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

% =========================================================================
function new_X = rescale(X, new_min, new_max)
% Rescales the values in X to fall in the range new_min to new_max

current_min = min(X(:));
current_max = max(X(:));

scale_factor = (current_max - current_min) / (new_max - new_min);
new_X = new_min + (X - current_min) / scale_factor;    
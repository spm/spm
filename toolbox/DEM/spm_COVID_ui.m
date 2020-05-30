function spm_COVID_ui
% FORMAT spm_COVID_ui
% Graphical user interface for DCM for COVID-19 (DEM_COVID)
% 
% For a guide to this interface, please use the help menu after launching
% the tool.
%__________________________________________________________________________
% Copyright (C) 2020 Wellcome Centre for Human Neuroimaging

% Peter Zeidman
% $Id: spm_COVID_ui.m 7866 2020-05-30 09:57:38Z karl $

% Add DEM toolbox to path
if ~isdeployed
    addpath(fullfile(spm('Dir'),'toolbox','DEM'));
end

% Check for the CSV files
required = {'time_series_covid19_confirmed_global.csv';
            'time_series_covid19_deaths_global.csv';
            'time_series_covid19_recovered_global.csv'};
all_found = all(cellfun(@(x)exist(x,'file'), required));

if ~all_found
    str = ['Please ensure that the data files are on the MATLAB path:\n' ...
           '* time_series_covid19_confirmed_global.csv\n' ...
           '* time_series_covid19_deaths_global.csv\n' ...
           '* time_series_covid19_recovered_global.csv\n' ...
           'They can be  downloaded from: ' ...
           '  https://github.com/CSSEGISandData/COVID-19/tree/master/' ...
           'csse_covid_19_data/csse_covid_19_time_series'];
    error(sprintf(str));
end

% Create GUI
[f,h] = construct_ui;

% Set up data and priors and display them
h = prepare_data_and_priors(h);

% Disable certain controls if compiled
if isdeployed
    set(h.menu_gen,'Enable','off');
end

% Store handles
guidata(f,h);

% -------------------------------------------------------------------------
% Creates the GUI
function [f,h] = construct_ui

m = [0.03 0.01];  % Margin size for panels (width and height)

% Work out position for the new window (left,bottom,width,height)
position = [0.05 0.05 0.3 0.8];

% Create figure
f = figure('Units','normalized','position',position,...
    'Name','DCM for COVID-19: Untitled model','NumberTitle','off',...
    'Toolbar','none','MenuBar','none','Visible','on','HandleVisibility','off');

% Create menu
h.menu_file   = uimenu(f,          'Label','&File','Tag','menu_file');
h.menu_new    = uimenu(h.menu_file,'Label','&New...','Tag','menu_new');
h.menu_open   = uimenu(h.menu_file,'Label','&Open Model...','Tag','menu_open');
h.menu_saveas = uimenu(h.menu_file,'Label','&Save As...','Tag','menu_saveas');
h.menu_quit   = uimenu(h.menu_file,'Label','&Quit','Tag','menu_quit');
h.menu_edit   = uimenu(f,          'Label','&Edit','Tag','menu_edit');
h.menu_gen    = uimenu(h.menu_edit,'Label','&Generative Model','Tag','menu_gen');
h.menu_help   = uimenu(f,          'Label','&Help','Tag','menu_help');
h.menu_guide  = uimenu(h.menu_help,'Label','&User Guide','Tag','menu_guide');

% Set panel width, heights and bottoms (ordered bottom to top)
pw = (1-(2*m(1)));
ph = [0.20 0.30 0.15 0.35] - m(2);
pb = cumsum([0 ph(1:end-1)]) + (m(2) .* [1 2 3 4]);

% Flip order of panel sizes to: top to bottom
ph = fliplr(ph);
pb = fliplr(pb);

% Add panels (top to bottom)
h.panel_priors   = create_panel(f,[m(1) pb(1) pw, ph(1)],'panel_priors',  '1. Set priors');
h.panel_invert   = create_panel(f,[m(1) pb(2) pw, ph(2)],'panel_invert',  '2. Fit models to data');
h.panel_results  = create_panel(f,[m(1) pb(3) pw, ph(3)],'panel_results', '3. View results');
h.panel_analysis = create_panel(f,[m(1) pb(4) pw, ph(4)],'panel_analysis','4. Bayesian model comparison');

% Add controls (priors panel)
h.uitable_priors  = create_table(h.panel_priors, [m(1) m(2) pw 1-(m(2)*2) ],'uitable_priors');

% Add controls (estimate panel)
h.button_estimate = create_button(h.panel_invert,[0.3 0.25 0.4 0.5],        'button_estimate','Estimate');

% Add controls (results panel)
h.edit_legend       = create_editbox(h.panel_results,   [m(1) m(2) pw 0.5],                                'edit_legend',       '');
h.label_view        = create_text(h.panel_results,      [m(1) 0.5+m(2) 0.3 0.2],                           'label_view',        'Result:');
h.popupmenu_view    = create_popupmenu(h.panel_results, [m(1)*2+0.3 0.5+m(2) 1-(m(1)*3+0.3) 0.2],          'popupmenu_view',    'Select result to view');
h.label_country     = create_text(h.panel_results,      [m(1) 0.5+m(2)+0.2+m(2) 0.3 0.2],                  'label_country',     'Default country:');
h.popupmenu_country = create_popupmenu(h.panel_results, [m(1)*2+0.3 0.5+m(2)+0.2+m(2) 1-(m(1)*3+0.3) 0.2], 'popupmenu_country', 'countries');
set(h.label_view,'HorizontalAlignment','right');
set(h.label_country,'HorizontalAlignment','right');
set(h.edit_legend','Enable','inactive');

% Add controls (Bayesian model comparison panel)
h.button_BMC     = create_button(h.panel_analysis,[0.10 0.5 0.35 0.4],'button_BMC','Selected country');
h.button_BMC_all = create_button(h.panel_analysis,[0.55 0.5 0.35 0.4],'button_BMC_all','All countries');
h.label_BMC      = create_text(h.panel_analysis,  [0.10 0.05 0.80 0.4],'label_BMC',...
    'Tip: in the file selector, hold down control (or command on Mac) to select multiple models');
set(h.label_BMC','HorizontalAlignment','left');

% Set callbacks
set(h.menu_new','Callback',         @callback_menu_new);
set(h.menu_open','Callback',        @callback_menu_open);
set(h.button_estimate','Callback',  @callback_button_estimate);
set(h.popupmenu_view','Callback',   @callback_popupmenu_view);
set(h.popupmenu_country','Callback',@callback_popupmenu_country);
set(h.button_BMC','Callback',       @callback_BMC);
set(h.button_BMC_all','Callback',   @callback_BMC_all);
set(h.menu_gen','Callback',         @callback_menu_edit);
set(h.menu_saveas','Callback',      @callback_menu_saveas);
set(h.menu_quit','Callback',        @callback_menu_quit);
set(h.menu_guide','Callback',       @callback_menu_guide);

% Disable controls that depend on results
set_results_controls_enabled(false,h);

% Store handles
guidata(f,h);

% -------------------------------------------------------------------------
% Stores the data and priors in the handles structure and resets the GUI
function h = prepare_data_and_priors(h)

% Read in the data
h = import_data(h);

% Get the default priors
[pE,pC,str] = spm_COVID_priors;
h.D.pE  = pE;
h.D.pC  = pC;
h.D.str = str;

% Update display
display_priors(h);
display_country_list(h);

% Reset menu items
set(h.popupmenu_view,'String','Select result to view');
set(h.edit_legend,'String','');

set_results_controls_enabled(false, h);

% -------------------------------------------------------------------------
% Reads in the time series data and stores in the handles structure
function handles = import_data(handles)

D = struct();

D.data    = DATA_COVID_JHU;
D.country = 'United Kingdom';

handles.D = D;

% -------------------------------------------------------------------------
% Copies edited priors from the GUI to the handles structure
function handles = update_priors_from_gui(handles)

% Read GUI
data = get(handles.uitable_priors,'Data');
pE = data(:,3); % expectation
pC = data(:,4); % variance

% Cell->matrix
pE = cell2mat(pE);
pC = cell2mat(pC);

% Log the expectations (while preventing infinity)
min_val           = 1e-6;
pE(pE < min_val)  = min_val;
pE                = log(pE);

% Matrix->structure
pE = spm_unvec(pE, handles.D.pE);
pC = spm_unvec(pC, handles.D.pC);

% Store
handles.D.pE = pE;
handles.D.pC = pC;

% -------------------------------------------------------------------------
% Displays the priors
function display_priors(handles)

if isfield(handles,'GCM')
    % individual estimation (not empirical Bayes)
    pE = handles.GCM{1}.M.pE; 
    pC = handles.GCM{1}.M.pC;
else
    % default priors
    pE = handles.D.pE; 
    pC = handles.D.pC;
end

varnames = fieldnames(pE);          % Matlab fields
labels   = handles.D.str.names(:);  % description

pE       = spm_vec(pE);    % prior expectation
if size(pC,2) > 1
    % covariance matrix -> vector of variance
    pC = diag(pC);    
else
    % structure or vector of variances -> vector
    pC = spm_vec(pC); 
end

% Un-log the prior expectation
pE = exp(pE);

% Set any very small numbers to zero for visual clarity
min_val = 1e-5;
pE(pE <= min_val) = 0;

set(handles.uitable_priors,'Data',[labels varnames num2cell([pE pC])]);
set(handles.uitable_priors,'ColumnName',{'Parameter description','Variable','Expectation','Variance'});
set(handles.uitable_priors,'ColumnEditable',[false false true true]);

% -------------------------------------------------------------------------
% Displays the list of countries
function display_country_list(handles)

country   = handles.D.country;        % default country
countries = {handles.D.data.country}; % list of countries

set(handles.popupmenu_country,'String',countries);

% Set default
idx = find(strcmp(country,countries));
if ~isempty(idx)
    set(handles.popupmenu_country,'Value',idx(1));
end

% -------------------------------------------------------------------------
% Prompts for a DCM file, loads it and stores it in handles
function [TF,handles] = load_model(hObject,handles)

[filename,pathname] = uigetfile('*.mat','Load model');

TF = ~isequal(filename,0) && ~isequal(pathname,0);
if ~TF
    return;
end

% Unpack
f = load(fullfile(pathname,filename));
handles.D.GCM = f.GCM;
handles.D.DCM = f.DCM;
handles.D.PEB = f.PEB;
handles.D.BMA = f.BMA;
handles.D.BMR = f.BMR;
handles.D.BPA = f.BPA;
handles.D.pE  = f.GCM{1}.M.pE;
handles.D.pC  = f.GCM{1}.M.pC;
handles.D.str = f.GCM{1}.M.str;

% Update figure title
fig = ancestor(hObject,'figure');
set(fig,'Name',filename);

% -------------------------------------------------------------------------
% Fits the model to the data
function estimate_model(hObject,handles)

% Read GUI priors table
handles = update_priors_from_gui(handles);

% Unpack
D    = handles.D;
pE   = D.pE;
pC   = D.pC;
str  = D.str;
data = D.data;

% Output structure - cell array of DCMs
GCM = cell(size(data(:)));
DCM = {};

for i = 1:numel(data)
    % Get data for one country
    Y = [data(i).death, data(i).cases];

    % Invert model
    [F,Ep,Cp,pE,pC] = spm_COVID(Y,pE,pC);

    % Create DCM structure
    dcm.M.pE  = pE;
    dcm.M.pC  = pC;
    dcm.M.str = str;
    dcm.Ep    = Ep;
    dcm.Cp    = Cp;
    dcm.F     = F;
    dcm.Y     = Y;    
    
    % Store
    GCM{i} = dcm;   
end

% Second level analysis
if numel(GCM) > 1
    % Build design matrix
    lat    = spm_vec([data.lat]);
    lon    = spm_vec([data.long]);
    lat    = lat*2*pi/360;
    lon    = lon*2*pi/360;
    X      = [];
    Xn     = {'const','log cell size'};
    for  i = 1:4
        X  = [X sin(i*lon) sin(i*lat)];
        Xn = {Xn{:}, sprintf('lat(%d)',i), sprintf('lon(%d)',i)};
    end
    X      = [log(spm_vec([data.pop])) X];
    X      = [ones(numel(data),1) X];
    X      = spm_orth(X,'norm');
    X(:,1) = 1;
    GLM.X      = X;
    GLM.Xnames = Xn;
    
    % Set GLM parameters
    GLM.alpha  = 1;
    GLM.beta   = 8;
    
    % Parametric empirical Bayes (with random effects in str.field)
    [PEB,DCM] = spm_dcm_peb(GCM,GLM,str.field);

    % Bayesian model averaging (over reduced models), testing for GLM effects
    [BMA,BMR] = spm_dcm_bmr_all(PEB,str.field);
    
    % Repeat country-specific inversions using empirical priors
    for i = 1:numel(DCM)
        % Variational Laplace
        [F,Ep,Cp] = spm_COVID(DCM{i}.Y,DCM{i}.M.pE,DCM{i}.M.pC);

        % Assemble prior and posterior estimates (and log evidence)
        DCM{i}.Ep = Ep;
        DCM{i}.Cp = Cp;
        DCM{i}.F  = F;
    end
    
    % Bayesian parameter averaging (over countries)
    BPA = spm_dcm_bpa(DCM,'nocd');
    
    % Store
    handles.D.PEB = PEB;
    handles.D.BMA = BMA;
    handles.D.BMR = BMR;
    handles.D.BPA = BPA;
end

% Store handles
handles.D.GCM = GCM;
handles.D.DCM = DCM;
guidata(hObject,handles);

% Save
timestamp = datestr(now);
timestamp = strrep(strrep(timestamp,':','-'),' ','_');
fn = sprintf('DCM_covid19_%s.mat',timestamp);
save(fn,'GCM','DCM','BMA','BPA','PEB','BMR');

% Update figure title
fig = ancestor(hObject,'figure');
set(fig,'Name',fn);

% Display results
set_results_controls_enabled(true,handles);
display_results(hObject, handles);

% -------------------------------------------------------------------------
% Prompts the user to select files then performs Bayesian model comparison
function compare_models(handles, country_idx)

if nargin < 2
    country_idx = [];
end

% Select file
[file,path] = uigetfile('DCM_*.mat','Select models to compare','MultiSelect','on');

% Validate
if ischar(file)
    errordlg('Please select at least two models to compare','Model comparison');
    return;
elseif iscell(file) && ~isempty(file)
    % Multiple selected - good
else
    return;
end

% Get free energy of each model
nm = size(file,2);
F  = nan(1,nm);
for i = 1:nm
    model = load(fullfile(path,file{i}));
    
    if isempty(country_idx)
        % all countries - use the PEB
        F(i) = model.PEB.F;
    else
        % specific country - use the DCM
        F(i) = model.DCM{country_idx}.F;
    end
end

% FFX posterior over models
P    = sum(F,1);
P    = P - max(P);
P    = exp(P);
post = P/sum(P);

% Set names for the models
names = strrep(file, '_', '\_');
names = strrep(names,'.mat','');

% Display free energy
spm_figure('GetWin','Results');
spm_clf;
subplot(2,3,2:3);
barh(F -min(F));
title('Log Bayes Factor','FontSize',16);
set(gca,'YTick',1:nm,'YTickLabel',names,'FontSize',12);
axis square, box off

% Display probability
subplot(2,3,5:6);
barh(post);
title('Posterior probability','FontSize',16);
set(gca,'YTick',1:nm,'YTickLabel',names,'FontSize',12);
axis square, box off
% -------------------------------------------------------------------------
% Switches on or off controls that depend on having an estimated model
function set_results_controls_enabled(TF,handles)
if TF
    str = 'on';
else
    str = 'off';
end
set(handles.popupmenu_view','Enable',str);
set(handles.menu_saveas','Enable',str);
% -------------------------------------------------------------------------
function callback_menu_new(hObject, eventdata)
% Resets the GUI with a new model

% Get handles
h = guidata(hObject);

% Reset title
fig = ancestor(hObject,'figure');
set(fig,'Name','DCM for COVID-19: Untitled model');

% Clear stored objects
if isfield(h,'D')
    h = rmfield(h,'D');
end

% Set up prequesits for modelling
h = prepare_data_and_priors(h);

guidata(hObject,h);

% -------------------------------------------------------------------------
function callback_menu_open(hObject, eventdata)
% Prompts for a filename, loads the model and displays it
handles = guidata(hObject);

% Load
[TF,handles] = load_model(hObject, handles);
if TF
    % Store and display
    guidata(hObject, handles);
    display_priors(handles);
    display_country_list(handles);
    display_results(hObject, handles);
    set_results_controls_enabled(true,handles);
end
% -------------------------------------------------------------------------
function callback_menu_saveas(hObject, eventdata)
% Prompts for a filename and saves

% Retrieve structures to save
handles = guidata(hObject);
DCM = handles.D.DCM;
GCM = handles.D.GCM;
PEB = handles.D.PEB;
BMA = handles.D.BMA;
BMR = handles.D.BMR;
BPA = handles.D.BPA;

% Create default filename
timestamp = datestr(now);
timestamp = strrep(strrep(timestamp,':','-'),' ','_');
defaultname = sprintf('DCM_covid19_%s.mat',timestamp);

% Prompt
[filename,pathname] = uiputfile('*.mat','Save model as', defaultname);

% Save
if ~isequal(filename,0) && ~isequal(pathname,0)
    save(fullfile(pathname,filename),'GCM','DCM','BMA','BPA','PEB','BMR');
end
% -------------------------------------------------------------------------
function callback_button_estimate(hObject, eventdata)
handles = guidata(hObject);
estimate_model(hObject, handles);
% -------------------------------------------------------------------------
function callback_popupmenu_country(hObject, eventdata)
handles = guidata(hObject);
if isfield(handles,'D') && isfield(handles.D,'DCM')
    display_results(hObject, handles);
end
% -------------------------------------------------------------------------
function callback_popupmenu_view(hObject, eventdata)
handles = guidata(hObject);
display_results(hObject, handles);
% -------------------------------------------------------------------------
function callback_menu_edit(hObject, eventdata)
open('spm_COVID_gen.m');
% -------------------------------------------------------------------------
function callback_BMC(hObject, eventdata)
handles     = guidata(hObject);
country_idx = get(handles.popupmenu_country,'Value');
compare_models(handles,country_idx);
% -------------------------------------------------------------------------
function callback_BMC_all(hObject, eventdata)
handles = guidata(hObject);
compare_models(handles);
% -------------------------------------------------------------------------
function callback_menu_quit(hObject, eventdata)
fig = ancestor(hObject,'figure');
close(fig);
% -------------------------------------------------------------------------
function callback_menu_guide(hObject, eventdata)
web('https://www.fil.ion.ucl.ac.uk/~pzeidman/covid19/guide.html');

% -------------------------------------------------------------------------
% Creates and styles a uipanel
function h = create_panel(parent,position,tag,title)
h = uipanel('Tag',tag,'Title',title, 'FontSize', 12, ...
            'Parent',parent,'Units','normalized','Position',position);
        
% -------------------------------------------------------------------------
% Creates and styles a uitable
function h = create_table(parent,position,tag)
h = uitable('Tag',tag,'FontSize', 12, ...
            'Parent',parent,'Units','normalized','Position',position);
   
% -------------------------------------------------------------------------
% Creates and styles a uibutton
function h = create_button(parent,position,tag,str)
h = uicontrol('Style','pushbutton','Tag',tag,'String',str, 'FontSize', 12, ...
            'Parent',parent,'Units','normalized','Position',position);
    
% -------------------------------------------------------------------------
% Creates and styles an editable text box
function h = create_editbox(parent,position,tag,str)
h = uicontrol('Style','edit','Tag',tag,'HorizontalAlignment','left',...
            'String',str, 'FontSize', 12, ...
            'Parent',parent,'Units','normalized','Position',position,...
            'Max',2);
        
% -------------------------------------------------------------------------
% Creates and styles a text label
function h = create_text(parent,position,tag,str)
h = uicontrol('Style','text','Tag',tag,'String',str, 'FontSize', 12, ...
            'Parent',parent,'Units','normalized','Position',position);
        
% -------------------------------------------------------------------------
% Creates and styles a popup menu
function h = create_popupmenu(parent,position,tag,str)
h = uicontrol('Style','popupmenu','Tag',tag,'String',str, 'FontSize', 12, ...
            'Parent',parent,'Units','normalized','Position',position);

% -------------------------------------------------------------------------
% Displays the selected analysis in a separate figure window
function display_results(hObject, handles)

% Set constants
FLU = [1692,28330]; % death rate for seasonal flu (per season)

% Get results
DCM = handles.D.DCM;
GCM = handles.D.GCM;
PEB = handles.D.PEB;
BMA = handles.D.BMA;
BPA = handles.D.BPA;

% Get specific values for the selected country
cidx = get(handles.popupmenu_country,'Value');     % country index
M.T  = 180;                                        % six-month period
cY    = DCM{cidx}.Y;                               % empirical data
cEp   = DCM{cidx}.Ep;                              % posterior expectations
cCp   = DCM{cidx}.Cp;                              % posterior covariances

% Get data and priors
data = handles.D.data;
str  = handles.D.str;
pE   = handles.D.pE;
pC   = handles.D.pC;

% Plot types
PLOT_SECOND_LEVEL_FX = 1;
BPA_OVER_COUNTRIES   = 2;
COUNTRY_DIFFERENCES  = 3;
COUNTRY_PREDICTION   = 4;
COUNTRY_LATENT       = 5;
COUNTRY_SENSITIVITY  = 6;
SOCIAL_DISTANCING    = 7;
HERD_IMMUNITY        = 8;
REPRODUCTION_RATIO   = 9;

% Plot type descriptions
plot_str{PLOT_SECOND_LEVEL_FX} = 'Second level effects (PEB)';
plot_str{BPA_OVER_COUNTRIES}   = 'Bayesian parameter average (BPA) over countries';
plot_str{COUNTRY_DIFFERENCES}  = 'Differences between countries';
plot_str{COUNTRY_PREDICTION}   = 'Selected country: projections';
plot_str{COUNTRY_LATENT}       = 'Selected country: latent states';
plot_str{COUNTRY_SENSITIVITY}  = 'Selected country: sensitivity analysis (slow)';
plot_str{SOCIAL_DISTANCING}    = 'Selected country: social distancing';
plot_str{HERD_IMMUNITY}        = 'Selected country: herd immunity';
plot_str{REPRODUCTION_RATIO}    = 'Selected country: reproduction ratio (slow)';

% Populate drop-down menu
idx = get(handles.popupmenu_view,'Value');
set(handles.popupmenu_view,'String',plot_str);
set(handles.popupmenu_view,'Value',idx);

% Current plot type
plot_type = idx;

% Prepare results figure
spm_figure('GetWin','Results');
spm_clf;
handles.axes_plot = gca;

% Plot
switch plot_type
    case PLOT_SECOND_LEVEL_FX
        
        help_txt = [...
        '(between country effects). This figure shows the relationship between ' ...
        'certain parameters of the generative model and the explanatory variables ' ...
        'in a general linear model of between country effects. The examples are ' ...
        'based upon a ranking of the absolute value of the second level parameter; ' ...
        'namely, the contribution of an explanatory variable to a model parameter. ' ...
        'The lower panel shows the (absolute) parameters in image format']; 
        
        % assemble parameters
        P     = [];             % posterior expectations
        for i = 1:numel(DCM)
            P(:,i) = spm_vec(DCM{i}.Ep);
        end
        P     = P(PEB.Pind,:);
        Pname = str.names(PEB.Pind);

        % find largest absolute (second level) effects and plot
        Ep      = abs(PEB.Ep);
        Ep(:,1) = 0;
        Sp      = sort(Ep(:),'descend');

        for i = 1:4
            [I,J] = find(Ep == Sp(i));
            subplot(2,4,i), plot(PEB.M.X(:,J),P(I,:),'.','MarkerSize',32,'Color',[0.8 0.8 1])
            xlabel(PEB.Xnames{J}),ylabel([ 'log ' Pname{I}])
            title(Pname{I},'FontSize',16),axis square, box off
        end

        % GLM (second level) parameters
        subplot(2,4,5:8)
        i = PEB.Pind;
        imagesc(Ep), title('Parameters of GLM','FontSize',16)
        set(gca,'XTick',1:numel(PEB.Xnames) ,'Xticklabel',PEB.Xnames)
        set(gca,'YTick',1:numel(str.names(i)),'Yticklabel',str.names(i))
        set(gca,'XTickLabelRotation',90)
        axis square, box off   
        
    case BPA_OVER_COUNTRIES
        help_txt = [...
        '(Bayesian parameter averages). This figure reports the Bayesian parameter ' ...
        'averages over countries following a hierarchical or parametric empirical ' ...
        'Bayesian analysis that tests for - and applies shrinkage priors to - ' ...
        'posterior parameter estimates for each country. The upper panel shows the ' ...
        'parameters as estimated in log space, while the lower panel shows the ' ...
        'same results for the corresponding scale parameters (scale parameters are ' ...
        'nonnegative parameters). The blue bars report posterior expectations, ' ...
        'while the thin red bars are prior expectations. The pink bars denote 90% ' ...
        'Bayesian confidence or credible intervals. One can interpret these ' ...
        'parameters as the average value for any given parameter of the generative ' ...
        'model, to which a random (country specific) effect is added to generate ' ...
        'the ensemble dynamics for each country. In turn, these ensemble ' ...
        'distributions determine the likelihood of various outcome measures under ' ...
        'larger number (i.e., Poisson) assumptions.'];

        Ep = BPA.Ep;                                  % posterior expectations                         
        Cp = BPA.Cp;                                  % posterior covariances

        subplot(2,1,1)
        spm_plot_ci(Ep,Cp), hold on, bar(spm_vec(pE),1/4), hold off
        xlabel('Parameter (see priors table for key)','FontSize',12);
        ylabel('log parameters','FontSize',16)
        axis square, box off, grid on

        subplot(2,1,2)
        spm_plot_ci(Ep,Cp,[],[],'exp')
        set(gca,'yLim',[0 32])
        xlabel('Parameter (see priors table for key)','FontSize',12);
        ylabel('Parameters','FontSize',16)
        axis square, box off, grid on
        
    case COUNTRY_DIFFERENCES
        help_txt = [ ...
            '(differences among countries). This figure reports the differences among ' ...
            'countries in terms of selected parameters of the generative model, ' ...
            'ranging from the size of a cell (i.e., the effective size of an infected ' ...
            'population), through to the probability of dying when in critical care. '...
            'Interesting country specific differences here include an apparent '...
            'attenuation of social distancing responses, relative to other countries, '...
            'in the United States and Australia. The blue bars represent the posterior '...
            'expectations, while the pink bars are 90% Bayesian credible intervals. '...
            'Notice that these intervals are not symmetrical about the mean because '...
            'scale parameters are plotted here - as opposed to the log parameters. The '...
            'next figure illustrates the predictions - in terms of new deaths and '...
            'cases - based upon these parameter estimates.'];

        % assemble parameters
        P     = [];                                   % posterior expectations
        C     = [];                                   % posterior variances
        for i = 1:numel(DCM)
            P(:,i) = spm_vec(DCM{i}.Ep);
            C(:,i) = diag(DCM{i}.Cp);
        end

        % report selected parameters (see spm_COVID_priors)
        p     = [2,4,5,7,8,9,11,13,15,16,17,19];
        for i = 1:length(p)

            spm_figure('GetWin','Results');
            
            % posterior density
            subplot(4,3,i)
            Ep   = P(p(i),:);
            Cp   = C(p(i),:);
            spm_plot_ci(Ep',Cp,[],[],'exp'), hold on
            title(str.names{p(i)},'FontSize',16)
            xlabel('country'), axis square, box off

            % country with greatest map estimate (and United Kingdom)
            [d,j] = max(exp(Ep));
            text(j,d,data(j).country,'FontSize',8);
            [d,j] = min(exp(Ep));
            text(j,d,data(j).country,'FontSize',8);
            j     = find(ismember({data.country},'United Kingdom'));
            text(j,exp(Ep(j)),'*','FontSize',12,'Color','r','HorizontalAlignment','center');

        end
        
    case COUNTRY_PREDICTION
        
        help_txt = [ ...
        '(predicted outcomes). This figure provides an example of predicted new ' ...
        'deaths and cases (and recovery rates and its critical care unit ' ...
        'occupancy) for an exemplar country; here, the United Kingdom. The panels ' ...
        'on the left shows the predicted outcomes as a function of weeks. The blue ' ...
        'line corresponds to the expected trajectory, while the shaded areas are ' ...
        '90%Bayesian credible intervals. The black dots represent empirical data, ' ...
        'upon which the parameter estimates are based. The lower right panel shows ' ...
        'the parameter estimates for the country in question. As in previous ' ...
        'figures, the prior expectations are shown as in bars over the posterior ' ...
        'expectations (and credible intervals). The upper right panel illustrates ' ...
        'the equivalent expectations in terms of cumulative deaths. The key point ' ...
        'to take from this figure is the quantification of uncertainty inherent in ' ...
        'the credible or confidence intervals. In other words, uncertainty about ' ...
        'the parameters propagates through to uncertainty in predicted outcomes. ' ...
        'This uncertainty changes over time because of the non-linear relationship ' ...
        'between model parameters and ensemble dynamics. By model design, one can ' ...
        'be certain about the final states; however, uncertainty about cumulative ' ...
        'death rates itself accumulates. The mapping from parameters, through ' ...
        'ensemble dynamics to outcomes is mediated by latent or hidden states. The ' ...
        'trajectory of these states is illustrated in the subsequent figure. ' ...
        'Public Health England estimates that on average 17,000 people have died ' ...
        'from the flu in England annually between 2014/15 and 2018/19. However, ' ...
        'the yearly deaths vary widely, from a high of 28,330 in 2014/15 to a low ' ...
        'of 1,692 in 2018/19. Public Health England does not publish a mortality ' ...
        'rate for the flu.'];

        spm_COVID_ci(cEp,cCp,cY);
        
        % add seasonal flu rates
        subplot(2,2,2), hold on
        x   = get(gca,'XLim');
        plot(x,[FLU(1) FLU(1)],'-.r',x,[FLU(2) FLU(2)],'-.r')
        spm_axis tight

    case COUNTRY_LATENT
        
        help_txt = [ ...
        '(latent causes of observed consequences). The upper panels reproduce the ' ...
        'expected trajectories of the previous figure, for an example country ' ...
        '(here the United Kingdom) in. Here, the expected death rate is shown in ' ...
        'blue, new cases in red, predicted recovery rate in orange and CCU ' ...
        'occupancy in puce. The black dots correspond to empirical data. The lower ' ...
        'four panels show the evolution of latent (ensemble) dynamics, in terms of ' ...
        'the expected probability of being in various states. The first (location) ' ...
        'panel shows that after about 5 to 6 weeks, there is sufficient evidence ' ...
        'for the onset of an episode to induce social distancing, such that the ' ...
        'probability of being found at work falls, over a couple of weeks to ' ...
        'negligible levels. At this time, the number of infected people increases ' ...
        '(to about 30%) with a concomitant probability of being infectious a few ' ...
        'days later. During this time, the probability of becoming immune ' ...
        'increases monotonically and saturates, within this cell, at about 20 ' ...
        'weeks. Clinically, the probability of becoming symptomatic rises to about ' ...
        '20%, with a small probability of developing acute respiratory distress ' ...
        'and, possibly death. In terms of testing, there is a progressive increase ' ...
        'in the number of people tested, with a concomitant decrease in those ' ...
        'untested or waiting for their results. Interestingly, initially the ' ...
        'number of negative tests increases monotonically, while the proportion of ' ...
        'positive tests catches up during the peak of the episode. Under these ' ...
        'parameters, the entire episode lasts for about 12 weeks or three months. ' ...
        'The increase in (herd) immunity is interesting and will become important ' ...
        'later. One might ask to what extent these trajectories depend upon ' ...
        'different model parameters. This is quantified in the next figure.'];
        
        [Z,X] = spm_COVID_gen(cEp,M,1:3);
        spm_COVID_plot(Z,X,cY);
        
    case COUNTRY_SENSITIVITY
        help_txt = [
        '(sensitivity analysis). These panels show the change in outcome measures ' ...
        '(upper panel: death rate. lower panel: new cases). The bar charts are the ' ...
        'derivatives of outcomes with respect to each of the parameters. Positive ' ...
        'values (on the right) exacerbate new cases when increased, while, ' ...
        'conversely, negative values (on the left) decrease new cases. As one ' ...
        'might expect, increasing social distancing, bed availability and the ' ...
        'probability of survival outside critical care, tend to decrease death ' ...
        'rate. Interestingly, increasing both the period of symptoms and ARDS ' ...
        'decreases overall death rate, because there is more time to recover to an ' ...
        'asymptomatic state. The lower panel shows the second order derivatives or ' ...
        'sensitivity. The next figure focuses on the effects of social distancing ' ...
        'as a way of ameliorating the impact on deaths.'];

        
        % sensitivity analysis in terms of partial derivatives
        %--------------------------------------------------------------------------
        [ddY,dY] = spm_diff(@(P,M,U)spm_COVID_gen(P,M,U),cEp,M,1,[1,1]);
        Np       = spm_length(cEp);

        % cumulative effects over time
        %--------------------------------------------------------------------------
        DY    = sum(dY);
        for i = 1:Np
            DDY{i} = sum(ddY{i});
        end
        DDY   = spm_cat(DDY');

        % plot results
        %--------------------------------------------------------------------------
        spm_figure('GetWin','Results');
        subplot(2,1,1)
        bar(DY)
        set(gca,'XTick',1:Np,'Xticklabel',str.names,'FontSize',8)
        ylabel('First-order sensitivity','FontSize',16), box off
        camorbit(90,0),axis square

        subplot(2,1,2)
        imagesc(DDY)
        set(gca,'YTick',1:Np,'Yticklabel',str.names,'FontSize',8)
        title('Second-order sensitivity','FontSize',16), box off
        axis square
        
    case SOCIAL_DISTANCING
        help_txt = [
        '(the effects of social distancing). This figure uses the same format as '...
        'Figure 9. However, here trajectories are reproduced under different '...
        'levels of social distancing; from zero through to 4 (in 16 steps). This '...
        'parameter is the exponent applied to the probability of not being '...
        'infected. In other words, it scores the sensitivity of social distancing '...
        'to the prevalence of the virus in the population. In this example (based '...
        'upon posterior expectations for the United Kingdom), death rates (per '...
        'day) and underlying latent states of the population decrease '...
        'progressively with social distancing. The cumulative death rate is shown '...
        'as a function of social distancing in the upper right panel. The vertical '...
        'line corresponds to the posterior expectation of the social distancing '...
        'exponent for this country. In the next figure, we repeat this analysis '...
        'but looking at the effect of herd immunity.'];

        % increase social distancing exponent from 0 to 4
        %--------------------------------------------------------------------------
        P     = cEp;                                 % expansion point
        sde   = linspace(0,4,16);                   % range of social distancing
        P.Rin = BPA.Ep.Rin;                         % adjust contacts at home
        S     = sde;
        for i = 1:numel(sde)

            % social distancing exponent
            %----------------------------------------------------------------------
            P.sde = cEp.sde + log(sde(i) + eps);
            [Y,X] = spm_COVID_gen(P,M,1);
            S(i)  = sum(Y(:,1));

            % plot results and hold graph
            %----------------------------------------------------------------------
            spm_figure('GetWin','Results');
            spm_COVID_plot(Y,X)
            for j = 1:6, subplot(3,2,j), hold on, end
        end

        % cumulative deaths as a function of social distancing
        %--------------------------------------------------------------------------
        subplot(3,2,2), hold off
        plot(sde,S,[1 1]*exp(cEp.sde),[min(S) max(S)],'-.')
        title('Social distancing','FontSize',16),
        xlabel('social distancing exponent')
        ylabel('cumulative deaths')
        axis square,box off

        disp('lifes saved'), disp(max(S) - min(S))
        
    case HERD_IMMUNITY
        help_txt = [ ...
        '(herd immunity). This figure reproduces the format of the previous '...
        'figure. However, here, we increased the initial proportion of the cell '...
        '(i.e., population) who were initially immune. Increasing the initial '...
        'immunity dramatically decreases death rates with a fall in the cumulative '...
        'deaths from several thousand to negligible levels with a herd immunity of '...
        'about 70%. The dashed line in the upper panel shows the equivalent deaths '...
        'over the same time period due to seasonal flu (based upon 2014/19 '...
        'figures). This death rate would require an initial or herd immunity of '...
        'about 60%. It is interesting to return to Figure 6 and identify at what '...
        'point - during the course of the infection episode - this level of herd '...
        'immunity is obtained.'...
        '-------------------------------------------------------------------------- '...
        'The Department for Transport (DfT) has announced there were 1,784 '...
        'reported road deaths in 2018, compared to 1,793 reported in 2017 - a 1% '...
        'fall. There were 25,511 people seriously injured in reported road traffic '...
        'accidents in 2018, compared to 24,831 in 2017 - a 3%year-on-year '...
        'increase '...
        '-------------------------------------------------------------------------- '];

        % progressively increase initial immunity
        %--------------------------------------------------------------------------
        m     = linspace(0,1,16);
        P     = cEp;
        S     = m;
        for i = 1:numel(m)
            P.m   = log(m(i) + 1e-6);
            [Y,X] = spm_COVID_gen(P,M,1);
            S(i)  = sum(Y);

            % plot results and hold graph
            %----------------------------------------------------------------------
            spm_figure('GetWin','Results');
            spm_COVID_plot(Y,X)
            for j = 1:6, subplot(3,2,j), hold on, end

        end

        % plot
        %--------------------------------------------------------------------------
        subplot(3,2,2), hold off
        S(S < 0) = 0;
        plot(m,S,m,FLU(1)*m.^0,'-.r',m,FLU(2)*m.^0,'-.r')
        title('Herd immunity','FontSize',16), 
        xlabel('proportion immune')
        ylabel('cumulative deaths')
        axis square,box off

    case REPRODUCTION_RATIO
        help_txt = [ ...
        '(basic reproduction ratio). This figure plots the predicted death rates ' ...
        'for the country in question above the concomitant fluctuations ' ...
        'in the basic reproduction rate (R0) and herd immunity. The blue lines ' ...
        'represent the posterior expectations while the shaded areas corprespond to ' ...
        '90% credible intervals.'];

        U       = [1,4,5];
        spm_COVID_ci(DCM{cidx}.Ep,DCM{cidx}.Cp,[],U)

        % add seasonal flu rates
        subplot(2,2,2), hold on
        x       = get(gca,'XLim');
        plot(x,[FLU(1) FLU(1)],'-.r',x,[FLU(2) FLU(2)],'-.r')

end

% Finish by displaying help text
set(handles.edit_legend,'String',help_txt);

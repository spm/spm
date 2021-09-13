function [DCM] = DEM_COVID_LTLA(LA)
% FORMAT [DCM] = DEM_COVID_LTLA(LA)
% LA - local authority
%
% Demonstration of COVID-19 modelling
%__________________________________________________________________________
%
% This demonstration routine fixed multiple regional death by date and new
% cases data and compiles estimates of latent states for local
% authorities.
%
% Technical details about the dynamic causal model used here can be found
% at https://www.fil.ion.ucl.ac.uk/spm/covid-19/.
%__________________________________________________________________________
% Copyright (C) 2020 Wellcome Centre for Human Neuroimaging

% Karl Friston
% $Id: DEM_CLIMATE_India.m 8151 2021-09-13 09:12:37Z karl $


cd('C:\Users\karl\Dropbox\Climate')
clear all, close all, clc

% download from web options
%--------------------------------------------------------------------------
options = weboptions('ContentType','table');
options.Timeout = 20;

% load data
%==========================================================================
D        = readtable('data.csv');

% get new cases by (lower tier) local authority
%--------------------------------------------------------------------------
vnames   = D.Properties.VariableNames;
sj       = find(ismember(vnames,'State_name'));
state    = unique(D(:,sj));

vars        = vnames(4:end);
vars([2 4 8 20 23 25 26]) = [];

%% organise via State
%--------------------------------------------------------------------------
Y     = struct([]);
for i = 1:numel(state)
    
    % get state
    %----------------------------------------------------------------------
    s   = find(ismember(D(:,sj),state(i,1)));
    disp(state(i,1))
    
    
    % create data structure
    %----------------------------------------------------------------------    
    for v = 1:numel(vars)
        
        Y(v,i).state = table2array(state{i,1});
        Y(v,i).type  = vars{v};
        Y(v,i).unit  = 'Normalised';
        Y(v,i).U     = v;
        Y(v,i).Y     = table2array(D(s,find(ismember(vnames,vars{v}))));
        d     = s;
        for j = 1:numel(s)
            year     = table2array(D(s(j),find(ismember(vnames,'Year'))));
            month    = table2array(D{s(j),find(ismember(vnames,'season'))});
            year     = num2str(year);
            if any(ismember({'NaN','all'},month))
                d(j) = datenum(['06-' year],'mm-yyyy');
            end
            if any(ismember({'Kharif'},month))
                d(j) = datenum(['03-' year],'mm-yyyy');
            end
            if any(ismember({'Rabi'},month))
                d(j) = datenum(['09-' year],'mm-yyyy');
            end
        end
        Y(v,i).date  = d;        
    end
end

%% Select states and data to averge
%==========================================================================
d0    = datestr(min(spm_vec(Y.date)),'dd-mm-yyyy');
for i = 1:size(Y,2)
   nY(i) = numel(spm_COVID_Y(Y(:,i),d0,0));
end
Y(:,nY < max(nY)) = [];


%% create grand average
%==========================================================================
clear A
for i = 1:size(Y,1)
    
    dates      = unique(spm_vec({Y(i,:).date}));
    
    A(i).state = 'Average';
    A(i).type  = Y(i).type;
    A(i).U     = Y(i).U;
    A(i).unit  = Y(i).unit;
    A(i).date  = dates;
    A(i).Y     = NaN(size(dates));
    
    for d = 1:numel(dates)
        y = [];
        for s = 1:size(Y,2)
            t = find(ismember(Y(i,s).date,dates(d)),1);
            y = [y Y(i,s).Y(t)];
        end
        s = ~(isnan(y) | y == 0);
        if any(s)
            A(i).Y(d,1) = mean(y(s));
        end
    end
end
A     = spm_COVID_Y(A,min(dates),1);

%% graphics
%==========================================================================
spm_figure('GetWin','Time series (average) I'); clf;
for i = 1:12
    subplot(4,3,i)
    plot(A(i).date,A(i).Y)
    datetick('x','yyyy')
    title(A(i).type)
end
spm_figure('GetWin','Time series (average) II'); clf;
for i = 1:7
    subplot(4,3,i)
    plot(A(i + 12).date,A(i + 12).Y)
    datetick('x','yyyy')
    title(A(i + 12).type)
end

% invert average
%----------------------------------------------------------------------
Y   = A;


%% fit each regional dataset
%==========================================================================
for r = 1:1 % size(Y,2)
    
    fprintf('%d out of %d\n',r,size(Y,2));
    
    % remove NANs, smooth and sort by date
    %----------------------------------------------------------------------
    y     = spm_COVID_Y(Y(:,r),d0,1);
        
    % data structure with vectorised data and covariance components
    %----------------------------------------------------------------------
    xY.y  = spm_vec(y.Y);
    xY.Q  = spm_Ce([y.n]);
    for i = 1:numel(y)
        hE(i) = 4 - log(var(y(i).Y));
    end
    
    %% get and set priors
    %----------------------------------------------------------------------
    spm_figure('GetWin','States'); clf;

    [pE,pC] = spm_CLIMATE_priors;

    M.date = d0;
    M.T    = '01-01-2020';
    u      = [2 3];
    [Y,X]  = spm_CLIMATE_gen(pE,M,u);
    spm_CLIMATE_plot(Y,X,u);
    
    %% model specification
    %======================================================================
    M.date = d0;
    M.Nmax = 32;                   % maximum number of iterations
    M.G    = @spm_CLIMATE_gen;     % generative function
    M.pE   = pE;                   % prior expectations (parameters)
    M.pC   = pC;                   % prior covariances  (parameters)
    M.hE   = hE + 2;               % prior expectation  (log-precision)
    M.hC   = 1/512;                % prior covariances  (log-precision)
    M.T    = y;                    % original data structure and times
    U      = [y.U];                % outputs to model
    
    % model inversion with Variational Laplace (Gauss Newton)
    %----------------------------------------------------------------------
    [Ep,Cp,Eh,F] = spm_nlsi_GN(M,U,xY);
    
    % save prior and posterior estimates (and log evidence)
    %----------------------------------------------------------------------
    DCM(r).M  = M;
    DCM(r).Ep = Ep;
    DCM(r).Eh = Eh;
    DCM(r).Cp = Cp;
    DCM(r).Y  = y;
    DCM(r).xY = xY;
    DCM(r).F  = F;
    
    % now-casting for this region and date
    %======================================================================
    spm_figure('GetWin','States'); clf;
    %----------------------------------------------------------------------
    M.T     = '01-01-2030';
    [Y,X,T] = spm_CLIMATE_gen(pE,M,[1 2 3]);
    spm_CLIMATE_plot(Y,X,[1 2 3]);
    
    % plot incidence
    %----------------------------------------------------------------------
    spm_figure('GetWin','States'); clf;
    subplot(2,1,1)
    spm_CLIMATE_ci(Ep,Cp,A,1,M)

    
end








function [DCM] = DEM_COVID_S
% FORMAT [DCM] = DEM_COVID_S
%
% Demonstration of COVID-19 modelling with stratified populations
%__________________________________________________________________________
%
% This demonstration routine uses a stratified population by age to fit
% death by date according to age bins. In brief, this uses the same kind of
% DCM for each age group; and the accompanying population densities
% are coupled via contact matrices; in other words, the number of people
% from another group I expect to be in contact with perday. In addition,
% some of the clinical and epidemiological parameters are group specific
% using prespecified profiles encoded in R. the parameters of the contact
% matrices are optimised and a reasonably uninformative priors.
%
% Technical details about the dynamic causal model used here can be found
% at https://www.fil.ion.ucl.ac.uk/spm/covid-19/.
%
% The (annotated) open source code creating these graphics is
% DEM_COVID_DASH.m
%__________________________________________________________________________
% Copyright (C) 2020 Wellcome Centre for Human Neuroimaging

% Karl Friston
% $Id: DEM_COVID_UTLA.m 7891 2020-07-07 16:34:13Z karl $


% NHS postcode data
%==========================================================================
% xy   = readmatrix('NHS_Postcode.csv','range',[2 1  2643867 2]);
% pc   = readmatrix('NHS_Postcode.csv','range',[2 5  2643867 5], 'OutputType','char');
% cc   = readmatrix('NHS_Postcode.csv','range',[2 12 2643867 12],'OutputType','char');
% G.X  = xy(:,1);
% G.Y  = xy(:,2);
% G.PC = pc;
% G.LT = cc;

% NHS provider code to post code
%--------------------------------------------------------------------------
% E    = importdata('etr.xlsx');
% E    = E.textdata(:,[1,10]);

% save('COVID_CODES','G','E')
% load codes
%--------------------------------------------------------------------------
load COVID_CODES

% retrieve recent data from https://coronavirus.data.gov.uk
%--------------------------------------------------------------------------
url  = 'https://c19downloads.azureedge.net/downloads/csv/';
websave('coronavirus-cases_latest.csv',[url,'coronavirus-cases_latest.csv']);

% retrieve recent data from https://www.england.nhs.uk
%--------------------------------------------------------------------------
URL  = 'https://www.england.nhs.uk/statistics/wp-content/uploads/sites/2/2020/07/';
for i = 0:4
    try
        dstr = datestr(datenum(date) - i,'dd-mmm-yyyy');
        if strcmp(dstr(1),'0'),dstr = dstr(2:end); end
        url  = [URL 'COVID-19-total-announced-deaths-' dstr '-1.xlsx'];
        websave('COVID-19-total-announced-deaths.xlsx',url);
        break
    end
end

% load (ONS Pillar 1) testing and (PHE) death-by-date data
%--------------------------------------------------------------------------
C    = importdata('coronavirus-cases_latest.csv');
D    = importdata('COVID-19-total-announced-deaths.xlsx');

% get death by date from each NHS trust
%--------------------------------------------------------------------------
DN   = datenum(D.textdata.Tab4DeathsByTrust(15,6:end - 4),'dd-mmm-yy');
NHS  = D.textdata.Tab4DeathsByTrust(18:end,3);
DA   = D.data.Tab4DeathsByTrust(3:end,2:end - 4);

% get new cases by (lower tier) local authority
%--------------------------------------------------------------------------
AreaCode = C.textdata(2:end,2);
AreaType = C.textdata(2:end,3);
AreaDate = C.textdata(2:end,4);
AreaName = C.textdata(2:end,1);

j        = find(ismember(AreaType,'Lower tier local authority'));
AreaCode = AreaCode(j);
AreaDate = AreaDate(j);
AreaName = AreaName(j);
AreaCase = C.data(j,4);

% organise via NHS trust
%--------------------------------------------------------------------------
clear D C
k     = 1;
for i = 1:numel(NHS)
    
    % get NHS Trust code
    %----------------------------------------------------------------------
    j = find(ismember(E(:,1),NHS(i)));
    if numel(j)
        
        % get postcode of NHS trust
        %------------------------------------------------------------------
        PC = E(j,2);
        
        % get local authority code of postcode
        %------------------------------------------------------------------
        g  = find(ismember(G.PC,PC));
        LA = G.LT(g);
        
        % get local authority code of NHS trust
        %------------------------------------------------------------------
        j  = find(ismember(AreaCode,LA));
        l  = find(ismember(G.LT,LA));
        
        if numel(j)
            
            % get deaths
            %--------------------------------------------------------------
            D(k).deaths = DA(i,:);
            D(k).date   = DN;
            D(k).name   = AreaName(j(1));
            D(k).NHS    = NHS(i);
            D(k).code   = LA;
            D(k).PC     = PC(1);
            D(k).X      = G.X(l);
            D(k).Y      = G.Y(l);
            
            % get cumulative cases
            %--------------------------------------------------------------
            CN    = datenum(AreaDate(j),'yyyy-mm-dd');
            for n = 1:numel(DN)
                [d,m]         = min(abs(CN - DN(n)));
                D(k).cases(n) = AreaCase(j(m));
            end
            
            % next region
            %--------------------------------------------------------------
            k = k + 1
        end
    end
end

% finding unique local authorities
%--------------------------------------------------------------------------
Area  = unique([D.code]);
k     = 1;
for i = 1:numel(Area)
    
    % accumulate deaths within each local authority
    %----------------------------------------------------------------------
    j     = find(ismember([D.code],Area(i)));
    DD(k) = D(j(1));
    for n = 2:numel(j)
        DD(k).deaths = DD(k).deaths + D(j(n)).deaths;
    end
    k     = k + 1;
end
D    = DD;

% find local authorities not accounted for
%--------------------------------------------------------------------------
UniqueArea = unique(AreaCode);
j          = unique(find(~ismember(UniqueArea,Area)));
UniqueArea = UniqueArea(j);

Ncode = cell2mat([D.code]');
Ncode = str2num(Ncode(:,3:end));
for i = 1:numel(UniqueArea)
    
    % find closest area in D
    %----------------------------------------------------------------------
    ncode = cell2mat(UniqueArea(i));
    ncode = str2num(ncode(:,3:end));
    [d,k] = min(abs(Ncode - ncode));
    
    % get local authority code of NHS trust
    %------------------------------------------------------------------
    j  = find(ismember(AreaCode,UniqueArea(i)));
    l  = find(ismember(G.LT,UniqueArea(i)));
    
    % supplement area
    %----------------------------------------------------------------------
    D(k).code(end + 1) = UniqueArea(i);
    D(k).X             = [D(k).X; G.X(l)];
    D(k).Y             = [D(k).Y; G.Y(l)];
    
    % add cumulative cases
    %--------------------------------------------------------------
    CN    = datenum(AreaDate(j),'yyyy-mm-dd');
    for n = 1:numel(DN)
        [d,m]         = min(abs(CN - DN(n)));
        D(k).cases(n) = D(k).cases(n) + AreaCase(j(m));
    end
end


% regional images
%--------------------------------------------------------------------------
n    = 512;
i    = ceil(n*G.X/7e5);
j    = ceil(n*G.Y/7e5);
i    = i(isfinite(j));
j    = j(isfinite(j));
UK   = sparse(i,j,1,n,2*n);

% create image array
%--------------------------------------------------------------------------
clear G DD
for k = 1:numel(D)
    i        = ceil(n*D(k).X/7e5);
    j        = ceil(n*D(k).Y/7e5);
    i        = i(isfinite(j));
    j        = j(isfinite(j));
    D(k).X   = mean(i);
    D(k).Y   = mean(j);
    LA       = spm_conv(double(logical(sparse(i,j,1,n,n))),1,1);
    ENG(:,k) = spm_vec(LA);
end
ENG = bsxfun(@rdivide,ENG,sum(ENG,2) + eps);

% plot map
%--------------------------------------------------------------------------
GRAPHICS = 1;
if GRAPHICS
    spm_figure('GetWin','UK - regional data'); clf;
    subplot(2,1,1)
    imagesc(1 - log(UK + 1)'), axis xy image off
end

for k = 1:numel(D)
    
    % prepend 16 days to timeseries
    %----------------------------------------------------------------------
    D(k).cases  = [zeros(1,16), D(k).cases];
    D(k).deaths = [zeros(1,16), D(k).deaths];
    D(k).date   = [(DN(1) - flip(1:16)) DN'];
    
    % assemble and smooth data matrix
    %----------------------------------------------------------------------
    s        = 7;                        % seven day average
    Y        = [spm_conv(D(k).deaths',s,0), gradient(spm_conv(D(k).cases,s))'];
    Y(Y < 0) = 0;
    D(k).YY  = Y;
    
    if GRAPHICS
        
        % show data
        %------------------------------------------------------------------
        subplot(2,1,1), hold on
        plot(D(k).X,D(k).Y,'.')
        axis image, axis off, box off
        title([D(k).name,D(k).PC])
        
        subplot(2,2,3), hold on
        plot(D(k).date,D(k).YY(:,1))
        datetick, axis square, title('Deaths')
        
        subplot(2,2,4), hold on
        plot(D(k).date,D(k).YY(:,2))
        datetick, axis square, title('New cases')
        drawnow
    end
    
end


% fit each regional dataset
%==========================================================================
for r = 1:numel(D)
    
    % get (Gaussian) priors over model parameters
    %----------------------------------------------------------------------
    [pE,pC] = spm_COVID_priors;
    
    % priors for this analysis
    %----------------------------------------------------------------------
    pE.N   = log(1/4);                    % population of local authority
    pC.N   = 1/16;
    pE.n   = 0;                           % initial number of cases (n)
    
    % variational Laplace (estimating log evidence (F) and posteriors)
    %======================================================================
    Y      = D(r).YY;
    
    % complete model specification
    %----------------------------------------------------------------------
    M.date = datestr(DN(1),'dd-mm-yyyy'); % date of first time point
    M.G    = @spm_COVID_gen;              % generative function
    M.FS   = @(Y)sqrt(Y);                 % feature selection  (link function)
    M.pE   = pE;                          % prior expectations (parameters)
    M.pC   = pC;                          % prior covariances  (parameters)
    M.hE   = [2 0];                       % prior expectation  (log-precision)
    M.hC   = 1/512;                       % prior covariances  (log-precision)
    M.T    = size(Y,1);                   % number of samples
    U      = [1 2];                       % outputs to model
    
    % model inversion with Variational Laplace (Gauss Newton)
    %----------------------------------------------------------------------
    [Ep,Cp]   = spm_nlsi_GN(M,U,Y);
    
    % save prior and posterior estimates (and log evidence)
    %----------------------------------------------------------------------
    DCM(r).M  = M;
    DCM(r).Ep = Ep;
    DCM(r).Cp = Cp;
    DCM(r).Y  = Y;
    
    % now-casting for this region and date
    %======================================================================-
    spm_figure('GetWin',D(r).name{1}); clf;
    %----------------------------------------------------------------------
    [Y,X] = spm_COVID_gen(DCM(r).Ep,DCM(r).M,[1 2]);
    spm_COVID_plot(Y,X,DCM(r).Y);
    
    
    %----------------------------------------------------------------------
    % Y(:,4) - effective reproduction ratio (R)
    % Y(:,5) - population immunity (%)
    % Y(:,8) - prevalence of infection (%)
    % Y(:,9) - number of infected at home, untested and asymptomatic
    %----------------------------------------------------------------------
    Y       = spm_COVID_gen(DCM(r).Ep,DCM(r).M,[4 5 8 9]);
    
    DR(:,r) = Y(:,1);                           % Reproduction ratio
    DI(:,r) = Y(:,2);                           % Prevalence of immunity
    DP(:,r) = Y(:,3);                           % Prevalence of infection
    DC(:,r) = Y(:,4);                           % Infected, asymptomatic people
    DT(:,r) = 1e4*exp(Ep.N)*Y(:,3)/exp(Ep.Tin); % New infections today
    
    
    % supplement with table of posterior expectations
    %----------------------------------------------------------------------
    subplot(3,2,2), cla reset, axis([0 1 0 1])
    title(D(r).name,'Fontsize',16)
    
    str      = sprintf('Effective population %.2f million',exp(Ep.N));
    text(0,0.9,str,'FontSize',12,'Color','b')
    
    str      = sprintf('Reproduction ratio %.2f',DR(end,r));
    text(0,0.8,str,'FontSize',12,'Color','b')
    
    str      = sprintf('Infected, asymptomatic people %.0f',DC(end,r));
    text(0,0.7,str,'FontSize',12,'Color','b')
    
    str      = sprintf('New infections today %.0f',DT(end,r));
    text(0,0.6,str,'FontSize',12,'Color','b')
    
    str      = sprintf('Prevalence of infection %.2f%s',DP(end,r),'%');
    text(0,0.5,str,'FontSize',12,'Color','b')
    
    str = sprintf('Prevalence of immunity %.1f%s',DI(end,r),'%');
    text(0,0.4,str,'FontSize',12,'Color','b')
    
    str = {'The prevalences refer to the effective population' ...
        'as estimated from the new cases and deaths (shown as' ...
        'dots on the upper left panel)'};
    text(0,0.0,str,'FontSize',9,'Color','k')
    
    spm_axis tight, axis off
    
end

save COVID_LA

load COVID_LA

str = {'Reproduction ratio',...
       'Prevalence of infection (%)'...
       'New infections per day'};
DD    = {DR,DP,DT};
for j = 1:numel(DD)
    
    % save
    %----------------------------------------------------------------------
    spm_figure('GetWin',str{j}); clf; colormap('pink')
    subplot(2,1,1)
    
    Y     = DD{j}';
    m     = max(Y(:));
    Y     = 64*Y/m;
    for i = 1:size(Y,2)
        
        % image
        %------------------------------------------------------------------
        E   = reshape(ENG*Y(:,i),n,n);
        image(64 - 2*E'), axis xy image off, drawnow
        dat = datestr(D(1).date(i),'dd-mmm-yy');
        text(128,256,{str{j},dat},...
            'HorizontalAlignment','Center',...
            'FontWeight','bold','FontSize',12);
        
        % save
        %------------------------------------------------------------------
        MOV(i) = getframe(gca);
    end
    
    % ButtonDownFcn
    %----------------------------------------------------------------------
    h = findobj(gca,'type','image');
    set(h(1),'Userdata',{MOV,16})
    set(h(1),'ButtonDownFcn','spm_DEM_ButtonDownFcn')
    
    subplot(2,1,2)
    imagesc(64 - E'), axis xy image off, drawnow
    dat = datestr(D(1).date(i),'dd-mmm-yy');
    title({str{j},dat},'FontSize',12)
    
    
    % List four most affected authorities
    %----------------------------------------------------------------------
    [y,k] = sort(DD{j}(end,:),'descend');
    for i = 1:4
        if y(i) > 32
            Atsr = sprintf('%s (%.0f)',D(k(i)).name{1},y(i));
        else
            Atsr = sprintf('%s (%.2f)',D(k(i)).name{1},y(i));
        end
        text(8 + D(k(i)).X,D(k(i)).Y,Atsr,...
            'HorizontalAlignment','Left',...
            'FontWeight','bold','FontSize',10);
    end
    
end


return

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
% Copyright (C) 2008-2022 Wellcome Centre for Human Neuroimaging

% get data
%==========================================================================
%
% regional populations (from Wikipedia)
%--------------------------------------------------------------------------
% South East 9,133,625
% London  8,908,081
% North West  5,479,093
% East of England  6,201,214
% West Midlands  5,900,757
% South West  5,599,735
% Yorkshire and the Humber  5,479,615
% East Midlands  4,804,149
% North East  2,657,909
%--------------------------------------------------------------------------
Pop(1) = 6.201;          %     {'East Of England'         }
Pop(2) = 8.908;          %     {'London'                  }
Pop(3) = 5.900 + 4.804;  %     {'Midlands'                }.
Pop(4) = 7.292 + 5.479;  %     {'North East And Yorkshire'}
Pop(5) = 5.479;          %     {'North West'              }
Pop(6) = 9.133;          %     {'South East'              }
Pop(7) = 5.599;          %     {'South West'              }

% retrieve recent data from https://coronavirus.data.gov.uk
%--------------------------------------------------------------------------
url  = 'https://c19downloads.azureedge.net/downloads/csv/';
websave('coronavirus-cases_latest.csv',[url,'coronavirus-cases_latest.csv']);

% retrieve recent data from https://www.england.nhs.uk
%--------------------------------------------------------------------------
url  = 'https://www.england.nhs.uk/statistics/wp-content/uploads/sites/2/2020/07/';
dstr = datestr(datenum(date) - 1,'dd-mmm-yyyy');
if strcmp(dstr(1),'0'),dstr = dstr(2:end); end
url  = [url 'COVID-19-total-announced-deaths-' dstr '-1.xlsx'];
websave('COVID-19-total-announced-deaths.xlsx',url);

% load data
%--------------------------------------------------------------------------
C    = importdata('coronavirus-cases_latest.csv');
D    = importdata('COVID-19-total-announced-deaths.xlsx');

DN   = D.textdata.Tab1DeathsByRegion(15,4:end - 4);
DA   = D.data.Tab3DeathsByAge(3:7,2:end - 4);

% convert data strings to numbers
%--------------------------------------------------------------------------
for i = 1:numel(DN)
    dn(i) = datenum(DN(i),'dd-mmm-yy');
end

% match death rates with cumulative cases
%--------------------------------------------------------------------------
j       = find(ismember(C.textdata(1,:),'Area type'));
i       = find(ismember(C.textdata(:,j),'Region'));
regions = unique(C.textdata(i,1));

% get cumulative cases for each region
%--------------------------------------------------------------------------
for r = 1:numel(regions)
    i     = find(ismember(C.textdata(:,1),regions{r}));
    cy{r} = C.data(i - 1,4);
    cn{r} = datenum(C.textdata(i,4),'yyyy-mm-dd');
end

% and associate with death rates
%--------------------------------------------------------------------------
for i = 1:numel(dn)
    for r = 1:numel(regions)
        d       = abs(cn{r} - dn(i));
        [d,j]   = min(d);
        CY(r,i) = cy{r}(j);
    end
end

% prepend 16 days to timeseries
%--------------------------------------------------------------------------
ns   = size(DA,1);
CR   = [zeros(1,16), sum(CY,1)];
DA   = [zeros(ns,16), DA];
dn   = [(dn(1) - flip(1:16)) dn];

% assemble and smooth data matrix
%--------------------------------------------------------------------------
s        = 7;                        % seven day average
Y        = [spm_conv(DA',s,0), gradient(spm_conv(CR,s))'];
Y(Y < 0) = 0;

% get (Gaussian) priors over model parameters
%==========================================================================
[pE,pC] = spm_COVID_priors;

% specific priors for this analysis
%--------------------------------------------------------------------------
pE.N   = log(sum(Pop));               % population size
pE.n   = 8;                           % initial cases
pC.N   = 0;

% group specific parameters
%--------------------------------------------------------------------------
r      = [1 1 1 1 1]/4;               % proportion resistant 
pE.r   = log(r + exp(-8));
pC.r   = logical(r)/256;

% contact matrices that couple ensemble densities
%--------------------------------------------------------------------------
Rin = [1 1 1 0 0;                     % effective number of contacts: home
       1 1 1 0 0;
       1 1 1 0 0;
       0 0 1 2 0;
       0 1 0 0 0];
Rou = [15 1  0 0 0;                   % effective number of contacts: work
       0 16 16 4 0;
       0 16 16 4 0;
       0 16 16 4 0;
       0 1  0  0 8];
pE.Rin = log(Rin + exp(-8));
pE.Rou = log(Rou + exp(-8));
pC.Rin = logical(Rin)/256;
pC.Rou = logical(Rou)/256;


% variational Laplace (estimating log evidence (F) and posteriors)
%==========================================================================

% complete model specification
%--------------------------------------------------------------------------
M.date = datestr(dn(1),'dd-mm-yyyy'); % date of first time point
M.G    = @spm_COVID_S;                % generative function
M.FS   = @(Y)sqrt(Y);                 % feature selection  (link function)
M.pE   = pE;                          % prior expectations (parameters)
M.pC   = pC;                          % prior covariances  (parameters)
M.hE   = [0 -2];                      % prior expectation  (log-precision)
M.hC   = 1/512;                       % prior covariances  (log-precision)
M.T    = size(Y,1);                   % number of samples
U      = 1:(ns + 1);                  % outputs to model

% model inversion with Variational Laplace (Gauss Newton)
%--------------------------------------------------------------------------
[Ep,Cp] = spm_nlsi_GN(M,U,Y);

% save prior and posterior estimates (and log evidence)
%--------------------------------------------------------------------------
DCM.M  = M;
DCM.Ep = Ep;
DCM.Cp = Cp;
DCM.Y  = Y;

save COVID_S

% show posterior estimates
%======================================================================
spm_figure('GetWin',['UK ' datestr(dn(end))]); clf;
%----------------------------------------------------------------------
M.T   = 385;
[Z,X] = spm_COVID_S(Ep,M,U);
spm_COVID_plot(Z,X,Y,[],U)

return

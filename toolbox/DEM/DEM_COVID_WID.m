function DCM = DEM_COVID_WID
% FORMAT DCM = DEM_COVID_WID
%
% Demonstration of COVID-19 modelling using variational Laplace
%__________________________________________________________________________
%
% This routine addresses the global factors that underwrite mortality and
% morbidity (and economic cost) by comparing different countries in terms
% of the epidemiological and sociobehavioural parameters that best explain
% the trajectory of daily cases and deaths. a subset of countries with a
% large population and a low vaccination rate are modelled by optimising a
% subset of country specific (free) parameters that capture differences in
% exposure to the virus and subsequent responses; in terms of quarantine,
% containment, resistance to infection and so on. The remaining model
% parameters are assumed to be conserved over countries and are based on
% posterior estimates derived from comprehensive timeseries data from the
% United Kingdom. The between country analyses are based upon available
% timeseries from Our World in Data
%
% The predictive validity of this modelling is established in terms of the
% accuracy of long-term forecasting up to 6 months - in countries that are
% modelled with sufficient accuracy (at least 50% variance in death rates 
% explained). A canonical correlation analysis is then used to establish
% the key parameters or factors that account for differences in fatalities
% over countries. Finally, the model is used to assess the impact of an
% equable vaccination programme starting over the next month using scenario
% modelling. This scenario modelling considers the impact of equable
% vaccination on cumulative deaths and gross domestic product. The
% conclusions are based upon a subset of countries accounting for over 50%
% of the world's population.
%
% Please see the annotations in this script for further details at each
% section of the analysis.
%__________________________________________________________________________
% Copyright (C) 2020 Wellcome Centre for Human Neuroimaging

% Karl Friston
% Copyright (C) 2008-2022 Wellcome Centre for Human Neuroimaging

% set up and preliminaries
%==========================================================================

% get data from all available countries and select countries with either a
% large population, or an impoverished vaccination rate
%--------------------------------------------------------------------------
cd('C:\Users\karl\Dropbox\Coronavirus')
close all, clear all, clc

D         = DATA_WID_data;
% [i,j]   = sort([D.population],'descend');
% country = {D.country};
% country = country(j(1:64));
UPDATE  = 1;
country = [ ...
    {'India'              }
    {'United States'      }
    {'Indonesia'          }
    {'Brazil'             }
    {'Bangladesh'         }
    {'Russia'             }
    {'Mexico'             }
    {'Japan'              }
    {'Vietnam'            }
    {'Turkey'             }
    {'Germany'            }
    {'Thailand'           }
    {'United Kingdom'     }
    {'France'             }
    {'Italy'              }
    {'South Africa'       }
    {'South Korea'        }
    {'Colombia'           }
    {'Spain'              }
    {'Argentina'          }
    {'Ukraine'            }
    {'Canada'             }
    {'Poland'             }
    {'Morocco'            }
    {'Peru'               }
    {'Malaysia'           }
    {'Cote d''Ivoire'     }
    {'Australia'          }
    {'Taiwan'             }
    {'Sri Lanka'          }
    {'Chile'              }
    {'Romania'            }
    {'Guatemala'          }
    {'Ecuador'            }
    {'Cambodia'           }
    {'Zimbabwe'           }
    {'Tunisia'            }
    {'Bolivia'            }
    {'Belgium'            }
    {'Dominican Republic' }
    {'Czechia'            }
    {'Greece'             }
    {'Azerbaijan'         }
    {'Portugal'           }
    {'Hungary'            }
    {'Israel'             }
    {'Switzerland'        }
    {'Hong Kong'          }
    {'Paraguay'           }
    {'Bulgaria'           }
    {'Lebanon'            }
    {'Denmark'            }
    {'Norway'             }
    {'Singapore'          }
    {'New Zealand'        }
    {'Ireland'            }
    {'Croatia'            }
    {'Uruguay'            }
    {'Mongolia'           }
    {'Lithuania'          }
    {'Slovenia'           }
    {'Latvia'             }
    {'Bahrain'            }
    {'Trinidad and Tobago'}];
    
for r = 1:numel(country)
    ri(r) = find(ismember({D.country},country{r}));
end


%% free parameters of model (i.e., country specific effects)
%==========================================================================
free  = {'n','o','m','exp','sde','qua',...
         'res','inn','lim','ons','rat',...
         'pcr','mob','rel','rol','fol'};
reset = {'pcr','mob'};
lock  = {'rol','fol'};

% (empirical) priors from the United Kingdom
%--------------------------------------------------------------------------
PCM   = load('DCM_UK.mat','DCM');
PCM   = PCM.DCM;
pE    = PCM.Ep;
pC    = spm_zeros(PCM.M.pC);
for i = 1:numel(free)
    pC.(free{i})  = PCM.M.pC.(free{i});
end
for i = 1:numel(reset)
    pE.(reset{i}) = PCM.M.pE.(reset{i});
end

% lock age-specific parameters
%==========================================================================
PC    = diag(spm_vec(pC));
for i = 1:numel(lock)
    [na,np] = size(pE.(lock{i}));
    if any(any(pC.(lock{i}))) && (na > 1)
        j = spm_fieldindices(pE,lock{i});
        j = reshape(j,na,np);
        for k = 1:np
            ii = j(:,k);
            PC(ii,ii) = PC(ii(1),ii(1));
        end
        fprintf(lock{i})
    end
end

% Cycle over selected countries to estimate model parameters
%==========================================================================
UK    = find(ismember({D.country},'United Kingdom'));
for k = 1:numel(ri)
    
    fprintf('%d out of %d\n',k,numel(ri));
    r         = ri(k);
    disp(D(r).country)
    try, clear Y; end
    
    % create data structure
    %----------------------------------------------------------------------
    Y(1).type = 'New cases';
    Y(1).unit = 'number/day';
    Y(1).U    = 2;
    Y(1).date = datenum(D(r).date);
    Y(1).Y    = D(r).cases;
    Y(1).h    = 0;
    
    Y(2).type = 'Daily deaths';
    Y(2).unit = 'number/day';
    Y(2).U    = 1;
    Y(2).date = datenum(D(r).date);
    Y(2).Y    = D(r).death;
    Y(2).h    = 2;
    
    Y(3).type = 'Positivity';
    Y(3).unit = 'percent';
    Y(3).U    = 23;
    Y(3).date = datenum(D(r).date);
    Y(3).Y    = D(r).positive * 100;
    Y(3).h    = 0;
    
    Y(4).type = 'New tests';
    Y(4).unit = 'number/day';
    Y(4).U    = 6;
    Y(4).date = datenum(D(r).date);
    Y(4).Y    = D(r).tests;
    Y(4).h    = 0;
    
    Y(5).type = 'Total vaccinated';
    Y(5).unit = 'number';
    Y(5).U    = 36;
    Y(5).date = datenum(D(r).date);
    Y(5).Y    = D(r).vaccine/1e6;
    Y(5).h    = 0;
    
    Y(6).type = 'R-ratio';
    Y(6).unit = 'ratio';
    Y(6).U    = 4;
    Y(6).date = datenum(D(r).date);
    Y(6).Y    = D(r).ratio;
    Y(6).h    = 0;
    
    % remove NANs, smooth and sort by date
    %----------------------------------------------------------------------
    [Y,S] = spm_COVID_Y(Y,'01-02-2020',16);
    
    % data structure with vectorised data and covariance components
    %----------------------------------------------------------------------
    xY.y  = spm_vec(Y.Y);
    xY.Q  = spm_Ce([Y.n]);
    hE    = spm_vec(Y.h);
    
    % country specific priors
    %======================================================================
    
    % population
    %----------------------------------------------------------------------
    N(4)    = D(r).aged_70_older/100;
    s       = (1 - N(4) - 1/2)/(70 - D(r).median_age);
    N(3)    = (70 - 35)*s;
    N(2)    = (1 - N(4) - N(3))*(35 - 15)/35;
    N(1)    = (1 - N(4) - N(3))*(15 - 0)/35;
    pE.N    = log(N(:)*D(r).population/1e6);
    
    % initial seeding
    %----------------------------------------------------------------------
    onset   = find(S(1:256,1) > max(S(1:256,1))/8,1) - 52;
    pE.n    = PCM.Ep.n - onset/8;
    pE.ons  = log(exp(PCM.Ep.ons) + onset);
    
    
    % testing rates
    %----------------------------------------------------------------------
    rel     = max(D(r).cases/D(r).population)/max(D(UK).cases/D(UK).population);
    pE.lim  = PCM.Ep.lim + log(rel);
    
    % vaccination rates
    %----------------------------------------------------------------------
    rel         = D(r).vrate;
    pE.rol(:,1) = PCM.Ep.rol(:,1) + log(rel);

    % hospital capacity
    %----------------------------------------------------------------------
    rel         = D(r).hospital_beds_per_thousand/D(UK).hospital_beds_per_thousand;
    if rel < 1
        hos     = erf(exp(PCM.Ep.hos))*rel;
        pE.hos  = log(erfinv(hos));
    elseif isnan(rel)
        pE.hos  = PCM.Ep.hos + log(1/512);
    end
    
    % initialise based on previous inversions
    %======================================================================
    if UPDATE
        try
            load DCM_OWID DCM
        catch
            UPDATE = 0;
        end
    end
    clear M
    M.P   = pE;
    if UPDATE
        
        % use previous inversion
        %------------------------------------------------------------------
        field = fieldnames(pE);
        for i = 1:numel(field)
            if any(any(pC.(field{i})))
                try
                    if all(size(pE.(field{i})) == size(DCM(k).Ep.(field{i})))
                        M.P.(field{i}) = DCM(k).Ep.(field{i});
                    end
                    fprintf('+')
                end
            end
        end
    end  
    
    
    % model specification
    %======================================================================
    M.Nmax = 16;                   % maximum number of iterations
    M.date = '01-02-2020';         % start date
    M.G    = @spm_SARS_gen;        % generative function
    M.FS   = @(Y)real(sqrt(Y));    % feature selection  (link function)
    M.pE   = pE;                   % prior expectations (parameters)
    M.pC   = PC;                   % prior covariances  (parameters)
    M.hE   = hE + 2;               % prior expectation  (log-precision)
    M.hC   = 1/512;                % prior covariances  (log-precision)
    M.T    = Y;                    % data structure
    
    U      = [Y.U];                % outputs to model
    
    % model inversion with Variational Laplace (Gauss Newton)
    %----------------------------------------------------------------------
    [Ep,Cp,Eh,F] = spm_nlsi_GN(M,U,xY);
    
    % save in DCM structure
    %----------------------------------------------------------------------
    DCM(k).country = D(r).country;
    
    DCM(k).M  = M;
    DCM(k).Ep = Ep;
    DCM(k).Eh = Eh;
    DCM(k).Cp = Cp;
    DCM(k).F  = F;
    DCM(k).S  = S;
    DCM(k).U  = U;
    DCM(k).D  = D(r);
    
    save('DCM_OWID.mat','DCM')
    
    % posterior predictions
    %======================================================================
    spm_figure('GetWin',DCM(k).country); clf;
    M.T     = datenum(date) + 512 - datenum('01-02-2020','dd-mm-yyyy');
    v       = [1 2 23];
    u       = [find(U == 1) find(U == 2) find(U == 23)];   
    [Z,X]   = spm_SARS_gen(DCM(k).Ep,M,v);
    spm_SARS_plot(Z,X,DCM(k).S(:,u),v)
    
    % data fits
    %----------------------------------------------------------------------
    spm_figure('GetWin',[DCM(k).country ' data fits']); clf;
    for i = 1:numel(U)
        subplot(4,2,i)
        spm_SARS_ci(Ep,Cp,S(:,i),U(i),M);
        title(Y(i).type,'FontSize',14), ylabel(Y(i).unit) 
    end
    
end

return



%% comparative analyses
%==========================================================================
% The following sections implement a sequence of analyses to identify
% models that are fit for purpose in explaining different countries, two
% establish the predictive validity of the ensuing models and then use them
% for scenario modelling, with a special focus on the impact of equable
% vaccination

% load DCM and identify those countries that have been modelled
%--------------------------------------------------------------------------
clear, load DCM_OWID
for r = 1:numel(DCM)
    try
        i(r) = numel(DCM(r).Ep.s);
    end
end
r     = find(i == 5);
DCM   = DCM(r);


%% check accuracy
%==========================================================================
% Run through the countries and assess the variance explained. If the model
% can predict at least 50% of the variability in daily deaths, then retain
% the country.
%--------------------------------------------------------------------------
spm_figure('GetWin','Accuracy of fits'); clf;
NT    = size(DCM(1).S,1);
for r = 1:numel(DCM)
    
    % unpack model and posterior expectations
    %----------------------------------------------------------------------
    M   = DCM(r).M;                                 % model (priors)
    Ep  = DCM(r).Ep;                                % posterior expectation
    Cp  = DCM(r).Cp;                                % posterior covariances
    S   = DCM(r).S;                                 % smooth timeseries
    U   = DCM(r).U;                                 % indices of outputs
    
    M.T  = NT;
    Z    = real(spm_SARS_gen(DCM(r).Ep,M,[1,15]));
    i    = isfinite(S(:,2));
    c    = corrcoef(S(i,2),Z(i,1));
    A(r) = c(1,2)^2;
    
    %  save fatalities
    %----------------------------------------------------------------------
    D(:,:,r) = Z;
    disp(r)
    
end

% illustrate the accuracy overall countries
%--------------------------------------------------------------------------
subplot(2,2,1), bar(A)
subplot(2,2,2), histogram(A,8)


%% retain accurate fits an estimate the total population accounted for by
% these models
%--------------------------------------------------------------------------
r     = find(A > .8);
DCM   = DCM(r);
A     = A(r);
D     = D(:,:,r);
for r = 1:numel(DCM) 
    NT(r)     = DCM(r).D.population;
    county{r} = DCM(r).D.country;
end

% print the total population accounted for by the subset of countries
%--------------------------------------------------------------------------
fprintf('total population: %.2f billion (%.0f%s of total population)\n\n',sum(NT)*1e-9,sum(NT)/7.9e7,'%')

% assign each country a colour
%--------------------------------------------------------------------------
for r = 1:numel(DCM)
    DCM(r).color = spm_softmax(randn(3,1))';
end

% illustrate data fits (best model)
%--------------------------------------------------------------------------
[d,r] = max(A)
Y     = DCM(r).M.T;
spm_figure('GetWin',[DCM(r).D.country ' best data fit']); clf;
for i = 1:numel(DCM(r).U)
    subplot(4,2,i)
    spm_SARS_ci(DCM(r).Ep,DCM(r).Cp,DCM(r).S(:,i),DCM(r).U(i),M);
    title(DCM(r).M.T(i).type,'FontSize',14), ylabel(Y(i).unit)
end

spm_figure('GetWin',[DCM(r).D.country ' Latent states']); clf;
[Z,X]   = spm_SARS_gen(DCM(r).Ep,M,[1 15]);
spm_SARS_plot(Z,X,DCM(r).S(:,2),[1 15])

% illustrate data fits (worst model)
%--------------------------------------------------------------------------
[d,r] = min(A)
r     = find(ismember(county,'India')); %%%
Y     = DCM(r).M.T;
spm_figure('GetWin',[DCM(r).D.country ' worst data fit']); clf;
for i = 1:numel(DCM(r).U)
    subplot(4,2,i)
    spm_SARS_ci(DCM(r).Ep,DCM(r).Cp,DCM(r).S(:,i),DCM(r).U(i),M);
    title(DCM(r).M.T(i).type,'FontSize',14), ylabel(Y(i).unit)
end

spm_figure('GetWin',[DCM(r).D.country ' Latent states']); clf;
[Z,X]   = spm_SARS_gen(DCM(r).Ep,M,[1 15]);
spm_SARS_plot(Z,X,DCM(r).S(:,2),[1 15])



%%  mortality burden
%==========================================================================
spm_figure('GetWin','Global burden'); clf;

DD    = 0;
DC    = 0;
for i = 1:numel(DCM)
    
    DC = DC + squeeze(D(:,2,i));
    subplot(2,2,2)
    plot(DC), hold on
    xlabel('days'),ylabel('cases'),title('Estimated deaths/day')
    axis square, box off, a = axis;
    
    DD = DD + squeeze(D(:,1,i));
    subplot(2,2,1)
    plot(DD), hold on
    xlabel('days'),ylabel('cases'),title('Reported deaths/day')
    axis square, box off, axis(a)
    
    DC = DC + squeeze(D(:,2,i));
    subplot(2,2,4)
    plot(cumsum(DC)), hold on
    xlabel('days'),ylabel('cases'),title('Cumulative estimate')
    axis square, box off, a = axis;
    
    DD = DD + squeeze(D(:,1,i));
    subplot(2,2,3)
    plot(cumsum(DD)), hold on
    xlabel('days'),ylabel('cases'),title('Cumulative reports')
    axis square, box off, axis(a)
    
end

% print estimated reported deaths
% https://www.nature.com/articles/d41586-022-00104-8
%--------------------------------------------------------------------------
fprintf('cumulative deaths estimated: %.2f million \n\n',sum(DC)*1e-6)
fprintf('cumulative deaths  reported: %.2f million \n\n',sum(DD)*1e-6)
fprintf('Ratio: %.2f \n\n',sum(DC)/sum(DD))

%% predictive validity
%==========================================================================
% This section quantifies the predictive validity of the models in relation
% to cumulative deaths over several months into the future. This is
% evaluated by withholding the past few months of data, estimating model
% parameters, and comparing the ensuing predictions with what actually
% happened. By repeating this for several countries one can get an idea of
% the predictive validity as a function of the period of forecasting.
%--------------------------------------------------------------------------

% find countries with the largest fatality rates
%--------------------------------------------------------------------------
clear deaths
for r = 1:numel(DCM)
    d         = DCM(r).D.death;
    deaths(r) = sum(d(isfinite(d)));
end
[deaths,jr]   = sort(deaths,'descend');

% jr  = find(ismember(county,'United Kingdom'));
clear Q
dd    = (0:3)*60 + du;
for j = 1:16
    
    % number of days and forecasting period (six months)
    %----------------------------------------------------------------------
    r     = jr(j);
    Nd    = size(DCM(r).S,1);
    du    = Nd - 6*30;
    d0    = datenum(DCM(r).M.date,'dd-mm-yyyy');
   
    % withhold most recent data
    %----------------------------------------------------------------------
    M     = DCM(r).M;
    U     = DCM(r).U;
    Y     = M.T;
    for i = 1:numel(Y)
        d         = find(Y(i).date < (du + d0));
        Y(i).Y    = Y(i).Y(d);
        Y(i).date = Y(i).date(d);
    end
    
    % create data structure for inversion
    %----------------------------------------------------------------------
    Y      = spm_COVID_Y(Y,'01-02-2020',0);
    xY.y   = spm_vec(Y.Y);
    xY.Q   = spm_Ce([Y.n]);
    hE     = spm_vec(Y.h);
    
    % replace reduced data in the model
    %----------------------------------------------------------------------
    M      = rmfield(M,'P');
    M.Nmax = 8;
    M.hE   = hE + 2;
    M.T    = Y;
    
    M.P   = M.pE;
    field = fieldnames(M.pE);
    for i = 1:numel(field)
        if all(size(M.pE.(field{i})) == size(DCM(r).Ep.(field{i})))
            M.P.(field{i}) = DCM(r).Ep.(field{i});
        end
        fprintf('+')
    end
    
    % model inversion with Variational Laplace (Gauss Newton)
    %----------------------------------------------------------------------
    [Ep,Cp,Eh,F] = spm_nlsi_GN(M,U,xY);
    
    % assess predictive accuracy
    %======================================================================
    spm_figure('GetWin',[DCM(r).D.country ' predictive accuracy']); clf;

    % data fits and predictions
    %----------------------------------------------------------------------
    M.T   = Nd;
    d     = 1:du;
    for i = 1:2
        subplot(3,1,i)
        [S,CS,s,C,t] = spm_SARS_ci(Ep,Cp,DCM(r).S(:,i),U(i),M); hold on
        plot(t(d),DCM(r).S(d,i),'or')
    end

    
    % cummulative deaths
    %----------------------------------------------------------------------
    subplot(3,1,3), cla
    
    S        = DCM(r).S(:,2); S(isnan(S)) = 0;
    CS       = cumsum(S);
    cs       = cumsum(s);
    Q(:,:,j) = [CS,cs];
    plot(t,CS,'k',t,cs)
    for i = 1:numel(dd)
        hold on
        plot(t([dd(i) dd(i)]),[0 max(cs)],':r')
        hold off
    end
    legend({'outcome','prediction'}), legend boxoff
    xlabel('date'), datetick('x','dd-mmm','keeplimits','keepticks')
    title('Forecast of cumulative deaths'), spm_axis tight
    
    % save temp copy
    %----------------------------------------------------------------------
    save OWID_Q Q
 
end

% quantify variance explained over countries
%--------------------------------------------------------------------------
dd(1) = [];

spm_figure('GetWin','Predictive accuracy'); clf;

for i = 1:numel(dd)
    subplot(3,3,i)

    plot(squeeze(Q(dd(i),1,:)),squeeze(Q(dd(i),2,:)),'.','Markersize',16)
    axis image, box off, title(sprintf('%i days',dd(i) - du))
    xlabel('recorded deaths'),ylabel('predicted deaths')
end

% correlation 
%--------------------------------------------------------------------------
for i = 1:size(Q,1)
    [r,p] = corrcoef(squeeze(Q(i,1,:)),squeeze(Q(i,2,:)));
    rho(i) = r(1,2);
    pho(i) = p(1,2);   
end

subplot(3,1,2), cla, hold off
plot(rho), hold on, plot(pho,':')
plot(ones(size(pho))*0.05,':k'), plot([du du],[0 1],'-.k')
box off, 
title('Correlation (and p - value)'), xlabel('days'), ylabel('correlation')


% percentage error
%--------------------------------------------------------------------------
subplot(3,1,3), cla, hold off
Q(Q(:) < 1) = 1;
e    = squeeze(Q(:,1,:) - Q(:,2,:))./squeeze(Q(:,1,:));
e    = abs(e)*100;
plot(e,'.c','MarkerSize',1/32),hold on
spm_plot_ci(mean(e,2)',var(e')), hold off
box off
set(gca,'YLim',[0 32])
title('Percentage error'), xlabel('days'), ylabel('percent')


%% canonical correlation analysis
%==========================================================================
% This subsection evaluates the role of vaccination and other model
% parameters in shaping the fatality. To evaluate the effect of model
% parameters on the epidemiological trajectory, it uses canonical
% correlation analysis: where the response variable is the trajectory over
% time from each country and the expansionary variable is a vector of
% (time invariant) epidemiological parameters that vary over countries.to
% interpret the relative contribution of different parameters in a
% quantitative sense, they are Euclidean normalised.

%% assemble response variable (deaths per hundred thousand)
%--------------------------------------------------------------------------
clear Y M
M.T   = size(S,1); 
for r = 1:numel(DCM)
    try
       Y(:,r) = 1e5*DCM(r).S(:,2)/DCM(r).D.population;
    end
    Y(isnan(Y(:))) = 0;
    % Y(:,r) = 1e5*spm_SARS_gen(DCM(r).Ep,M,15)/DCM(r).D.population;
    disp(r)
end

% transpose: one row per country
%--------------------------------------------------------------------------
Y     = spm_conv(Y(1:end,:),16,0);
Y     = Y';
EY    = Y;
Y     = spm_detrend(Y);

%% assemble expansion variables
%==========================================================================
parameters = { ...
        '1: population', ... 
        '2: youth', ... 
        '3: initial exposure', ... 
        '4: spread: containment', ... 
        '5: spread: quarantine', ... 
        '6: lockdown: sensitivity', ... 
        '7: lockdown: duration', ... 
        '8: resistance', ... 
        '9: testing rates', ... 
        '10: hospitalisation', ... 
        '11: vaccination: rate', ... 
        '12: vaccination: delay'}; 

X     = [];                                   % design matrix
X0    = [];                                   % confounds
for r = 1:numel(DCM)
        
    % explanatory variables
    %----------------------------------------------------------------------
    Ep = DCM(r).Ep;
    u  = [ ...
        log(sum(exp(Ep.N))), ...               % 1: population
        log(exp(Ep.N(1))/sum(exp(Ep.N))), ...  % 2: youth
        Ep.o, ...                              % 3: initial exposure
        Ep.m, ...                              % 4: spread: containment
        Ep.exp, ...                            % 5: spread: quarantine
        mean(Ep.sde), ...                      % 6: lockdown: sensitivity
        mean(Ep.qua), ...                      % 7: lockdown: duration
        mean(Ep.res), ...                      % 8: resistance
        mean(Ep.lim), ...                      % 9: testing rates
        mean(Ep.hos), ...                      % 10: hospitalisation
        mean(Ep.rol(:,1)),...                  % 11: vaccination: rate
        mean(Ep.rol(:,2))];                    % 12: vaccination: delay
    
    % constant and confounding variables
    %----------------------------------------------------------------------
    u0 = 1;
    
    % acummulate variables for this country
    %----------------------------------------------------------------------
    X  = [X;  real(u)];
    X0 = [X0; real(u0)];
    
end

% canonical correlation analysis
%--------------------------------------------------------------------------
X      = spm_en(spm_detrend(X));
CVA    = spm_cva(Y,X,X0,[],2)

spm_figure('GetWin','CVA'); clf;
spm_cva_ui('results',CVA)

subplot(2,2,1)
xlabel('days'), ylabel('excess death rate')
box off

subplot(2,2,3)
xlabel('predicted'), ylabel('observed'), title('Canonical correlations')
box off

subplot(2,2,4)
plot(EY',':c'), hold on, plot(mean(EY),'k')
xlabel(' days'), ylabel('deaths/day/100,000'), title('Average fatality rates')
axis  square, box off



%% unpack canonical variance and vectors (parameters)
%--------------------------------------------------------------------------
spm_figure('GetWin','CVA Vectors and variates'); clf;
subplot(2,1,1), cla
for r = 1:numel(DCM)
    x = CVA.w(r,1);
    y = CVA.w(r,2);

    plot(x,y,'.','Color',DCM(r).color,'MarkerSize',16), hold on
    text(x + 1/4,y,DCM(r).D.country,'Color',DCM(r).color,'FontSize',10);
end
xlabel('first canonical variate')
ylabel('second canonical variate')
title('Canonical variates')
axis square, box off

subplot(2,1,2), bar(CVA.W(:,1:2)','grouped')
xlabel('canonical vector')
ylabel('weight')
title('Canonical vectors')
axis square, box off, legend(parameters), legend boxoff






%% scenario modelling of equable vaccination rates
%==========================================================================
% In the subsequent sections, we drill down on the effect of equable
% vaccination by simulating future scenarios under different levels of
% vaccination rates and quantifying the effects in terms of deaths, cases
% and economic burden. First, we quantify the effect of introducing equable
% vaccination rates by setting the follow-up rollout of vaccination to
% prior expectations commencing 32 days after the current date
%--------------------------------------------------------------------------
spm_figure('GetWin','vaccination rates'); clf;

d      = datenum(date) - datenum(DCM(1).M.date,'dd-mm-yyyy') % duration of epidemi
pE.fol = log([0.0001 (d + 512) 32;
              0.001  (d + 128) 32;
              0.01   (d + 64 ) 32;
              0.01   (d + 32 ) 32]);

Z     = 0;
E     = 0;
for r = 1:numel(DCM)
        
    % unpack model and posterior expectations
    %----------------------------------------------------------------------
    M   = DCM(r).M;                                 % model (priors)
    Ep  = DCM(r).Ep;                                % posterior expectation
    Cp  = DCM(r).Cp;                                % posterior covariances
    S   = DCM(r).S;                                 % smooth timeseries
    U   = DCM(r).U;                                 % indices of outputs
    
    M.T = datenum(date) + 512 - datenum('01-02-2020','dd-mm-yyyy');
    z   = real(spm_SARS_gen(Ep,M,22));
    Z   = Z + z*DCM(r).D.population/100;
 
    % with equable rollout
    %----------------------------------------------------------------------
    Ep.fol = pE.fol;
    
    e   = real(spm_SARS_gen(Ep,M,22));
    E   = E + e*DCM(r).D.population/100;
    
    % graphics
    %----------------------------------------------------------------------
    subplot(2,1,1), plot(1:M.T,Z,1:M.T,E)
    xlabel('days'), ylabel('cumulative doses')
    title('Cumulative first doses (global)','FontSize',14)
    legend({'current','equable'}), box off
    
    subplot(2,1,2)
    col = DCM(r).color ;
    plot(1:M.T,z,'Color',col), hold on
    title('Current rollout (country-specific)','FontSize',14)
    xlabel('days'), ylabel('cumulative doses')
    i  = find(ismember(U,22));
    
    plot(S(:,i),'.','Color',col)
    text(64,256 - r*4,DCM(r).D.country,'Color',col,'FontSize',8)
    
    plot(1:M.T,e,'-.','Color',col), hold on
    title('Equable rollout (country-specific)','FontSize',14)
    xlabel('days'), ylabel('cumulative doses')
    set(gca,'YLim',[0 100]),box off, drawnow

end

% display results
%--------------------------------------------------------------------------
% https://www.who.int/docs/default-source/coronaviruse/act-accelerator/covax/costs-of-covid-19-vaccine-delivery-in-92amc_08.02.21.pdf
%--------------------------------------------------------------------------
dif = (E(end) - Z(end));
fprintf('total doses: %.2f billion \n',Z(end)/1e9)
fprintf('extra doses: %.2f billion (%.1f%s)\n\n',dif/1e9,100*dif/Z(end),'%' )
fprintf('extra cost: $%.2f billion (@ $1.8 per dose)\n\n',dif*1.8/1e9 )



%% scenario modelling - fatalities
%==========================================================================
% now repeat the procedure but focusing on fatalities
%--------------------------------------------------------------------------
spm_figure('GetWin','new deaths'); clf;
Z     = 0;
E     = 0;
for r = 1:numel(DCM)
    
    % unpack model and posterior expectations
    %----------------------------------------------------------------------
    M   = DCM(r).M;                                 % model (priors)
    Ep  = DCM(r).Ep;                                % posterior expectation
    Cp  = DCM(r).Cp;                                % posterior covariances
    S   = DCM(r).S;                                 % smooth timeseries
    U   = DCM(r).U;                                 % indices of outputs
    
    M.T = datenum(date) + 512 - datenum('01-02-2020','dd-mm-yyyy');
        
    % certified deaths and 28 day death rates
    %----------------------------------------------------------------------
    subplot(3,1,1), set(gca,'ColorOrderIndex',1)
    
    % certified and 28 day deaths
    %----------------------------------------------------------------------
    [~,~,y] = spm_SARS_ci(Ep,Cp,[],15,M); hold on
    Z = Z + y;
    title('Daily deaths','FontSize',14)
    
    % with equable rollout
    %----------------------------------------------------------------------
    Ep.fol = pE.fol;
    
    % certified deaths
    %----------------------------------------------------------------------
    subplot(3,1,1), set(gca,'ColorOrderIndex',2)
    [~,~,e] = spm_SARS_ci(Ep,Cp,[],15,M);

    E  = E + e;
    
    % lives saved (per country)
    %----------------------------------------------------------------------
    subplot(3,1,2), hold on
    plot(1:M.T,y - e,'Color',DCM(r).color)
    [i,j] = max(y - e);
    text(j,i,DCM(r).D.country,'Color',DCM(r).color,'FontSize',10)
    xlabel('days'), ylabel('Daily deaths')
    title('Daily lives saved (per country)','FontSize',14)
    spm_axis tight, box off

    % lives saved (globally)
    %----------------------------------------------------------------------
    subplot(3,1,3)
    plot(1:M.T,cumsum(Z),1:M.T,cumsum(E))
    xlabel('days'), ylabel('cumulative deaths')
    title('Cumulative deaths (global)','FontSize',14)
    legend({'current','equable'})
    spm_axis tight, box off
    
end

% display results
%--------------------------------------------------------------------------
disp('total deaths (millions)'), disp(sum(Z)*1e-6);
disp('lives saved  (millions)'), disp((sum(Z) - sum(E))*1e-6);

% $-QALY = $80,000 (green Book), 1 death = 7.6 QALY
%--------------------------------------------------------------------------
dif = (sum(Z) - sum(E));
fprintf('Cost ($-QALY = $80,000) of deaths (1 death = 7.6 QALY): %.0f x $608,000 = $%.2fB\n\n',dif, dif*608000/1e9)




%% scenario modelling - Hospital admissions
%==========================================================================
spm_figure('GetWin','Hospital admissions'); clf;
Z     = 0;
E     = 0;
for r = 1:numel(DCM)
    
    % unpack model and posterior expectations
    %----------------------------------------------------------------------
    M   = DCM(r).M;                                 % model (priors)
    Ep  = DCM(r).Ep;                                % posterior expectation
    Cp  = DCM(r).Cp;                                % posterior covariances
    S   = DCM(r).S;                                 % smooth timeseries
    U   = DCM(r).U;                                 % indices of outputs
    
    M.T = datenum(date) + 512 - datenum('01-02-2020','dd-mm-yyyy');
        
    % new cases
    %----------------------------------------------------------------------
    subplot(3,1,1), set(gca,'ColorOrderIndex',1)
    
    % rep aorted cases
    %----------------------------------------------------------------------
    [~,~,y] = spm_SARS_ci(Ep,Cp,[],16,M); hold on
    Z = Z + y;
    
    % with equable rollout
    %----------------------------------------------------------------------
    Ep.fol = pE.fol;
    
    % new cases
    %----------------------------------------------------------------------
    subplot(3,1,1), set(gca,'ColorOrderIndex',2)
    [~,~,e] = spm_SARS_ci(Ep,Cp,[],16,M); hold on

    E  = E + e;
   
    % cases avoided (per country)
    %----------------------------------------------------------------------
    subplot(3,1,2), hold on
    plot(1:M.T,y - e,'Color',DCM(r).color)
    [i,j] = max(y - e);
    text(j,i,DCM(r).D.country,'Color',DCM(r).color,'FontSize',10)
    xlabel('days'), ylabel('cumulative cases')
    title('Admissions avoided (per country)','FontSize',14)
    spm_axis tight, box off
    
    % Cases avoided (globally)
    %----------------------------------------------------------------------
    subplot(3,1,3)
    plot(1:M.T,cumsum(Z),1:M.T,cumsum(E))
    xlabel('days'), ylabel('cumulative cases')
    title('Cumulative admissions (global)','FontSize',14)
    legend({'current','equable'})
    spm_axis tight, box off
    
end

% display results
%--------------------------------------------------------------------------
disp('total cases (M)'),   disp(sum(Z)*1e-6);
disp('cases avoided (M)'), disp((sum(Z) - sum(E))*1e-6);

% hospital admissions
%--------------------------------------------------------------------------
dif   = sum(Z) - sum(E);
fprintf('Cost of hospital admissions: %.0f x $60,000 = Â£%.2fB\n\n',dif, dif*60000/1e9)


%% scenario modelling - economy
%==========================================================================
spm_figure('GetWin','economy'); clf;
Z     = 0;
E     = 0;
for r = 1:numel(DCM)
    
    % unpack model and posterior expectations
    %----------------------------------------------------------------------
    M   = DCM(r).M;                                 % model (priors)
    Ep  = DCM(r).Ep;                                % posterior expectation
    Cp  = DCM(r).Cp;                                % posterior covariances
    S   = DCM(r).S;                                 % smooth timeseries
    U   = DCM(r).U;                                 % indices of outputs
    
    M.T = datenum(date) + 512 - datenum('01-02-2020','dd-mm-yyyy');
        
    % certified deaths and 28 day death rates
    %----------------------------------------------------------------------
    subplot(3,1,1), set(gca,'ColorOrderIndex',1)
    
    % Gross domestic product
    %----------------------------------------------------------------------
    [~,~,y] = spm_SARS_ci(Ep,Cp,[],32,M); hold on
    title('Economic activity','FontSize',14)
    
    % with equable rollout
    %----------------------------------------------------------------------
    Ep.fol = pE.fol;
    
    % mobility
    %----------------------------------------------------------------------
    subplot(3,1,1), set(gca,'ColorOrderIndex',2)
    [~,~,e] = spm_SARS_ci(Ep,Cp,[],32,M); hold on
    
    y = (100 - y)*DCM(r).D.gdp_per_capita*DCM(r).D.population/365/100;
    e = (100 - e)*DCM(r).D.gdp_per_capita*DCM(r).D.population/365/100;
    y = y*1e-12;
    e = e*1e-12;

    if all(isfinite(y)) && all(isfinite(y))
        Z = Z + y;
        E = E + e;
    end
    
    % loss to the economy (per country)
    %----------------------------------------------------------------------
    subplot(3,1,2), hold on
    plot(1:M.T,y,'Color',DCM(r).color)
    [i,j] = max(y);
    text(j,i,DCM(r).D.country,'Color',DCM(r).color,'FontSize',10)
    
    xlabel('days'), ylabel('Trillion USD p.a.')
    title('loss of GDP p.a. (per country)','FontSize',14)
    spm_axis tight, box off
    
    % loss to the economy (globally)
    %----------------------------------------------------------------------
    subplot(3,1,3)
    plot(1:M.T,cumsum(Z),1:M.T,cumsum(E))
    xlabel('days'), ylabel('Trillion USD p.a.')
    title('Cumulative loss (global)','FontSize',14)
    legend({'current','equable'})
    spm_axis tight, box off
    
end

% display results
%--------------------------------------------------------------------------
disp('total cost (trillion USD)'), disp(sum(Z));
disp('savings with vaccination (billion USD)'), disp(1e3*(sum(Z) - sum(E)));

return





%% NOTES for using mixtures of explanatory variables based on SVD

% free parameters
%------------------------------------------------------------------------
fluct = {'pcr','mob','pro'};
free  = {'n','o','m','exp','sde','qua','inn','lim','ons','rat','rel'};
other = {'N','hos','rol'};

% explanatory variables (parameters)
%==========================================================================
factors = {'demographics','epidemiology','behaviour','spread','clinical','testing','vaccination'};

% % demographics
% %------------------------------------------------------------------------
% P.N   = 64;                   % (01) population size (millions)

demographics = {'N'};

% % epidemiological exposure
% %------------------------------------------------------------------------
% P.n   = exp(-8);              % (02) proportion of initial cases (cases)
% P.o   = 0.1;                  % (04) initial exposed proportion
% P.inn = 1;                    % (39) seasonal phase

epidemiology = {'n','o','inn'};

% % viral transmission (sociobehavioural)
% %------------------------------------------------------------------------
% P.sde = 4;                    % (07) time constant of lockdown
% P.qua = 64;                   % (08) time constant of unlocking

% % viral transmission (sociobehavioural)
% %------------------------------------------------------------------------
behaviour = {'sde','qua'};

% P.exp = 0.02;                 % (09) viral spreading (rate)
% P.m   = 0.1;                  % (05) relative eflux

% % viral transmission (sociobehavioural)
% %------------------------------------------------------------------------
spread = {'exp','m'};

% % clinical resources
% %------------------------------------------------------------------------
% P.hos = 1;                    % (10) admission rate (hospital) [erf]

clinical = {'hos'};

% % testing parameters
% %------------------------------------------------------------------------
% P.lim = [1/2 1   2   2]/1000; % (35) testing: capacity (P1, P2, & LFD)
% P.rat = [8   24  8   4];      % (36) testing: dispersion
% P.ons = [100 200 300 400];    % (37) testing: onset
% P.rel = 0.9;                  % (52) PCR testing of fatalities

testing = {'lim','rat','ons','rel'};

% % vaccination parameters
% %------------------------------------------------------------------------
% P.rol = [exp(-16) 365 32];    % (41) vaccination rollout (1st)

vaccination = {'rol'};


%% assemble explanatory variables (principal components)
%--------------------------------------------------------------------------
spm_figure('GetWin','factors'); clf;

X     = [];
C     = [];
nf    = numel(factors);
for f = 1:nf
    fields = eval(factors{f});
    x     = [];
    for i = 1:numel(fields)
        y     = [];
        for r = 1:numel(DCM)
            y(:,r) = real( spm_vec( DCM(r).Ep.(fields{i})(:,1) ) );
        end
        x   = [x; y];
    end
    
    % summarise with the principal singular variate for this factor
    %----------------------------------------------------------------------
    x       = x';
    x       = spm_detrend(x);
    x       = spm_en(x);
    [u,s,v] = spm_svd(x,.8);
    
    ii{f}   = size(X,2) + 1:size(u,2);
    X       = [X u];
    
    % plot
    %----------------------------------------------------------------------
    subplot(nf,2,2*f - 1), bar(diag(s)), title(factors{f}), axis square, box off
    subplot(nf,2,2*f), bar(v), title('weights'), axis square, box off

end


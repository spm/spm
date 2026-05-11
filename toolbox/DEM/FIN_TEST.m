function FIN_TEST
% Iterated simulations
% FORMAT FIN_TEST
%--------------------------------------------------------------------------

%% Iterated simulations
%==========================================================================
spm_clear
n     = 10;
m     = 1;
tab   = zeros(5,4,n,m);
DEM   = cell(n,m);
F     = cell(n,m);
for i = 1:n
    for j = 1:m

        % durations
        %------------------------------------------------------------------
        SIM.N  = fix((i - 1)*365/7);   % end point (weeks) [1]
        SIM.D  = 256;                  % depth of training data (weeks) [256]
        SIM.nT = 52;                   % duration of simulation (in weeks)
        SIM.dT = 4;                    % time between rebalancing (in weeks)
        SIM.Sr = -spm_invNcdf(.01);    % prior Sharpe ratio


        % specify number of indicator states and assets
        %------------------------------------------------------------------
        SIM.n  = 3;                    % number of indicator states
        SIM.m  = 6;                    % number of assets
        SIM.d  = 2;                    % order of detrending
        
        [t,f,d]      = DEM_FIN(SIM,1);   % evaluate
        tab(:,:,i,j) = table2array(t); % performance table
        DEM{i,j}     = d;              % generative model
        F{i,j}       = f               % ELBOs

    end
end


% bar chart results 
%--------------------------------------------------------------------------
VariableNames{1} = 'Annual RoR (%)';
VariableNames{2} = 'Volatility (%)';
VariableNames{3} = 'Sharpe ratio';
VariableNames{4} = 'Drawdown (%)';
RowNames = {'hold','ex post EV','ex post KL','ex ante EV','ex ante KL'};

for j = 1:m

    str   = sprintf('Annual performance - %i',j)
    spm_figure('GetWin',str); clf
    L     = full(spm_cat(F(:,j)));

    Dates = zeros(1,n);
    for i = 1:n
        Dates(i) = DEM{i}.G.date(end);
    end

    subplot(4,1,1)
    p = squeeze(tab(:,1,:,j));
    bar(p'), xlabel('year'), ylabel('percent')
    title(VariableNames{1})
    set(gca,'XTickLabel',datestr(Dates,'yyyy'))
    legend(RowNames)

    subplot(4,1,2)
    p = squeeze(tab(:,3,:,j));
    bar(p'), xlabel('year'), ylabel('ratio')

    subplot(4,2,5)
    p = squeeze(tab(:,1,:,j));
    p = minus(p,p(1,:));
    plot(p'), xlabel('year'), ylabel('percent difference')

    subplot(4,2,8)
    bar(-L(:,2)), xlabel('year'), ylabel('nats')
    title('Evidence for slaving')

    subplot(4,2,7)
    bar(L(:,1) - min(L(:,1))), xlabel('year'),ylabel('nats')
    title('ELBO')

    subplot(4,2,6)
    bar(-1./L(:,4),'c'), xlabel('year'), ylabel('weeks')
    hold on, plot([0,n],[4,4],'--r')
    title('Principal time constant')

    % Tabular results
    %----------------------------------------------------------------------
    avg = max(tab(:,:,:,j),[],3);
    Avg = array2table(avg);
    Avg.Properties.Description = '10 year maximum';
    Avg.Properties.VariableNames = VariableNames;
    Avg.Properties.RowNames      = RowNames
    
    avg = min(tab(:,:,:,j),[],3);
    Avg = array2table(avg);
    Avg.Properties.Description = '10 year minimum';
    Avg.Properties.VariableNames = VariableNames;
    Avg.Properties.RowNames      = RowNames
    
    avg = mean(tab(:,:,:,j),3);
    Avg = array2table(avg);
    Avg.Properties.Description = '10 year average';
    Avg.Properties.VariableNames = VariableNames;
    Avg.Properties.RowNames      = RowNames

end

return


% log evidence
%==========================================================================
load DEMFIN_FG.mat
L     = zeros(n,4,m);
for j = 1:m
    L(:,:,j) = full(spm_cat(F(1:n,j)));
end

% bar chart results
%--------------------------------------------------------------------------
spm_figure('GetWin','Marginal likelihood'); clf
B    = squeeze(tab(1,1,1:n,1:m));
R    = squeeze(tab(5,1,1:n,1:m));
E    = squeeze(L(:,1,:));
P    = squeeze(L(:,2,:));
T    = squeeze(L(:,4,:));

i    = find(abs(P(:)) > std(P(:))*6);
P(i) = 0;

x = ((1:m) + 16)*8;
y = 1:n;

D     = gradient(E);
S     = zeros(n,m + n*2);
for i = 1:n
    S(i,(1:m) + i*2) = D(i,:);
end

d0 = datenum('08-Nov-25','dd-mmm-yy');
dx = (16 + (1:(m + n*2)))*8;
dd = d0 - 7*dx;
plot(dd,sum(S))
datetick('x','mmm-yy','keeplimits')

subplot(4,3,1)
imagesc(x,y,E), ylabel('year'), xlabel('depth (weeks)')
title('Log evidence'), axis square

subplot(4,3,2)
imagesc(x,y,D), ylabel('year'), xlabel('depth (weeks)')
title('differential ELBO'), axis square

subplot(4,3,3)
imagesc(S), ylabel('nats/week'), xlabel('depth (8 weeks)')
title('differential ELBO'), axis square

subplot(4,1,2)
imagesc(flip(S,2)), ylabel('nats/week'), xlabel('depth (8 weeks)')
title('differential ELBO')

subplot(4,1,3)
d0 = datenum('08-Nov-25','dd-mmm-yy');
dx = 7*(16 + (1:(m + n*2))*8);
dd = d0 - dx;
plot(dd,sum(S)), title('Sum ELBO')
datetick('x','mmm-yy')

subplot(4,2,7), hold off
bar(x,mean(R)), ylabel('RoR'), xlabel('Depth')
hold on, plot(x,mean(B),'--r'), hold off
title('Average return')

subplot(4,2,8)
bar(x,sum(P < 0)), ylabel('Number +ve'),   xlabel('Depth')
title('enslavement')


% bar chart results
%==========================================================================
spm_figure('GetWin','log evidence'); clf

% log evidence
%--------------------------------------------------------------------------
load DEMFIN_all.mat
L     = zeros(n,4,m);
for j = 1:m
    L(:,:,j) = full(spm_cat(F(1:n,j)));
end

B    = squeeze(tab(1,1,1:n,1:m));
R    = squeeze(tab(5,1,1:n,1:m));
E    = squeeze(L(:,1,:));
P    = squeeze(L(:,2,:));
T    = squeeze(L(:,4,:));

i    = find(abs(P(:)) > std(P(:))*6);
P(i) = 0;
T    = -1./T;

x = ((1:m) + 7)*16;
y = 1:n;

subplot(4,2,1)
imagesc(x,y,-P),   ylabel('year'), xlabel('depth (weeks)')
title('log evidence for enslavement')

subplot(4,2,2)
imagesc(x,y,P < 0), ylabel('year'), xlabel('depth (weeks)')
title('Bayesian model selection')

subplot(4,2,3)
imagesc(x,y,T), ylabel('weeks'), xlabel('depth (weeks)')
title('Principal time constants')

subplot(4,2,4)
bar(mean(T),'c'), hold on, plot(T','.','MarkerSize',16)
ylabel('weeks'), xlabel('depth (weeks)')
title('Principal time constants')

subplot(4,2,5)
imagesc(R - B), ylabel('year'),   xlabel('depth (weeks)')
title('Average return over baseline')

subplot(4,2,6)
imagesc(R > 0), ylabel('year'), xlabel('depth (weeks)')
title('Years with positive return')

subplot(4,2,7), bar(mean(R)), spm_axis tight
hold on, plot(mean(B),'--r'), hold off
ylabel('percent per annum'), xlabel('depth (weeks)')
title('Average return')

subplot(4,2,8), bar(sum(R > 0)), spm_axis tight
ylabel('incidence'), xlabel('depth (weeks)')
title('Number of years with positive return')








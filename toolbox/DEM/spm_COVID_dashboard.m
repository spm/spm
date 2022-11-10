function spm_COVID_dashboard(DCM)
% Dashboard for coronavirus simulations
% FORMAT spm_COVID_plot(Y,X,Z)
% DCM.Ep
% DCM.M
% DCM.data
%
% This auxiliary routine plots the predicted prevalence of infection, the
% production rate and social distancing as a predicted timeline with
% annotated dates and statistics.
%__________________________________________________________________________
% Copyright (C) 2020 Wellcome Centre for Human Neuroimaging

% Karl Friston
% Copyright (C) 2008-2022 Wellcome Centre for Human Neuroimaging

% Plot outcomes
%==========================================================================
% The JBC will be responsible for setting the new COVID-19 Alert level to
% communicate the current level of risk clearly to the public. The alert
% levels are:
%
% Level 1 COVID-19 is not known to be present in the UK
%
% Level 2 COVID-19 is present in the UK, but the number of cases and
% transmission is low
%
% Level 3 A COVID-19 epidemic is in general circulation
%
% Level 4 A COVID-19 epidemic is in general circulation; transmission is
% high or rising exponentially
%
% Level 5 As level 4 and there is a material risk of healthcare services
% being overwhelmed
%
% The Government will engage with the devolved administrations to explore
% how the centre can operate most effectively across the UK, as it is
% established. Over time the Government will consider whether the JBC
% should form part of an extended infrastructure to address biosecurity
% threats to the UK, and whether the COVID-19 alert level system should be
% expanded to other potential infectious diseases.
%
% https://www.gov.uk/government/publications/our-plan-to-rebuild-the-uk-governments-covid-19-recovery-strategy/our-plan-to-rebuild-the-uk-governments-covid-19-recovery-strategy%--------------------------------------------------------------------------
% get dates and date strings
%--------------------------------------------------------------------------

% figure
%--------------------------------------------------------------------------
spm_figure('GetWin','Dashboard'); clf;

% get dates and date strings
%--------------------------------------------------------------------------
try
    d0 = DCM.date;                       % initial date
catch
    d0 = '25-Jan-2020';                  % for UK data
end
L0      = datenum('16-Mar-2020');        % official social distancing
L1      = datenum('10-May-2020');        % official relaxation
L2      = datenum('23-Mar-2020');        % official lockdown
n       = datenum(d0):datenum('1-Aug-2020');
dstr    = datestr(n,'mmmdd');
N       = numel(n);
DCM.M.T = N;
[Y,X]   = spm_COVID_gen(DCM.Ep,DCM.M,1:4);

% hard (threshold) strategy
%--------------------------------------------------------------------------
P    = spm_vecfun(DCM.Ep,@exp);             % posterior scale parameters
Prev = X{2}(:,2) + X{2}(:,3);            % prevalence of infection
Psde = spm_sigma(Prev,P.sde);

% plot hidden states and outcomes for this region or country
%==========================================================================

% thresholds
%--------------------------------------------------------------------------
U   = [0 .1 .2 .6 4];

% social distancing
%--------------------------------------------------------------------------
clf, subplot(3,1,1)
p   = X{2}(:,2)*100;
R   = Y(:,4);
D   = Y(:,1);
plot(n,1 - Psde,'LineWidth',2), hold on
plot(n,D/max(D),'b-.')
datetick('x','mmmdd')
axis([n(1), n(end) 0 1])
xlabel('time'),ylabel('proportion') ,box off
title('Social distancing','FontSize',16)

% prevalence of infection and reproduction rate
%--------------------------------------------------------------------------
subplot(3,1,2)
plot(n,p), hold on
plot(n,R,'r')
plot(n,ones(1,N),'-.r')
datetick('x','mmmdd')
axis([n(1), n(end) 0 8])
xlabel('time'),ylabel('proportion (%)') ,box off
title('Prevalence of infection','FontSize',16)

% loop over levels
%==========================================================================
for i = 1:numel(U)
    
    ui   = i/(16 + numel(U));
    
    % intervals for this level
    %----------------------------------------------------------------------
    if i == 1
        j  = find(p <= U(i + 1));
        Ul = U(i);
        Uu = U(i + 1);
    elseif i == numel(U)
        j  = find(p >= U(i));
        Ul = U(i);
        Uu = 8;
    else
        j  = find(p >= U(i) & p <= U(i + 1));
        Ul = U(i);
        Uu = U(i + 1);
    end
    
    % Timeline
    %----------------------------------------------------------------------
    subplot(3,1,1)
    for k = 1:numel(j)
        try
            fill(n(j(k) + [0 1 1 0]),[0 0 1 1],'r', ...
                'FaceAlpha',ui,'Edgecolor','none')
            jk   = j(k);
        end
    end
    
    % label
    %----------------------------------------------------------------------
    str = sprintf('%s (p = %.1f%s, R = %.2f)',dstr(jk,:),p(jk),'%',R(jk));
    text(n(jk),i/(1 + numel(U)),str,...
        'FontWeight','bold','HorizontalAlignment','right',...
        'Color','b','FontSize',10)
    
    % prevalence of infection
    %----------------------------------------------------------------------
    subplot(3,1,2)
    fill([n(1) n(end) n(end) n(1)],[Ul Ul Uu Uu],'b', ...
        'FaceAlpha',ui,'Edgecolor','none')
    
end
subplot(3,1,2)

plot([L0 L0],[0 8],'b')
j   = find(n == L0);
text(L0,6,dstr(j,:),'FontWeight','bold','FontSize',10,'Color','b')
plot([L1 L1],[0 8],'b:')
j   = find(n == L1);
text(L1,6,dstr(j,:),'FontWeight','bold','FontSize',10,'Color','b')
plot([L2 L2],[0 8],'b--')
j   = find(n == L2);
text(L2,7,dstr(j,:),'FontWeight','bold','FontSize',10,'Color','b')

leg{1}       = 'prevalence of infection';
leg{end + 1} = 'reproduction ratio';
leg{end + 1} = 'threshold';
for i = 1:numel(U)
    leg{end + 1} = sprintf('Level %d (%.1f%s)',i,U(i),'%');
end
leg{end + 1} = 'social distancing';
leg{end + 1} = 'relaxation';
leg{end + 1} = 'lockdown';
legend(leg),legend('boxoff')



function p = spm_sigma(x,u,s)
% reverse sigmoid function
% FORMAT p = spm_sigma(p,u)
% x    - probability
% u    - threshold
% u    - sensitivity (default four)
%
% p    - probability (0 < p < 1)
%
% This function is reverse sigmoid function that scales the input argument
% by the bias and flips the (biased) input. This provides a monotonically
% decreasing sigmoid function of the input that hits 50% at the threshold
% (u). The scaling ensures the probability at x = 0 is about one, for a
% suitably large sensitivity parameter s.
%--------------------------------------------------------------------------

% default sensitivity
%--------------------------------------------------------------------------
if nargin < 3, s = 4; end

% sigmoid function
%--------------------------------------------------------------------------
p = spm_phi(s*(u - x)/u);

return




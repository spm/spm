function spm_dcm_sessions
% Apply contrast vector to multiple DCM models
% FORMAT spm_dcm_sessions
%
% Contrasts are specified interactively and applied to a
% number of DCM models. This routine can be used, for example,
% to do Bayesian fixed or random effects analysis on 
% contrasts of DCM parameters.
%
% This function returns p-values for one-sided t-tests. The Bayesian 
% probabilites are p(effect_size > threshold) where the 
% threshold is specified by the user. If you wish to test for 
% effects being smaller than a threshold you can use negative 
% values when you specify the contrasts. p-values for two-sided
% tests are twice as large.
%
% In Bayesian fixed effects analysis the mean estimates from
% each DCM are weighted by their relative precision. Bayesian
% random effects analysis is based on the between-model variance.
% If the threshold is 0, and p is the random effects p-value 
% from classical inference then the Bayesian RFX probability value
% is 1-p. As usual, only the random effects procedures allow 
% you to make an inference about the population from which the 
% data (eg. subjects) are drawn.
%__________________________________________________________________________
% Copyright (C) 2004-2019 Wellcome Trust Centre for Neuroimaging

% Will Penny
% $Id: spm_dcm_sessions.m 7577 2019-04-24 08:59:56Z guillaume $


Finter = spm_figure('GetWin','Interactive');
Fgraph = spm_figure('GetWin','Graphics');

num_models = spm_input('Apply contrast to how many DCM models ? ','+1','r',[],1);
P      = spm_select(num_models,'^DCM.*\.mat$','select DCM*.mat files');

% Get contrast
str    = 'contrast for';
D      = spm_input(str,1,'b',{'A','B','C'});
load(deblank(P(1,:)));
con    = spm_dcm_contrasts(deblank(P(1,:)),D);

% Get threshold
str   = 'Threshold';
DCM.T = spm_input(str,1,'e',0,[1 1]);
T     = DCM.T;

% Get mean and variance of effect in each model
con_error = 0;
for model = 1:num_models,
    load(deblank(P(model,:)));
    mean_con(model) = full(con')*full(spm_vec(DCM.Ep));
    var_con(model)  = full(con')*full(DCM.Cp)*full(con);   
    disp(sprintf('Model %d, Contrast Mean=%1.4f, Contrast Variance=%1.4f',model,mean_con(model),var_con(model)));
    if var_con(model) == 0
        disp(sprintf('Error in spm_dcm_sessions: this contrast not valid for model %d',model));
        con_error=1;
    end
end

% Invalid contrast specified
%--------------------------------------------------------------------------
if con_error
    return
end

%- Bayesian fixed effects analysis 
%--------------------------------------------------------------------------
precision_con = 1./var_con;
precision     = sum(precision_con);
v   = 1/precision;
c   = sum((precision_con.*mean_con)/precision);
x   = c + [-32:32]*sqrt(v)*6/32;
p   = 1/sqrt(2*pi*v)*exp(-[x - c].^2/(2*v));
PP  = 1 - spm_Ncdf(T,c,v);
disp(' ');
disp(sprintf('Bayesian Fixed Effect Mean = %1.4f', full(c)));
disp(sprintf('Bayesian Fixed Effect Variance = %1.4f', full(v)));

figure(Fgraph)
subplot(2,1,1)
plot(x,p,[1 1]*T,[0 max(p)],'-.');
title({'Bayesian Fixed Effects: Posterior density of contrast',...
        sprintf('P(contrast > %0.2f) = %.1f%s',T,PP*100,'%')},...
       'FontSize',12)
xlabel('contrast');
ylabel('probability density');

i    = find(x >= T);
hold on
fill([x(i) fliplr(x(i))],[i*0 fliplr(p(i))],[1 1 1]*.8)
axis square, grid on
hold off

% Random effects analysis 
% - assuming an inverse gamma prior on the variance
% the posterior is a t-distribution 
%--------------------------------------------------------------------------
v    = std(mean_con)^2;
c    = mean(mean_con);
t    = (sqrt(num_models)*(c-T))/sqrt(v);
x    = c + [-32:32]*sqrt(v)*6/32;
p    = spm_Tpdf((sqrt(num_models)*(x-c))/sqrt(v),num_models-1);
PP   = spm_Tcdf(t,num_models-1);
disp(' ');
disp(sprintf('Bayesian Random Effect Mean = %1.4f', full(c)));
disp(sprintf('Bayesian Random Effect Variance = %1.4f', full(v)));

subplot(2,1,2)
plot(x,p,[1 1]*T,[0 max(p)],'-.');
title({'Bayesian Random Effects: Posterior density of contrast',...
        sprintf('P(contrast > %0.2f) = %.1f%s',T,PP*100,'%')},...
       'FontSize',12)
xlabel('contrast');
ylabel('probability density');

i    = find(x >= T);
hold on
fill([x(i) fliplr(x(i))],[i*0 fliplr(p(i))],[1 1 1]*.8)
axis square, grid on
hold off

% Get Classical RFX p-value
% Note: 1-p is equal to Bayesian RFX if T=0
%--------------------------------------------------------------------------
t  = (sqrt(num_models)*c)/sqrt(v);
p  = 1 - spm_Tcdf(t,num_models - 1);
disp(sprintf('Classical Random Effects p-value = %1.4f', p));
disp('Note: 1-p is equal to Bayesian RFX if threshold is zero');
spm_input('Thank you',1,'d');


function [] = spm_dcm_sessions ()
% Apply contrast vector to multiple DCM models
%
% FORMAT [] = spm_dcm_sessions ()
%

Finter = spm_figure('GetWin','Interactive');
Fgraph = spm_figure('GetWin','Graphics');

num_models = spm_input('Apply contrast to how many DCM models ? ','+1','r',[],1);
P     = spm_get(num_models,'DCM*.mat',{'select DCM*.mat files'});

% Get contrast
str     = 'contrast for';
D       = spm_input(str,1,'b',{'A','B','C'});
load(P{1});
con=spm_dcm_contrasts(P{1},D);

% Get threshold
str = 'Threshold';
DCM.T      = spm_input(str,1,'e',0,[1 1]);
T=DCM.T;

% Get mean and variance of effect in each model
for model=1:num_models,
    load(P{model});
    
    l     = DCM.M.l;
    m     = DCM.M.m;
    i       = find(con); 
    switch D,
    case 'A',
        j = 1;
    case 'B',
        j = 1 + l*l;
    case 'C',
        j = 1 + l*l + l*l*m;
    end
    C       = sparse(i + j,1,con(i),length(DCM.Ep),1);
    mean_con(model) = C'*DCM.Ep;
    var_con(model) = C'*DCM.Cp*C;   
    disp(sprintf('Model %d, Contrast Mean=%1.4f, Contrast Variance=%1.4f',model,mean_con(model),var_con(model)));
end


%- Bayesian fixed effects analysis 
%-----------------------------------------------------------
precision_con=1./var_con;
precision=sum(precision_con);
v=1/precision;
c   = sum((precision_con.*mean_con)/precision);
x    = c + [-32:32]*sqrt(v)*6/32;
p    = 1/sqrt(2*pi*v)*exp(-[x - c].^2/(2*v));
PP   = 1 - spm_Ncdf(T,c,v);
disp(' ');
disp(sprintf('Bayesian Fixed Effect Mean = %1.4f', full(c)));
disp(sprintf('Bayesian Fixed Effect Variance = %1.4f', full(v)));

figure(Fgraph)
subplot(2,1,1)
plot(x,p,[1 1]*T,[0 max(p)],'-.');
title({'Bayesian Fixed Effects: Posterior density of contrast',...
        sprintf('P(contrast > %0.2f) = %.1f%s',T,PP*100,'%')},...
    'FontSize',12)
xlabel('contrast')
ylabel('probability density')

i    = find(x >= T);
hold on
fill([x(i) fliplr(x(i))],[i*0 fliplr(p(i))],[1 1 1]*.8)
axis square, grid on
hold off

% Random effects analysis
v = std(mean_con)^2;
c    = mean(mean_con);
x    = c + [-32:32]*sqrt(v)*6/32;
p    = 1/sqrt(2*pi*v)*exp(-[x - c].^2/(2*v));
PP   = 1 - spm_Ncdf(T,c,v);
disp(' ');
disp(sprintf('Bayesian Random Effect Mean = %1.4f', full(c)));
disp(sprintf('Bayesian Random Effect Variance = %1.4f', full(v)));

subplot(2,1,2)
plot(x,p,[1 1]*T,[0 max(p)],'-.');
title({'Bayesian Random Effects: Posterior density of contrast',...
        sprintf('P(contrast > %0.2f) = %.1f%s',T,PP*100,'%')},...
    'FontSize',12)
xlabel('contrast')
ylabel('probability density')

i    = find(x >= T);
hold on
fill([x(i) fliplr(x(i))],[i*0 fliplr(p(i))],[1 1 1]*.8)
axis square, grid on
hold off


t=c/sqrt(v);
p= 1 - spm_tcdf(t,num_models-1);
disp(sprintf('Classical Random Effects p-value = %1.4f', p));

spm_input('Thank you',1,'d');

        
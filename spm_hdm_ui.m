function [Ep,Cp,K1,K2] = spm_hdm_ui(xSPM,SPM,hReg)
% user interface for hemodynamic model estimation
% FORMAT [Ep,Cp,K1,K2] = spm_hdm_ui(xSPM,SPM,hReg);
%
% xSPM   - structure containing specific SPM details
% SPM    - structure containing generic  SPM details
% hReg   - Handle of results section XYZ registry (see spm_results_ui.m)
%
% Ep     - conditional expectations of the hemodynamic model parameters
% Cp     - conditional  covariance  of the hemodynamic model parameters
% K1     - 1st order kernels
% K2     - 2nd order kernels
%          (see main body of routine for details of model specification)
%___________________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% Karl Friston
% $Id: spm_hdm_ui.m 4048 2010-08-26 16:29:36Z guillaume $


% get figure handles
%---------------------------------------------------------------------------
Finter = spm_figure('GetWin','Interactive');
header = get(Finter,'Name');
set(Finter,'Name','Hemodynamic modelling')

% inputs
%===========================================================================

% which session?
%---------------------------------------------------------------------------
s    = length(SPM.Sess);
if s > 1
    s = spm_input('which session',1,'n1',s,s);
end
Sess = SPM.Sess(s);

% 'causes' or imputs U
%---------------------------------------------------------------------------
spm_input('Input specification:...  ',1,'d');
U.dt = Sess.U(1).dt;
u    = length(Sess.U);
if u == 1 && length(Sess.U(1).name) == 1
    U.name = Sess.U(1).name;
    U.u    = Sess.U(1).u(33:end,1);
else
    U.name = {};
    U.u    = [];
    for  i = 1:u
    for  j = 1:length(Sess.U(i).name)
        str   = ['include ' Sess.U(i).name{j} '?'];
        if spm_input(str,2,'y/n',[1 0])
            U.u             = [U.u Sess.U(i).u(33:end,j)];
            U.name{end + 1} = Sess.U(i).name{j};
        end
    end
    end
end

%-Echo time (TE) of data acquisition
%--------------------------------------------------------------------------
TE    = 0.04;
TE_ok = 0;
while ~TE_ok
    TE = spm_input('Echo time, TE [s]', '+1', 'r', TE);
    if ~TE || (TE < 0) || (TE > 0.1)
        str = { 'Extreme value for TE or TE undefined.',...
            'Please re-enter TE (in seconds!)'};
        spm_input(str,'+1','bd','OK',[1],1);
    else
        TE_ok = 1;
    end
end

    
%-System outputs
%===========================================================================

% enforce adjustment w.r.t. all effects
%---------------------------------------------------------------------------
xY     = struct(    'Ic'        ,0,...  
                    'name'      ,'HDM',...
                    'Sess'      ,s);

% get region stucture
%---------------------------------------------------------------------------
xY     = struct('name', 'HDM', 'Sess', s);
[y xY] = spm_regions(xSPM,SPM,hReg,xY);

% place response and confounds in response structure
%---------------------------------------------------------------------------
y      = xY.u;
Y.y    = y;
Y.dt   = SPM.xY.RT;
Y.X0   = xY.X0;

%-Estimate
%===========================================================================
spm('Pointer','Watch')
spm('FigName','Estimation in progress');


% Model specification: m input; 4 states; 1 outout; m + 6 parameters
%---------------------------------------------------------------------------
% u(m) - mth stimulus function     (u)
%
% x(1) - vascular signal           log(s)
% x(2) - rCBF                      log(f)
% x(3) - venous volume             log(v)
% x(4) - deoyxHb                   log(q)
%
% y(1) - BOLD                      (y)
%
% P(1)       - signal decay               d(ds/dt)/ds)      half-life = log(2)/P(1) ~ 1sec
% P(2)       - autoregulation             d(ds/dt)/df)      2*pi*sqrt(1/P(1)) ~ 10 sec
% P(3)       - transit time               (t0)              ~ 1 sec
% P(4)       - exponent for Fout(v)       (alpha)           c.f. Grubb's exponent (~ 0.38)
% P(5)       - resting oxygen extraction  (E0)              ~ range 20 - 50%
% P(6)       - ratio of intra- to extra-  (epsilon)         ~ range 0.5 - 2
%              vascular components   
%              of the gradient echo signal   

% P(6 + 1:m) - input efficacies - d(ds/dt)/du)  ~0.3 per event
%--------------------------------------------------------------------------

% priors (3 modes of hemodynamic variation)
%--------------------------------------------------------------------------
m       = size(U.u,2);
[pE,pC] = spm_hdm_priors(m,3);

% model
%--------------------------------------------------------------------------
M.f     = 'spm_fx_hdm';
M.g     = 'spm_gx_hdm';
M.x     = [0 0 0 0]'; 
M.pE    = pE;    
M.pC    = pC;
M.m     = m;
M.n     = 4;
M.l     = 1;
M.N     = 64;
M.dt    = 24/M.N;
M.TE    = TE;

% nonlinear system identification
%--------------------------------------------------------------------------
[Ep,Cp,Ce,K0,K1,K2,M0,M1] = spm_nlsi(M,U,Y);

%-display results
%==========================================================================
t       = [1:M.N]*M.dt;
Fhdm    = spm_figure;
set(Fhdm,'name','Hemodynamic Modeling')


% display input parameters
%--------------------------------------------------------------------------
subplot(2,2,1)
P     = Ep(7:end);
C     = diag(Cp(7:end,7:end));
[i j] = max(abs(P));
spm_barh(P,C)
axis square
title({'stimulus efficacy'; 'with 90% confidence intervals'},'FontSize',10)
set(gca,'Ytick',[1:m],'YTickLabel',U.name,'FontSize',8)
str = {};
for i = 1:m
    str{end + 1} = U.name{i};
    str{end + 1} = sprintf('mean = %0.2f',P(i));
    str{end + 1} = '';
end
set(gca,'Ytick',[1:m*3]/3 + 1/2,'YTickLabel',str)
xlabel('relative efficacy per event/sec')


% display hemodynamic parameters
%---------------------------------------------------------------------------
subplot(2,2,3)
P     = Ep(1:6);
pE    = pE(1:6);
C     = diag(Cp(1:6,1:6));
spm_barh(P,C,pE)
title({ 'hemodynamic parameters'},'FontSize',10)
set(gca,'Ytick',[1:18]/3 + 1/2)
set(gca,'YTickLabel',{  'SIGNAL decay',...
             sprintf('%0.2f per sec',P(1)),'',...
            'FEEDBACK',...
             sprintf('%0.2f per sec',P(2)),'',...
            'TRANSIT TIME',...
             sprintf('%0.2f seconds',P(3)),'',...
            'EXPONENT',...
             sprintf('%0.2f',P(4)),'',...
            'EXTRACTION',...
             sprintf('%0.0f %s',P(5)*100,'%'),'',...
            'log SIGNAL RATIO',...
             sprintf('%0.2f %s',P(6),'%'),''},'FontSize',8)


% get display state kernels (i.e. state dynamics) 
%==========================================================================

% Volterra kernels of states
%--------------------------------------------------------------------------
[H0,H1] = spm_kernels(M0,M1,M.N,M.dt);

subplot(3,2,2)
plot(t,exp(H1(:,:,j)))
axis square
title({['1st order kernels for ' U.name{j}];...
       'state variables'},'FontSize',9)
ylabel('normalized values')
legend('s','f','v','q',0);
grid on


% display output kernels (i.e. BOLD response) 
%--------------------------------------------------------------------------
subplot(3,2,4)
plot(t,K1(:,:,j))
axis square
title({'1st order kernel';...
       'output: BOLD'},'FontSize',9)
ylabel('normalized flow signal')
grid on

subplot(3,2,6)
imagesc(t,t,K2(:,:,1,j,j))
axis square
title({'2nd order kernel';...
       'output: BOLD'},'FontSize',9)
xlabel({'time {seconds} for'; U.name{j}})
grid on


%-Reset title
%--------------------------------------------------------------------------
spm('FigName',header);
spm('Pointer','Arrow')
spm_input('Thank you',1,'d')

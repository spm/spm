function [Ep,Cp,Ce] = spm_hdm_ui(SPM,VOL,xX,xCon,xSDM,hReg)
% user interface for hemodynamic model estimation
% FORMAT [Ep Cp Ce] = spm_hdm_ui(SPM,VOL,xX,xCon,xSDM,hReg);
%
% SPM    - structure containing SPM, distribution & filtering detals
% VOL    - structure containing details of volume analysed
% xX     - Design Matrix structure
% xSDM   - structure containing contents of SPM.mat file
% xCon   - Contrast definitions structure (see spm_FcUtil.m for details)
% hReg   - Handle of results section XYZ registry (see spm_results_ui.m)
%
% Ep     - conditional expectations of the hemodynamic model parameters
% Cp     - conditional  covariance  of the hemodynamic model parameters
%          (see main body of routine for details of model specification)
% Ce     - ML estimate of error covariance
%___________________________________________________________________________
% %W% Karl Friston %E%

% get figure handles
%---------------------------------------------------------------------------
Finter = spm_figure('GetWin','Interactive');
header = get(Finter,'Name');
set(Finter,'Name','Hemodynamic modelling')


% inputs
%===========================================================================

% which session?
%---------------------------------------------------------------------------
s      = length(xSDM.Sess);
if s > 1
	s = spm_input('which session','+1','n',1,[1 1],s);
end
Sess   = xSDM.Sess{s};

% 'causes' or imputs U
%---------------------------------------------------------------------------
U.u    = [];
U.dt   = xX.dt;
n      = length(Sess.sf{1});
e      = spm_input('trial types?','+1','b','event|epoch',[1 0]);
d      = 0;
for  i = 1:length(Sess.sf)
	trial       = Sess.name{i};
	if e
		u   = Sess.sf{i}/U.dt;
	else
		str = sprintf('duration (scans): %s',trial);
		d   = spm_input(str,'+0','w',d,[1 1]);
		q   = max([1 d*xX.RT/xX.dt]);
		u   = full(Sess.sf{i});
		u   = conv(u,ones(q,1));
		u   = sparse(u(1:n));
	end
	U.u         = [U.u u];
	U.name{i}   = trial;
end


% system outputs
%===========================================================================

% enforce adjustment w.r.t. all effects high pass filtering
%---------------------------------------------------------------------------
xY     = struct(	'Ic'		,1,...	
 	   		'name'		,'HDM',...
 	   		'filter'	,'high',...
 	   		'Sess'		,s);


% get region stucture
%---------------------------------------------------------------------------
[y xY] = spm_regions(SPM,VOL,xX,xCon,xSDM,hReg,xY);


% serial correlations?
%---------------------------------------------------------------------------
h      = spm_input('model serial correlations','+1','b','yes|no',[2 1]);

%-append low frequency confounds
%---------------------------------------------------------------------------
n      = size(xY.y,1);
X0     = ones(n,1);
if ~strcmp(xY.filter,'none')
	X0   = [X0 full(xX.K{s}.KH)];
end

%-adjust and place in response variable structure
%---------------------------------------------------------------------------
y      = xY.u;
y      = y - X0*(pinv(X0)*y);
Y.y    = y;
Y.dt   = xX.RT;
Y.X0   = X0;
Y.Ce   = spm_Ce(n,h);


% estimate
%===========================================================================
spm('Pointer','Watch')
spm('FigName','Estimation in progress');


% Model specification: m input; 4 state; 1 outout; m + 5 parameters
%---------------------------------------------------------------------------
% u(m) - mth stimulus function     (u)
%
% x(1) - vascular signal           (s)
% x(2) - rCBF                      (f)
% x(3) - venous volume             (v)
% x(4) - deoyxHb                   (q)
%
% y(1) - BOLD                      (y)
%
% P(1)       - signal decay     - d(ds/dt)/ds)  half-life = log(2)/P(2) ~ 1sec
% P(2)       - autoregulation   - d(ds/dt)/df)  2*pi*sqrt(1/P(3)) ~ 10 sec
% P(3)       - transit time               (t0)  ~ 1 sec
% P(4)       - exponent for Fout(v)    (alpha)  c.f. Grubb's exponent (~ 0.38)
% P(5)       - resting oxygen extraction  (E0)  ~ range 20 - 50%
% P(5 + 1:m) - input efficacies - d(ds/dt)/du)  ~0.3 per event
%---------------------------------------------------------------------------

% priors
%---------------------------------------------------------------------------
m       = size(U.u,2);
[pE,pC] = spm_hdm_priors(m);

% model
%---------------------------------------------------------------------------
M.fx    = 'spm_fx_HRF';
M.lx    = 'spm_lambda_HRF';
M.x     = [0 1 1 1]';
M.pE    = pE;    
M.pC    = pC;
M.pD    = (Y.dt/4)^2;
M.m     = m;
M.n     = 4;
M.l     = 1;
M.N     = 64;
M.dt    = 24/M.N;

% nonlinear system identification
%---------------------------------------------------------------------------
[Ep,Cp,Ce,K0,K1,K2,M0,M1,L] = spm_nlsi(M,U,Y);

%-display results
%===========================================================================
t       = [1:M.N]*M.dt;
Fhdm    = spm_figure;
set(Fhdm,'name','Hemodynamic Modeling')


% display input parameters
%---------------------------------------------------------------------------
subplot(2,2,1)
P     = Ep(6:end);
C     = diag(Cp(6:end,6:end));
[i j] = max(P);
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
P     = Ep(1:5);
C     = diag(Cp(1:5,1:5));
spm_barh(P,C)
title({	'hemodynamic parameters'},'FontSize',10)
set(gca,'Ytick',[1:15]/3 + 1/2)
set(gca,'YTickLabel',{	'SIGNAL decay',...
			 sprintf('%0.2f per sec',P(1)),'',...
			'FEEDBACK',...
			 sprintf('%0.2f per sec',P(2)),'',...
			'TRANSIT TIME',...
			 sprintf('%0.2f seconds',P(3)),'',...
			'EXPONENT',...
			 sprintf('%0.2f',P(4)),'',...
			'EXTRACTION',...
			 sprintf('%0.0f %s',P(5)*100,'%'),''},'FontSize',8)


% get display state kernels (i.e. state dynamics) 
%===========================================================================

% Bilinear representation (2nd order)
%---------------------------------------------------------------------------
[M0,M1,L] = spm_bi_reduce(M,Ep);
L         = sparse([1:M.n],[1:M.n] + 1,1,M.n,length(M0));

% Volterra kernels of states
%---------------------------------------------------------------------------
[H0,H1] = spm_kernels(M0,M1,L,M.N,M.dt);

subplot(3,2,2)
plot(t,H1(:,:,j))
axis square
title({	['1st order kernels for ' U.name{j}];...
	'state variables: s, f, q and v'},'FontSize',9)
ylabel('normalized values')
grid on


% display output kernels (i.e. BOLD response) 
%---------------------------------------------------------------------------
subplot(3,2,4)
plot(t,K1(:,:,j))
axis square
title({	'1st order kernel';...
	'output: BOLD'},'FontSize',9)
ylabel('normalized flow signal')
grid on

subplot(3,2,6)
imagesc(t,t,K2(:,:,1,j,j))
axis square
title({	'2nd order kernel';...
	'output: BOLD'},'FontSize',9)
xlabel({'time {seconds} for'; U.name{j}})
grid on


%-Reset title
%---------------------------------------------------------------------------
spm('FigName',header);
spm('Pointer','Arrow')
spm_input('Thank you',1,'d')

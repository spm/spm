function spm_dcm_ui(SPM,VOL,xX,xCon,xSDM,hReg,Action)
% user interface for Dynamic Causal Modeling (DCM)
% FORMAT spm_dcm_ui(SPM,VOL,xX,xCon,xSDM,hReg,Action);
%
% SPM    - structure containing SPM, distribution & filtering detals
% VOL    - structure containing details of volume analysed
% xX     - Design Matrix structure
% xSDM   - structure containing contents of SPM.mat file
% xCon   - Contrast definitions structure (see spm_FcUtil.m for details)
% hReg   - Handle of results section XYZ registry (see spm_results_ui.m)
%
% saves DCM.mat
%
%	DCM.M      - model  specification structure (see spm_nlsi)
%	DCM.Y      - output specification structure (see spm_nlsi)
%	DCM.U      - input  specification structure (see spm_nlsi)
%	DCM.A      - intrinsic connection matrix
%	DCM.B      - input-dependent connection matrix
%	DCM.C      - nput connection matrix
%	DCM.pA     - pA - posterior probabilities
%	DCM.pB     - pB - posterior probabilities
%	DCM.pC     - pC - posterior probabilities
%	DCM.H1     - 1st order Volterra Kernels - neuronal
%	DCM.K1     - 1st order Volterra Kernels - hemodynamic
%	DCM.xh     - time points - kernels
%	DCM.xy     - time points - ouputs
%	DCM.xu     - time points - inputs
%	DCM.R      - residuals
%	DCM.y      - predicte responses
%	DCM.xY     - original response variable structures
%	DCM.T      - threshold for inference based on posterior p.d.f
%
%___________________________________________________________________________
% %W% Karl Friston %E%

% get figure handles
%---------------------------------------------------------------------------
Finter = spm_figure('GetWin','Interactive');
Fgraph = spm_figure('GetWin','Graphics');
header = get(Finter,'Name');
set(Finter,'Name','Dynamic Causal Modelling')
WS     = spm('WinScale');
global SWD


% options
%---------------------------------------------------------------------------
if nargin < 7
	str    = 'specify or review an analysis';
	if spm_input(str,1,'b',{'new','results'},[0 1]);
		Action = 'results';
	else
		Action = spm_input('specify',1,'b',{'VOIs','graph'});
	end
end


switch Action

% specify VOIs
%---------------------------------------------------------------------------
case 'VOIs'

	% display volume maxima
	%-------------------------------------------------------------------
	spm_list('List',SPM,VOL,[],[],'',hReg);

	% which session?
	%-------------------------------------------------------------------
	s         = length(xSDM.Sess);
	if s > 1
		s = spm_input('which session','+1','r',1,1,s);
	end

	% get VOIs
	%-------------------------------------------------------------------
	str   = ['xY = struct(	''Ic'',		 1,',...
 	   			'''filter'',	''high'',',...
 	   			'''Sess'',' 	 sprintf('%i',s) ',',...
 	      			'''def'',	''sphere'',',...
  	    			'''spec'',	 4);',...
		'spm_regions(SPM,VOL,xX,xCon,xSDM,hReg,xY);'];

	h(1)  = uicontrol(Finter,'String','save region',...
		'Position',[040 200 240 020].*WS,...
		'CallBack',str);
	str   = ['delete(get(gco,''Userdata'')),',...
		 'spm_input(''Thank you'',1,''d''),return'];
	h(2)  = uicontrol(Finter,'String','done',...
		'Position',[300 200 060 020].*WS,...
		'CallBack',str);
	set(h(2),'UserData',h);


% specify graph and estimate
%---------------------------------------------------------------------------
case 'graph'

	% name
	%===================================================================
	name  = spm_input('name for DCM_???.mat','+1','s');


	% outputs
	%===================================================================

	% get cell array of region structures
	%-------------------------------------------------------------------
	clear xY
	P     = spm_get([2 8],'VOI*.mat',{'select VOIs'});
	m     = length(P);
	for i = 1:m
		p     = load(P{i},'xY');
		xY(i) = p.xY;
	end

	% inputs
	%===================================================================

	% get 'causes' or imputs U
	%-------------------------------------------------------------------
	Sess  = xSDM.Sess{xY.Sess};
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


	% graph connections
	%===================================================================
	n     = size(U.u,2);
	a     = zeros(m,m);
	c     = zeros(m,n);
	b     = zeros(m,m,n);
	d     = uicontrol(Finter,'String','done',...
			'Position',[300 200 060 020].*WS);
	dx    = 20;


	%-intrinsic connections
	%-------------------------------------------------------------------
	spm_input('Specify intrinsic connections from',1,'d')
	spm_input('to',3,'d')

	for i = 1:m
		str    = sprintf('%s %i',xY(i).name,i);
		h1(i)  = uicontrol(Finter,'String',str,...
				'Style','text',...
				'HorizontalAlignment','right',...
				'Position',[080 356-dx*i 080 020].*WS);
		h2(i)  = uicontrol(Finter,'String',sprintf('%i',i),...
				'Style','text',...
				'Position',[180+dx*i 356 020 020].*WS);
	end
	for i = 1:m
		for j = 1:m
			h3(i,j) = uicontrol(Finter,...
				'Position',[180+dx*j 360-dx*i 020 020].*WS,...
				'Style','radiobutton');
			if i == j
				set(h3(i,j),'Value',1,...
					'enable','off',...
					'BackgroundColor',[0.5 0.5 0.5]);
			end
			
		end
	end
	drawnow

	% wait for 'done'
	%-----------------------------------------------------------
	while(1)
		waitforbuttonpress;
		if strcmp(get(gco,'Type'),'uicontrol')
			if strcmp(get(gco,'String'),'done')
				for i = 1:m
				for j = 1:m
					a(i,j) = get(h3(i,j),'Value');
				end
				end
				delete([h1(:); h2(:); h3(:)])
				break
			end
		end
	end


	%-effects of causes
	%-------------------------------------------------------------------
	for k = 1:n

		% buttons and labels
		%-----------------------------------------------------------
		str   = sprintf(...
			'Effects of %-12s on regions... and connections',...
			 U.name{k});
		spm_input(str,1,'d')

		dx    = 20;
		for i = 1:m
			h1(i)  = uicontrol(Finter,'String',xY(i).name,...
				'Style','text',...
				'Position',[080 356-dx*i 080 020].*WS);
			h2(i)  = uicontrol(Finter,...
				'Position',[160 360-dx*i 020 020].*WS,...
				'Style','radiobutton');
		end
		for i = 1:m
		    for j = 1:m
			h3(i,j) = uicontrol(Finter,...
				'Position',[220+dx*j 360-dx*i 020 020].*WS,...
				'Style','radiobutton');
		    end
		end
		drawnow

		% wait for 'done'
		%-----------------------------------------------------------
		while(1)
			waitforbuttonpress;
			if strcmp(get(gco,'Type'),'uicontrol')
			if strcmp(get(gco,'String'),'done')

			% get c
			%--------------------------------------------------
			for i = 1:m
				c(i,k)   = get(h2(i)  ,'Value');
			end

			% get b ensuring 2nd order effects are allowed
			%--------------------------------------------------
			for i = 1:m
				for j = 1:m
					b(i,j,k) = get(h3(i,j),'Value');
					if i == j & ~c(i,k)
						b(i,j,k) = 0;
					end
				end
			end
			delete([h1(:); h2(:); h3(:)])
			spm_input('Thank you',1,'d')
			break

			end
			end
		end
	end
	delete(d)

	% Model specification:
	%===================================================================
	spm('Pointer','Watch')
	spm('FigName','Estimation in progress');

	% confounds
	%-------------------------------------------------------------------
	v     = size(xY(1).u,1);
	X0    = [ones(v,1) full(xX.K{xY(1).Sess}.KH)];


	% response variables
	%-------------------------------------------------------------------
	n     = length(xY);
	for i = 1:n
		Y.y(:,i)    = xY(i).u;
		Y.name{i}   = xY(i).name;
	end
	Y.dt  = xX.RT;
	Y.X0  = X0;
	Y.Ce  = spm_Ce(v*ones(1,n));


	% input variables
	%-------------------------------------------------------------------
	[u m] = size(U.u);


	% prior contraints - expectations
	%-------------------------------------------------------------------
	A     = -eye(n)*log(2)/0.5;		% 500ms neuronal half life
	B     = zeros(n,n,m);
	C     = zeros(n,m);
	D     = zeros(n,1);

	% prior contraints - variances
	%-------------------------------------------------------------------
	T     = 1/4;				% threshold T1/2 = log(2)/T
	q     = (T/2)^2;
	a     = a - diag(diag(a));
	a     = a*q;
	b     = b*q;
	c     = c*q;
	d     = zeros(n,1);

	% model specification and nonlinear system identification
	%-------------------------------------------------------------------
	M.fx  = 'spm_fx_dcm';
	M.lx  = 'spm_lx_dcm';
	M.x   = sparse(n*5,1);
	M.pE  =      [A(:); B(:); C(:); D(:)];
	M.pC  = diag([a(:); b(:); c(:); d(:)]);
	M.m   = size(U.u,2);
	M.n   = size(M.x,1);
	M.l   = size(Y.y,2);
	M.N   = 128;
	M.dt  = 16/M.N;
	[Ep,Cp,Ce,H0,H1,H2,K0,K1,K2] = spm_nlsi(M,U,Y);


	% predicted responses
	%-------------------------------------------------------------------
	y     = spm_nlsi_int(Ep,M,U,v);
	R     = Y.y - y;
	R     = R - X0*(pinv(X0)*(R));


	% assemble model estimation structure
	%-------------------------------------------------------------------
	pp         = 1 - spm_Ncdf(M.pE + T, Ep,diag(Cp));
	qp         = 1 - spm_Ncdf(M.pE + T,-Ep,diag(Cp));
	pp         = max([pp qp],[],2);
	[ A  B  C] = spm_bl_reshape(Ep,m,n,0);
	[pA pB pC] = spm_bl_reshape(pp,m,n,0);

	DCM.M      = M;
	DCM.Y      = Y;
	DCM.U      = U;
	DCM.A      = A;
	DCM.B      = B;
	DCM.C      = C;
	DCM.pA     = pA;
	DCM.pB     = pB;
	DCM.pC     = pC;
	DCM.H1     = H1(:,[1:n],:);
	DCM.K1     = K1(:,[1:n],:);
	DCM.xh     = [1:M.N]*M.dt;
	DCM.xy     = [1:v]*Y.dt;
	DCM.xu     = [1:u]*U.dt;
	DCM.R      = R;
	DCM.y      = y;
	DCM.xY     = xY;
	DCM.T      = T;


	%-Save and reset title
	%-------------------------------------------------------------------
	save(fullfile(SPM.swd,['DCM_' name]),'DCM');
	spm('FigName',header);
	spm('Pointer','Arrow')
	return


% review results
%---------------------------------------------------------------------------
case 'results'

	%-display model details
	%-------------------------------------------------------------------
	Fdcm  = spm_figure;
	set(Fdcm,'name','Dynamic Causal Modeling')

	%-get results
	%-------------------------------------------------------------------
	P     = spm_get(1,'DCM*.mat',{'select DCM_???.mat'});
	load(P{:})
	m     = size(DCM.U.u,2);
	l     = size(DCM.Y.y,2);

	while(1)
	% get options
	%-------------------------------------------------------------------
	str   = {	'Selected regions',...
			'Inputs and outputs',...
			'Predicted outputs',...
			'Connectivity matrices',...
			'Kernels',...
			'Intrinsic connections'};
	for i = 1:m
		str{end + 1} = ['Effects of ' DCM.U.name{i}];
	end
	str{end + 1} = 'quit';
	OPT   = spm_input('display',1,'m',str);
	figure(Fdcm)
	spm_figure('Clear',Fdcm)

 	% Selected regions
	%-------------------------------------------------------------------
	if OPT == 1


		% transverse
		%-----------------------------------------------------------
		subplot(2,2,3)
		title('transverse')
		u = [[0 1 0];[-1 0 0]]';
		spm_dcm_display(DCM.xY,[],[],[],u,32);

		% sagittal
		%-----------------------------------------------------------
		subplot(2,2,1)
		title('sagittal')
		u = [[0 1 0];[0 0 1]]';
		spm_dcm_display(DCM.xY,[],[],[],u,32);

		% coronal
		%-----------------------------------------------------------
		subplot(2,2,2)
		title('coronal')
		u = [[1 0 0];[0 0 1]]';
		spm_dcm_display(DCM.xY,[],[],[],u,32);

		% table
		%-----------------------------------------------------------
		subplot(2,2,4)
		title(str{OPT},'FontSize',12)
		y = 0;
		line([0 4],[y y])
		y = y - 1;
		text(0.0,y,['Name'],    'FontSize',10)
		text(1.0,y,['Voxels'],  'FontSize',10)
		text(2.0,y,['Location (mm)'],'FontSize',10)
		y = y - 1;
		line([0 4],[y y],'LineWidth',4)
		y = y - 1;
		for i = 1:length(DCM.xY)
			name = DCM.xY(i).name;
			N    = length(DCM.xY(i).s);
			L    = DCM.xY(i).xyz;
			r    = DCM.xY(i).spec;
			text(0,y,name,                      	   'FontSize',8)
			text(1,y,sprintf('%0.0f',N),               'FontSize',8)
			text(2,y,sprintf('%-4.0f %-4.0f %-4.0f',L),'FontSize',8)
			y = y - 1;
		end
		line([0 4],[y y])
		axis off square


 	% inputs and outputs
	%-------------------------------------------------------------------
	elseif OPT == 2

		% assemble graphics data
		%-----------------------------------------------------------
		for i = 1:length(DCM.xY)
			G(i).y  = DCM.y(:,i) + DCM.R(:,i);
			G(i).x  = DCM.xy;
		end

		% inputs
		%-----------------------------------------------------------
		for i = 1:m
			subplot(2*m,1,i)
			plot(DCM.xu,DCM.U.u(:,i))
			title(['Inputs - ' DCM.U.name{i}])
			axis tight
			set(gca,'FontSize',10)
			ylabel('event density {Hz}')
		end
		xlabel('time {seconds}')

		% outputs
		%-----------------------------------------------------------
		subplot(2,1,2)
		title('Outputs')
		spm_dcm_display(DCM.xY,[],[],G);


	% Predicted outputs
	%-------------------------------------------------------------------
	elseif OPT == 3

		% assemble graphics data
		%-----------------------------------------------------------
		for i = 1:length(DCM.xY)
			G(i).y  = [DCM.y(:,i) DCM.y(:,i) + DCM.R(:,i)];
			G(i).x  = DCM.xy;
		end

		% outputs
		%-----------------------------------------------------------
		subplot(1,1,1)
		title(str{OPT})
		spm_dcm_display(DCM.xY,[],[],G,[],24,[],1/4);


	% Connectivity matrices
	%-------------------------------------------------------------------
	elseif OPT == 4

		% intrinsic interactions
		%-----------------------------------------------------------
		subplot(4,2,1)
		imagesc(DCM.A)
		axis image
		set(gca,'YTick',[1:l],'YTickLabel',DCM.Y.name)
		set(gca,'XTick',[1:l],'XTickLabel',DCM.Y.name)
		title('Intrinsic connections','FontSize',10)
		set(gca,'FontSize',10)
		ylabel('to')

		% intrinsic interactions - probabilities
		%-----------------------------------------------------------
		subplot(4,2,2)
		bar(DCM.pA)
		set(gca,'XTick',[1:l],'XTickLabel',DCM.Y.name)
		title(sprintf('P(A > %0.2f)',DCM.T),'FontSize',10)
		set(gca,'FontSize',10)
		axis square

		% extrinsic effects
		%-----------------------------------------------------------
		subplot(4,2,3)
		imagesc(DCM.C)
		axis image
		set(gca,'YTick',[1:l],'YTickLabel',DCM.Y.name)
		set(gca,'XTick',[1:m],'XTickLabel',DCM.U.name)
		title('Extrinsic connections','FontSize',10)
		set(gca,'FontSize',10)
		xlabel('from')
		ylabel('to')

		% intrinsic interactions - probabilities
		%-----------------------------------------------------------
		subplot(4,2,4)
		bar(DCM.pC)
		set(gca,'XTick',[1:l],'XTickLabel',DCM.Y.name)
		title(sprintf('P(C > %0.2f)',DCM.T),'FontSize',10)
		set(gca,'FontSize',10)
		xlabel('recipient region')
		axis square

		% input effects
		%-----------------------------------------------------------
		for i = 1:m

			subplot(2*m,2,2*(i - 1) + 2*m + 1)
			imagesc(DCM.B(:,:,i))
			set(gca,'YTick',[1:l],'YTickLabel',DCM.Y.name)
			ylabel(['Effect of ' DCM.U.name{i}])
			axis image

			subplot(2*m,2,2*(i - 1) + 2*m + 2)
			bar(DCM.pB(:,:,i))
			set(gca,'XTick',[1:l],'XTickLabel',DCM.Y.name)
			ylabel(sprintf('P(B > %0.2f)',DCM.T))
			set(gca,'FontSize',10)
			axis square
		end
		xlabel('recipient region')


	% Kernels
	%-------------------------------------------------------------------
	elseif OPT == 5

		% input effects
		%-----------------------------------------------------------
		for i = 1:m

			% input effects - neuronal
			%---------------------------------------------------
			x     = DCM.xh;
			y     = DCM.H1(:,:,i);
			subplot(m,3,3*(i - 1) + 1)
			plot(x,y)
			set(gca,'XLim',[0 8])
			axis square
			title(['neuronal responses to ' DCM.U.name{i}])
			grid on
			xlabel('time {seconds}')
			for j = 1:l
				text(x(j*8),y(j*8,j),DCM.Y.name{j},...
					'FontSize',8,...
					'HorizontalAlignment','Center')
			end

			% input effects - simulated EEG
			%---------------------------------------------------
			y     = gradient(y);
			subplot(m,3,3*(i - 1) + 2)
			plot(x,y)
			set(gca,'XLim',[0 5])
			axis square ij
			title('simulated ERP')
			grid on
			xlabel('time {ms}')
			for j = 1:l
				text(x(j*8),y(j*8,j),DCM.Y.name{j},...
					'FontSize',8,...
					'HorizontalAlignment','Center')
			end

			% input effects - hemodynamic
			%---------------------------------------------------
			x     = DCM.xh;
			y     = DCM.H1(:,:,i);
			k     = DCM.K1(:,:,i);
			subplot(m,3,3*(i - 1) + 3)
			plot(x,k,x,y,':')
			set(gca,'XLim',[0 16])
			axis square
			title('hemodynamic responses')
			grid on
			xlabel('time {seconds}')
			for j = 1:l
				text(x(j*16),k(j*16,j),DCM.Y.name{j},...
					'FontSize',8,...
					'HorizontalAlignment','Center')
			end
		end



	% Intrinsic connections
	%-------------------------------------------------------------------
	elseif OPT == 6

		% Intrinsic effects
		%-----------------------------------------------------------
		a        = DCM.pA;
		a(:,:,2) = DCM.A;
		spm_dcm_display(DCM.xY,a)
		title({	str{OPT};...
			sprintf('P(|connection| > %0.2f)',DCM.T);...
			'connection strength'},'FontSize',10)

	% Intrinsic connections
	%-------------------------------------------------------------------
	elseif OPT >= 7 & OPT < 7 + m

		% input effects
		%-----------------------------------------------------------
		i        = OPT - 6;
		b        = DCM.pB(:,:,i);
		b(:,:,2) = DCM.B(:,:,i);
		c        = DCM.pC(:,i);
		c(:,2)   = DCM.C(:,i);
		spm_dcm_display(DCM.xY,b,c)
		title({	str{OPT};...
			sprintf('P(|connection| > %0.2f)',DCM.T);...
			'connection strength'},'FontSize',10)

	% quit
	%-------------------------------------------------------------------
	elseif OPT == 7 + m

		%-Reset title and delete figure
		%-----------------------------------------------------------
		delete(Fdcm)
		spm('FigName',header);
		spm('Pointer','Arrow')
		spm_input('Thank you',1,'d')
		return

	%-end if OPT
	%-------------------------------------------------------------------
	end

	%-end while
	%-------------------------------------------------------------------
	end
end

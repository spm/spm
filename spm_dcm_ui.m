function [DCM] = spm_dcm_ui(Action)
% User interface for Dynamic Causal Modeling (DCM)
% FORMAT [DCM] = spm_dcm_ui('specify');
% FORMAT [DCM] = spm_dcm_ui('estimate');
% FORMAT [DCM] = spm_dcm_ui('review');
% FORMAT [DCM] = spm_dcm_ui('compare');
% FORMAT [DCM] = spm_dcm_ui('average');
%
% Specify a new model
% Estimate a specified model
% Review a previously estimated model
% Compare two or more estimated models
%
% saves DCM_???.mat
%
%	DCM.M      - model  specification structure (see spm_nlsi)
%	DCM.Y      - output specification structure (see spm_nlsi)
%	DCM.U      - input  specification structure (see spm_nlsi)
%	DCM.Ep     - posterior expectations (see spm_nlsi)
%	DCM.Cp     - posterior covariances (see spm_nlsi)
%	DCM.A      - intrinsic connection matrix
%	DCM.B      - input-dependent connection matrix
%	DCM.C      - input connection matrix
%	DCM.pA     - pA - posterior probabilities
%	DCM.pB     - pB - posterior probabilities
%	DCM.pC     - pC - posterior probabilities
%	DCM.vA     - vA - variance of parameter estimates
%	DCM.vB     - vB - variance of parameter estimates
%	DCM.vC     - vC - variance of parameter estimates
%	DCM.H1     - 1st order Volterra Kernels - hemodynamic
%	DCM.H2     - 1st order Volterra Kernels - hemodynamic
%	DCM.K1     - 1st order Volterra Kernels - neuronal
%	DCM.K1     - 1st order Volterra Kernels - neuronal
%	DCM.R      - residuals
%	DCM.y      - predicted responses
%	DCM.xY     - original response variable structures
%	DCM.T      - threshold for inference based on posterior p.d.f
%   DCM.Ce     - Estimated observation noise covariance
%   DCM.v      - Number of scans
%   DCM.n      - Number of regions
%
%___________________________________________________________________________
%
% DCM is a causal modelling procedure for dynamical systems in which
% causality is inherent in the differential equations that specify the model.
% The basic idea is to treat the system of interest, in this case the brain,
% as an input-state-output system.  By perturbing the system with known
% inputs, measured responses are used to estimate various parameters that
% govern the evolution of brain states.  Although there are no restrictions
% on the parameterisation of the model, a bilinear approximation affords a
% simple re-parameterisation in terms of effective connectivity.  This
% effective connectivity can be latent or intrinsic or, through bilinear
% terms, model input-dependent changes in effective connectivity.  Parameter
% estimation proceeds using fairly standard approaches to system
% identification that rest upon Bayesian inference.
% 
% Dynamic causal modelling represents a fundamental departure from
% conventional approaches to modelling effective connectivity in
% neuroscience.  The critical distinction between DCM and other approaches,
% such as structural equation modelling or multivariate autoregressive
% techniques is that the input is treated as known, as opposed to stochastic.
% In this sense DCM is much closer to conventional analyses of neuroimaging
% time series because the causal or explanatory variables enter as known
% fixed quantities.  The use of designed and known inputs in characterising
% neuroimaging data with the general linear model or DCM is a more natural
% way to analyse data from designed experiments.  Given that the vast
% majority of imaging neuroscience relies upon designed experiments we
% consider DCM a potentially useful complement to existing techniques.  
%___________________________________________________________________________
% Copyright (C) 2005 Wellcome Department of Imaging Neuroscience

% Karl Friston
% $Id: spm_dcm_ui.m 155 2005-05-15 11:46:28Z klaas $



% get figure handles
%---------------------------------------------------------------------------
Finter = spm_figure('GetWin','Interactive');
Fgraph = spm_figure('GetWin','Graphics');
header = get(Finter,'Name');
set(Finter,'Name','Dynamic Causal Modelling')
WS     = spm('WinScale');


% options
%---------------------------------------------------------------------------
if ~nargin
	str    = 'Action: ';
	%Action = spm_input(str,1,'b',{'specify','estimate','review','compare'});
    
    % Use pull-down menu
    Actions={'specify','estimate','review','compare','average','quit'};
    selected = spm_input(str,1,'m',Actions);
    Action=Actions{selected};
end


switch Action

% specify graph and estimate
%---------------------------------------------------------------------------
case 'specify'

	% get design and directory
	%===================================================================
	swd   = spm_str_manip(spm_select(1,'^SPM\.mat$','Select SPM.mat'),'H');
	load(fullfile(swd,'SPM.mat'))
	cd(swd)

	% name
	%===================================================================
	name  = spm_input('name for DCM_???.mat','+1','s');

	% outputs
	%===================================================================

	% get cell array of region structures
	%-------------------------------------------------------------------
	P     = spm_select([2 8],'^VOI.*\.mat$',{'select VOIs'});
	m     = size(P,1);
	for i = 1:m
		p     = load(P(i,:),'xY','-mat');
		xY(i) = p.xY;
	end

	% inputs
	%===================================================================

	% get 'causes' or inputs U
	%-------------------------------------------------------------------
	spm_input('Input specification:...  ',1,'d');
	Sess   = SPM.Sess(xY(1).Sess);
	U.dt   = Sess.U(1).dt;
	u      = length(Sess.U);
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
		pause(0.01)
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
		spm_input(str,1,'d');
        
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
                if a(i,j)==1
                    % If there is an intrinsic connection
                    % allow it to be modulated
                    h3(i,j) = uicontrol(Finter,...
                        'Position',[220+dx*j 360-dx*i 020 020].*WS,...
                        'Style','radiobutton');
                end
            end
        end
		drawnow

		% wait for 'done'
		%-----------------------------------------------------------
		set(gcf,'CurrentObject',h2(1))
		while(1)
			pause(0.01)
			if strcmp(get(gco,'Type'),'uicontrol')
			if strcmp(get(gco,'String'),'done')

			% get c
			%--------------------------------------------------
			for i = 1:m
				c(i,k)   = get(h2(i)  ,'Value');
			end

			% get b allowing any 2nd order effects 
			%--------------------------------------------------
            for i = 1:m
                for j = 1:m
                    if a(i,j)==1
                        b(i,j,k) = get(h3(i,j),'Value');
                    end
                end
            end
			delete([h1(:); h2(:); h3(find(a==1))])
			spm_input('Thank you',1,'d')
			break

			end
			end
		end
	end
	delete(d)

    % confounds (NB: the data have been filtered and whitened)
	%-------------------------------------------------------------------
	v     = size(xY(1).u,1);
	X0    = xY(1).X0;

	% response variables
	%-------------------------------------------------------------------
	n     = length(xY);
	Y.dt  = SPM.xY.RT;
	Y.X0  = X0;
	for i = 1:n

		% regional responses
		%-----------------------------------------------------------
		Y.y(:,i)  = xY(i).u;
		Y.name{i} = xY(i).name;
	end

	% error components (one for each region) - i.i.d. (because of W)
	%-------------------------------------------------------------------
	Y.Ce  = spm_Ce(ones(1,n)*v);

    DCM.a=a;
    DCM.b=b;
    DCM.c=c;
    DCM.U=U;
    DCM.Y=Y;
    DCM.xY=xY;
    DCM.v=v;
    DCM.n=n;
    
    
    %-Save and reset title
	%-------------------------------------------------------------------
	if str2num(version('-release'))>=14,
		save(fullfile(swd,['DCM_' name]),'-V6','DCM');
	else
		save(fullfile(swd,['DCM_' name]),'DCM');
	end;
	spm('FigName',header);
	spm('Pointer','Arrow')
	return
    
case 'estimate',
	
    
    spm_dcm_estimate;

	return


% review results
%---------------------------------------------------------------------------
case 'review'

	%-display model details
	%-------------------------------------------------------------------
	set(Finter,'name','Dynamic Causal Modeling')

	%-get results
	%-------------------------------------------------------------------
	P     = spm_select(1,'^DCM.*\.mat$','select DCM_???.mat');
	load(P);
	m     = DCM.M.m;
	n     = DCM.M.n;
	l     = DCM.M.l;

	%-get threshold and recompute posterior probabilities
	%-------------------------------------------------------------------
	str        = 'Threshold {Hz}';
	DCM.T      = spm_input(str,1,'e',0,[1 1]);
	pp         = 1 - spm_Ncdf(DCM.T,abs(DCM.Ep),diag(DCM.Cp));
	[pA pB pC] = spm_dcm_reshape(pp,m,l,1);
	DCM.pA     = pA;
	DCM.pB     = pB;
	DCM.pC     = pC;

	while(1)
	% get options
	%-------------------------------------------------------------------
	str   = {};
	for i = 1:m
		str{i} = ['    Effects of ' DCM.U.name{i} '    '];
	end

	str   = {str{:}		'    Intrinsic connections',...
				'Contrast of connections',...
				'location of regions',...
				'Inputs',...
				'Outputs',...
				'1st order Kernels',...
				'quit'};

	OPT   = spm_input('display',2,'m',str);
	figure(Fgraph)
	spm_figure('Clear',Fgraph)


	% Inputs
	%-------------------------------------------------------------------
	if OPT <= m 

		% input effects
		%-----------------------------------------------------------
		subplot(2,1,1)
		i        = OPT;
		b        = DCM.pB(:,:,i);
		b(:,:,2) = DCM.B(:,:,i);
		c        = DCM.pC(:,i);
		c(:,2)   = DCM.C(:,i);
		spm_dcm_display(DCM.xY,b,c)
		title({	str{OPT};...
			sprintf('P(|connection| > %0.2f)',DCM.T);...
			'connection strength'},'FontSize',12)


		% direct effects - connections
		%-----------------------------------------------------------
		subplot(4,2,5)
		bar(DCM.C(:,i))
		title('C - direct effects {Hz}','FontSize',12)
		set(gca,'XTick',[1:l],'XTickLabel',DCM.Y.name)
		axis square
		grid on

		% direct effects - probabilities
		%-----------------------------------------------------------
		subplot(4,2,6)
		bar(DCM.pC(:,i))
		title('C {probabilities}','FontSize',12)
		set(gca,'XTick',[1:l],'XTickLabel',DCM.Y.name)
		ylabel(sprintf('P(A > %0.2f)',DCM.T))
		axis square
		grid on

		% modulatory effects - connections
		%-----------------------------------------------------------
		subplot(4,2,7)
		bar3(DCM.B(:,:,i),0.4)
		set(gca,'YTick',[1:l],'YTickLabel',DCM.Y.name)
		set(gca,'XTick',[1:l],'XTickLabel',DCM.Y.name)
		title('B - modulatory effects {Hz}','FontSize',12)
		xlabel('from')
		ylabel('to')
		axis square

		% modulatory effects - probabilities
		%-----------------------------------------------------------
		subplot(4,2,8)
		bar3(DCM.pB(:,:,i),0.4)
		title('B {probabilities}','FontSize',12)
		set(gca,'XTick',[1:l],'XTickLabel',DCM.Y.name)
		set(gca,'YTick',[1:l],'YTickLabel',DCM.Y.name)
		xlabel('from')
		ylabel('to')
		zlabel(sprintf('P(A > %0.2f)',DCM.T))
		axis square
		grid on


	% Intrinsic connections
	%-------------------------------------------------------------------
	elseif OPT == 1 + m

		% Intrinsic effects
		%-----------------------------------------------------------
		subplot(2,1,1)
		a        = DCM.pA;
		a(:,:,2) = DCM.A;
		spm_dcm_display(DCM.xY,a)
		title({	str{OPT};...
			sprintf('P(|connection| > %0.2f)',DCM.T);...
			'connection strength'},'FontSize',10)

		% intrinsic interactions
		%-----------------------------------------------------------
		subplot(2,2,3)
		bar3(DCM.A,0.4)
		title('A - Intrinsic connections','FontSize',12)
		set(gca,'XTick',[1:l],'XTickLabel',DCM.Y.name)
		set(gca,'YTick',[1:l],'YTickLabel',DCM.Y.name)
		xlabel('from')
		ylabel('to')
		grid on
		axis square

		% intrinsic interactions - probabilities
		%-----------------------------------------------------------
		subplot(2,2,4)
		bar3(DCM.pA,0.4)
		title('A {probabilities}','FontSize',12)
		set(gca,'XTick',[1:l],'XTickLabel',DCM.Y.name)
		set(gca,'YTick',[1:l],'YTickLabel',DCM.Y.name)
		xlabel('from')
		ylabel('to')
		zlabel(sprintf('P(A > %0.2f)',DCM.T))
		grid on
		axis square


	% contrast
	%-------------------------------------------------------------------
	elseif OPT == 2 + m

		%-get contrast
		%-----------------------------------------------------------
		str     = 'contrast for';
		D       = spm_input(str,1,'b',{'A','B','C'});

		switch D

			case 'A' % intrinsic connections
			%---------------------------------------------------
            C       = spm_dcm_contrasts(P(:),'A');
			i       = find(C); j = 1;
			C       = sparse(i + j,1,C(i),length(DCM.Ep),1);
            
			case 'B' % modulatory inputs
			%---------------------------------------------------
            C       = spm_dcm_contrasts(P(:),'B');
			i       = find(C); j = 1 + l*l;
			C       = sparse(i + j,1,C(i),length(DCM.Ep),1);

			case 'C' % direct (driving) inputs
			%---------------------------------------------------
            C       = spm_dcm_contrasts(P(:),'C');
			i       = find(C); j = 1 + l*l + l*l*m;
			C       = sparse(i + j,1,C(i),length(DCM.Ep),1);

		end

		%-posterior density and inference
		%-----------------------------------------------------------
		c    = C'*DCM.Ep;
		v    = C'*DCM.Cp*C;
		x    = c + [-32:32]*sqrt(v)*6/32;
		p    = full(1/sqrt(2*pi*v)*exp(-[x - c].^2/(2*v)));  % conversion to full necessary to account for sparse matrices bug in MATLAB 6.5.0 R13 
		PP   = 1 - spm_Ncdf(DCM.T,c,v);

		figure(Fgraph)
		subplot(2,1,1)
		plot(x,p,[1 1]*DCM.T,[0 max(p)],'-.');
		title({'Posterior density of contrast',...
		sprintf('P(contrast > %0.2f) = %.1f%s',DCM.T,PP*100,'%')},...
			'FontSize',12)
		xlabel('contrast')
		ylabel('probability density')

		i    = find(x >= DCM.T);
		hold on
		fill([x(i) fliplr(x(i))],[i*0 fliplr(p(i))],[1 1 1]*.8)
		axis square, grid on
		hold off

		%-contrast
		%-----------------------------------------------------------
		[A B C] = spm_dcm_reshape(C,m,l,1);
		subplot(4,2,5)
		imagesc(A)
		title('contrast {A}','FontSize',12)
		set(gca,'YTick',[1:l],'YTickLabel',DCM.Y.name)
		set(gca,'XTick',[1:l],'XTickLabel',DCM.Y.name)
		axis image

		subplot(4,2,6)
		imagesc(C)
		title('contrast {C}','FontSize',12)
		set(gca,'XTick',[1:m],'XTickLabel',DCM.U.name)
		set(gca,'YTick',[1:l],'YTickLabel',DCM.Y.name)
		axis image

		for i = 1:m
			subplot(4,m,i + 3*m)
			imagesc(B(:,:,i))
			title(['contrast {B}-' DCM.U.name{i}],'FontSize',12)
			set(gca,'XTick',[1:l],'XTickLabel',DCM.Y.name)
			set(gca,'YTick',[1:l],'YTickLabel',DCM.Y.name)
			axis image
		end



 	% location of regions
	%-------------------------------------------------------------------
	elseif OPT == 3 + m


		% transverse
		%-----------------------------------------------------------
		subplot(2,2,3)
		title('transverse')
		u = [[0 1 0];[-1 0 0]]';
		spm_dcm_display(DCM.xY,[],[],u,32);

		% sagittal
		%-----------------------------------------------------------
		subplot(2,2,1)
		title('sagittal')
		u = [[0 1 0];[0 0 1]]';
		spm_dcm_display(DCM.xY,[],[],u,32);

		% coronal
		%-----------------------------------------------------------
		subplot(2,2,2)
		title('coronal')
		u = [[1 0 0];[0 0 1]]';
		spm_dcm_display(DCM.xY,[],[],u,32);

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


 	% inputs
	%-------------------------------------------------------------------
	elseif OPT == 4 + m

		% graph
		%-----------------------------------------------------------
		x     = [1:length(full(DCM.U.u))]*DCM.U.dt;
		for i = 1:m
			subplot(m,1,i)
			plot(x,full(DCM.U.u(:,i)))
			title(['Input - ' DCM.U.name{i}],'FontSize',12)
			%axis tight
			ylabel('event density {Hz}')
			if i==m	
				xlabel('time {seconds}')
			end
			%axis square tight
			grid on
            % Change axis so we can see boxcars properly
            u_min=min(full(DCM.U.u(:,i)));
            u_max=max(full(DCM.U.u(:,i)));
            u_range=u_max-u_min;
            axis([0 max(x) u_min-0.1*u_range u_max*1.1]);
		end


	% outputs
	%-------------------------------------------------------------------
	elseif OPT == 5 + m

		% graph
		%-----------------------------------------------------------
		x  = [1:length(DCM.y)]*DCM.Y.dt;;
		for i = 1:l
            subplot(l,1,i);
            plot(x,DCM.Y.y(:,i));
            if i==l
		xlabel('time {seconds}');
	    end
 			grid on
            hold on
            X0=DCM.Y.X0;
            
            % Remove zero columns from X0 if there are any 
            % (although there should'nt be !)
            % Otherwise we won't be able to evaluate inv(X0'X0)
            ncol_X0=size(X0,2);
            new_X0=[];
            for col_X0=1:ncol_X0,
                if ~(length(find(X0(:,col_X0)==0))==DCM.v)
                    new_X0=[new_X0 X0(:,col_X0)];
                end
            end
            X0=new_X0;
         
            % Note: adding on X0*beta to DCM prediction y
            % is equivalent to orthogonalising y WRT X0
            % because data have already been orth wrt X0
            % (overall_prediction=DCM.y(:,i)+X0*beta;)
            overall_prediction=DCM.y(:,i)-X0*inv(X0'*X0)*X0'*DCM.y(:,i);
            plot(x,overall_prediction,'r');
            title([DCM.Y.name{i} ': data and model predictions' ],'FontSize',12);
            axis normal
        end
        disp(' ');
        
        
	% Kernels
	%-------------------------------------------------------------------
	elseif OPT == 6 + m

		% input effects
		%-----------------------------------------------------------
		x     = [1:DCM.M.N]*DCM.M.dt;
		d     = 2/DCM.M.dt;
		for i = 1:m

			% input effects - neuronal
			%---------------------------------------------------
			y     = DCM.K1(:,:,i);
			subplot(m,2,2*(i - 1) + 1)
			plot(x,y)
			set(gca,'XLim',[0 16])
			axis square
			title(['neuronal responses to ' DCM.U.name{i}])
			grid on
			xlabel('time {seconds}')
			for j = 1:l
				text(x(j*d),y(j*d,j),DCM.Y.name{j},...
					'FontSize',8,...
					'HorizontalAlignment','Center')
			end

			% input effects - hemodynamic
			%---------------------------------------------------
			y     = DCM.K1(:,:,i);
			k     = DCM.H1(:,:,i);
			subplot(m,2,2*(i - 1) + 2)
			plot(x,k,x,y,':')
			set(gca,'XLim',[0 16])
			axis square
			title('hemodynamic responses')
			grid on
			xlabel('time {seconds}')
			for j = 1:l
				text(x(j*d),k(j*d,j),DCM.Y.name{j},...
					'FontSize',8,...
					'HorizontalAlignment','Center')
			end
		end

	% quit
	%-------------------------------------------------------------------
	elseif OPT == 7 + m

		%-Reset title and delete figure
		%-----------------------------------------------------------
		spm_clf
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
    

% compare different models
%---------------------------------------------------------------------------
case 'compare',
    
    num_models = spm_input('Number of DCM models to compare','+1','r',[],1);
	P     = spm_select(num_models,'^DCM.*\.mat$','select DCM*.mat files');
    
    % load all models and compute their evidence
    for model_index=1:num_models,
	    load(P(model_index,:),'DCM','-mat');
        model_VOIs{model_index}   = [];
        % concatenate names of all VOIs for each model
        for k   = 1:size(DCM.xY,2),
            model_VOIs{model_index}    = [model_VOIs{model_index} DCM.xY(k).name];
        end
        evidence(model_index)       = spm_dcm_evidence(DCM);
        aic(model_index)            = evidence(model_index).aic_overall;
        bic(model_index)            = evidence(model_index).bic_overall;
    end
    
    % model comparison is only valid if y is identical:
    % check that all models refer to the same set of VOIs
    for model_index=1:num_models,
        if ~strcmp(model_VOIs{1}, model_VOIs{model_index})
        	str = {	'Selected models contain different sets of VOIs!',...
            		'Model comparison only valid for models with identical VOIs.',...
                    'Procedure aborted.'};
	        spm_input(str,1,'bd','OK',[1],1);
            return
        end
    end
    
    % Get and plot posterior probabilities of models assuming
    % each model has equal probability a-priori 
    maic=aic-mean(aic);
    pm=exp(maic)/sum(exp(maic));

    figure(Fgraph);
    spm_clf
    subplot(2,1,1);
    bar(pm);
    Vscale=axis;
    Vscale(4)=1;
    axis(Vscale);
    set(gca,'FontSize',18);
    ylabel('p(y|m)');  % given the flat prior p(m), p(y|m)=p(m|y)
    xlabel('m');
    title('Posterior probabilities of models from AIC');
    
    mbic=bic-mean(bic);
    pm=exp(mbic)/sum(exp(mbic));

    subplot(2,1,2);
    bar(pm);
    Vscale=axis;
    Vscale(4)=1;
    axis(Vscale);
    set(gca,'FontSize',18);
    ylabel('p(y|m)');  % given the flat prior p(m), p(y|m)=p(m|y)
    xlabel('m');
    title('Posterior probabilities of models from BIC');
    
    % Output model comparison details to MATLAB command window
    for ii=1:num_models,
        for jj=1:num_models,
            if ~(jj==ii)
                disp('---------------------------------------------------------------');
                    disp(sprintf('Model %d: %s',ii,P(ii)));
                    disp('          versus ');
                    disp(sprintf('Model %d: %s',jj,P(jj)));
                    diff=1;
                disp(' ');
                disp('All costs are in units of binary bits');
                disp(' ');
                for k=1:size(DCM.A,1),
                    nats=-(evidence(ii).region_cost(k)-evidence(jj).region_cost(k));
                    bits=nats/log(2);
                    disp(sprintf('Region %s: relative cost  = %1.2f, BF= %1.2f',DCM.Y.name{k},bits,2^(-bits)));
                end
                % AIC penalty
                nats=evidence(ii).aic_penalty-evidence(jj).aic_penalty;
                bits=nats/log(2);
                disp(sprintf('AIC Penalty = %1.2f, BF = %1.2f',bits,2^(-bits)));
                % BIC penalty
                nats=evidence(ii).bic_penalty-evidence(jj).bic_penalty;
                bits=nats/log(2);
                disp(sprintf('BIC Penalty = %1.2f, BF = %1.2f',bits,2^(-bits)));
                % AIC overall
                nats=-diff*(aic(ii)-aic(jj));
                bits=nats/log(2);
                bf_aic=2^(-bits);
                disp(sprintf('AIC Overall = %1.2f, BF = %1.2f',bits,bf_aic));
                % BIC overall
                nats=-diff*(bic(ii)-bic(jj));
                bits=nats/log(2);
                bf_bic=2^(-bits);
                disp(sprintf('BIC Overall = %1.2f, BF = %1.2f',bits,bf_bic));
                disp(' ');
                
                if (bf_bic > exp(1)) & (bf_aic > exp(1))
                    disp(sprintf('Consistent evidence in favour of model %d',ii));
                    disp(sprintf('Bayes factor >= %1.2f', min(bf_aic,bf_bic)));
                    disp(' ');
                elseif ((1/bf_bic) > exp(1)) & ((1/bf_aic) > exp(1))
                    disp(sprintf('Consistent evidence in favour of model %d',jj));
                     disp(sprintf('Bayes factor >= %1.2f', min(1/bf_aic,1/bf_bic)));
                    disp(' ');
                else
                    disp('No consistent evidence in favour of either model');
                    disp(' ');
                end
                disp('---------------------------------------------------------------');
                
            end
        end
    end
    
    
    spm_input('Model comparison details in MATLAB command window',1,'d');
    spm_input('Thank you',3,'d');
    return
    
case 'average',
    
    spm_dcm_average;
    spm_clf
    spm('FigName',header);
    spm('Pointer','Arrow')
    spm_input('Thank you',1,'d')
    return
    
case 'quit',
    spm_clf
    spm('FigName',header);
    spm('Pointer','Arrow')
    spm_input('Thank you',1,'d')
    return
end

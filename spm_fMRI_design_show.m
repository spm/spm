function [xX,Sess] = spm_fMRI_design_show(xX,Sess,s,i)
% Interactive review of fMRI design matrix
% FORMAT [xX,Sess] = spm_fMRI_design_show(xX,Sess)
%
% xX            - structure describing design matrix
% xX.X          - design matrix
% xX.dt         - time bin {secs}
% xX.RT         - Repetition time {secs}
% xX.iH         - vector of H partition (condition effects)      indices
% xX.iC         - vector of C partition (covariates of interest) indices
% xX.iB         - vector of B partition (block effects)          indices
% xX.iG         - vector of G partition (nuisance variables)     indices
% xX.Xnames     - cellstr of effect names corresponding to columns
%                 of the design matrix
%
% Sess{s}.BFstr    - basis function description string
% Sess{s}.DSstr    - Design description string
% Sess{s}.rep      - session replication flag
% Sess{s}.row      - scan   indices      for session s
% Sess{s}.col      - effect indices      for session s
% Sess{s}.name{i}  - of ith trial type   for session s
% Sess{s}.ind{i}   - column indices      for ith trial type {within session}
% Sess{s}.bf{i}    - basis functions     for ith trial type
% Sess{s}.sf{i}    - stick functions     for ith trial type
% Sess{s}.ons{i}   - stimuli onset times for ith trial type (secs)
% Sess{s}.pst{i}   - peristimulus times  for ith trial type (secs)
% Sess{s}.Pv{i}    - vector of paramters for ith trial type
% Sess{s}.Pname{i} - name   of paramters for ith trial type
%_______________________________________________________________________
% %W% Karl Friston %E%


% xX and Sess
%-----------------------------------------------------------------------
if nargin == 0
	load(spm_get(1,'fMRIDesMtx.mat','select SPM_fMRIDesMtx'))
elseif nargin<2
	error('insufficient arguments')
end


% Do not proceed unless there are trials specified
%-----------------------------------------------------------------------
for j = 1:length(Sess)
    if ~length(Sess{j}.name)
        spm('alert*','User-specifed regressors only!',mfilename,sqrt(-1));
        return
    end
end


%-Defaults: Setup GUI 
%-----------------------------------------------------------------------
if nargin < 3
	s = 1;
	i = 1;

	%-Get Interactive window and delete any previous DesRepUI menu
	%---------------------------------------------------------------
	Finter = spm_figure('GetWin','Interactive');
	delete(findobj(get(Finter,'Children'),'flat','Tag','DesRepUI'))

	%-Add a scaled design matrix to the design data structure
	%---------------------------------------------------------------
	if ~isfield(xX,'nKX'), xX.nKX = spm_DesMtx('Sca',xX.X,xX.Xnames); end

	%-Draw menu
	%---------------------------------------------------------------
	hC     = uimenu(Finter,'Label','Explore fMRI design',...
		'Separator','on',...
		'Tag','DesRepUI',...
		'UserData',struct('xX',xX,'Sess',{Sess}),...
		'HandleVisibility','on');
	for j = 1:length(Sess)
		h     = uimenu(hC,'Label',sprintf('Session %.0f ',j),...
			'HandleVisibility','off');
		for k = 1:length(Sess{j}.name)
			cb = ['tmp = get(get(gcbo,''UserData''),',...
					         '''UserData''); ',...
				sprintf(['spm_fMRI_design_show(',...
					'tmp.xX,tmp.Sess,%d,%d);'],j,k)];
			uimenu(h,'Label',Sess{j}.name{k},...
	     	   	         'CallBack',cb,...
	     	   	         'UserData',hC,...
	     	   	         'HandleVisibility','off')
		end
	end
end


%-Graphics...
%=======================================================================

%-Get Graphics window
%-----------------------------------------------------------------------
Fgraph = spm_figure('GetWin','Graphics');
spm_results_ui('Clear',Fgraph,0)


% Display design matrix X
%-----------------------------------------------------------------------
axes('Position',[0.125,0.700,0.155,0.225])
if isfield(xX,'nKX')
	hDesMtxIm = image(xX.nKX*32+32);
else
	hDesMtxIm = imagesc(spm_en(xX.X));
end
xlabel('effect')
ylabel('scan')
title('Design Matrix','FontSize',16)

%-Setup callbacks to allow interrogation of design matrix
%-----------------------------------------------------------------------
set(hDesMtxIm,'UserData',struct('X',xX.X,'Xnames',{xX.Xnames}))
set(hDesMtxIm,'ButtonDownFcn','spm_DesRep(''SurfDesMtx_CB'')')



% Session subpartition
%-----------------------------------------------------------------------
axes('Position',[0.550,0.700,0.155,0.225])
sX   = xX.X(Sess{s}.row,Sess{s}.col);
imagesc(spm_en(sX)')
set(gca,'YTick',[1:size(sX,1)])
set(gca,'YTickLabel',xX.Xnames(Sess{s}.col)')
title({sprintf('Session %d',s) Sess{s}.DSstr})

% Collinearity
%-----------------------------------------------------------------------
axes('Position',[0.750,0.700,0.155,0.225])
imagesc(corrcoef(sX))
title('correlations')
axis off, axis square

% Trial-specific regressors - time domain
%-----------------------------------------------------------------------
rX    = sX(:,Sess{s}.ind{i});
axes('Position',[0.125,0.405,0.325,0.225])
plot(Sess{s}.row,rX)
xlabel('scan')
ylabel('regressor[s]')
title({['Regressors for ' Sess{s}.name{i}] })
axis tight

% Trial-specific regressors - frequency domain
%-----------------------------------------------------------------------
axes('Position',[0.580,0.405,0.325,0.225])
gX    = abs(fft(rX)).^2;
gX    = gX*diag(1./sum(gX));
q     = size(gX,1);
Hz    = [0:(q - 1)]/(q*xX.RT);
q     = 2:fix(q/2);
plot(Hz(q),gX(q,:))
xlabel('Frequency (Hz)')
ylabel('spectral density')
title('Frequency domain')
grid on
axis tight


% if trial (as opposed to trial x trial interaction)
%-----------------------------------------------------------------------
if length(Sess{s}.ons) >= i

	% Basis set and peristimulus sampling
	%---------------------------------------------------------------
	axes('Position',[0.125,0.110,0.325,0.225])
	t    = [1:size(Sess{s}.bf{i},1)]*xX.dt;
	pst  = Sess{s}.pst{i};
	plot(t,Sess{s}.bf{i},pst,0*pst,'.','MarkerSize',16)
	str  = sprintf('TR = %0.0fsecs',xX.RT);
	xlabel({'time (secs)' str sprintf('%0.0fms time bins',1000*xX.dt)})
	title({'Basis set and peristimulus sampling' Sess{s}.BFstr})
	axis tight
	grid on

	% if a paramteric variate is specified
	%---------------------------------------------------------------
	if length(Sess{s}.Pv{i})

		% onsets and parametric modulation
		%-------------------------------------------------------
		axes('Position',[0.580,0.110,0.325,0.225])
		plot(Sess{s}.ons{i},Sess{s}.Pv{i},'.','MarkerSize',8)
		title({'trial specific parameters' Sess{s}.Pname{i}})
		xlabel('time (secs}')
		ylabel(Sess{s}.Pname{i})
		grid on
	end
end

%-Pop up Graphics figure window
%-----------------------------------------------------------------------
figure(Fgraph);

function [xX,Sess] = spm_fMRI_design_show(xX,Sess,s,i)
% Interactive review of fMRI design matrix
% FORMAT [xX,Sess] = spm_fMRI_design_show(xX,Sess,s,i)
%
% xX            - structure describing design matrix
% xX.X          - design matrix
% xX.dt         - time bin {secs}
% xX.RT         - Repetition time {secs}
% xX.iH         - vector of H partition (condition effects)      indices,
% xX.iC         - vector of C partition (covariates of interest) indices
% xX.iB         - vector of B partition (block effects)          indices
% xX.iG         - vector of G partition (nuisance variables)     indices
% xX.Xnames     - cellstr of effect names corresponding to columns
%                 of the design matrix
%
% Sess{s}.U{i}  -  see spm_fMRI_design for session s, trial i.
%
%_______________________________________________________________________
% %W% Karl Friston %E%

% xX and Sess
%-----------------------------------------------------------------------
if nargin == 0
	load(spm_get(1,'fMRIDesMtx.mat','Select SPM_fMRIDesMtx.mat'));
elseif nargin < 2
	error('insufficient arguments')
end


% Do not proceed unless there are trials specified
%-----------------------------------------------------------------------
for j = 1:length(Sess)
    if ~length(Sess{j}.U)
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

	%-Create menu
	%---------------------------------------------------------------
	hC     = uimenu(Finter,'Label','Explore fMRI design',...
		'Separator','on',...
		'Tag','DesRepUI',...
		'UserData',struct('xX',xX,'Sess',{Sess}),...
		'HandleVisibility','on');
	for j = 1:length(Sess)
		h     = uimenu(hC,'Label',sprintf('Session %.0f ',j),...
			'HandleVisibility','off');
		for k = 1:length(Sess{j}.Fcname)
			cb = ['tmp = get(get(gcbo,''UserData''),',...
					         '''UserData''); ',...
				sprintf(['spm_fMRI_design_show(',...
					'tmp.xX,tmp.Sess,%d,%d);'],j,k)];
			uimenu(h,'Label',Sess{j}.Fcname{k},...
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


% Trial-specific regressors - time domain
%-----------------------------------------------------------------------
sX    = xX.X(Sess{s}.row,Sess{s}.col);
rX    = sX(:,Sess{s}.Fci{i});
subplot(2,2,1)
plot(Sess{s}.row,rX)
xlabel('scan')
ylabel('regressor[s]')
title({'Time domain',['regressors for ' Sess{s}.Fcname{i}]})
grid on
axis tight

% Trial-specific regressors - frequency domain
%-----------------------------------------------------------------------
subplot(2,2,2)
gX    = abs(fft(rX)).^2;
gX    = gX*diag(1./sum(gX));
q     = size(gX,1);
Hz    = [0:(q - 1)]/(q*xX.RT);
q     = 2:fix(q/2);
plot(Hz(q),gX(q,:))
patch([0 1 1 0]/128,[0 0 1 1]*max(max(gX)),[1 1 1]*.9)
xlabel('Frequency (Hz)')
ylabel('relative spectral density')
title({'Frequency domain','128 second High-pass filter'})
grid on
axis tight


% if trial (as opposed to trial x trial interaction)
%-----------------------------------------------------------------------
if length(Sess{s}.U) >= i

	% Basis set and peristimulus sampling
	%---------------------------------------------------------------
	subplot(2,2,3)
	dt   = Sess{s}.U{i}.dt;
	t    = [1:size(Sess{s}.bf,1)]*dt;
	pst  = Sess{s}.U{i}.pst;
	plot(t,Sess{s}.bf,pst,0*pst,'.','MarkerSize',16)
	str  = sprintf('TR = %0.2fsecs',xX.RT);
	xlabel({'time (secs)' str sprintf('%0.0fms time bins',1000*dt)})
	title({'Basis set and peristimulus sampling' Sess{s}.Bfname})
	axis tight
	grid on

	% if a paramteric variate is specified
	%---------------------------------------------------------------
	if length(Sess{s}.U{i}.P)

		% onsets and parametric modulation
		%-------------------------------------------------------
		subplot(2,2,4)
		plot(Sess{s}.U{i}.ons,Sess{s}.U{i}.P,'.','MarkerSize',8)
		title({'trial specific parameters' Sess{s}.U{i}.Pname{:}})
		xlabel('time {secs}')
		ylabel('parameter value')
		grid on
	end
end

%-Pop up Graphics figure window
%-----------------------------------------------------------------------
figure(Fgraph);

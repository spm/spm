function [X,Sess] = spm_fMRI_design_show(X,Sess,s,i)
% Interactive review of fMRI design matrix
% FORMAT [X,Sess] = spm_fMRI_design_show(X,Sess)
%
% X.dt    - time bin {secs}
% X.RT    - Repetition time {secs}
% X.xX    - regressors
% X.bX    - session effects
% X.Xname - names of subpartiton columns {1xn}
% X.Bname - names of subpartiton columns {1xn}
%
% Sess{s}.BFstr   - basis function description string
% Sess{s}.DSstr   - Design description string
% Sess{s}.row     - scan   indices      for session s
% Sess{s}.col     - effect indices      for session s
% Sess{s}.name{i} - of ith trial type   for session s
% Sess{s}.ind{i}  - column indices      for ith trial type {within session}
% Sess{s}.bf{i}   - basis functions     for ith trial type
% Sess{s}.ons{i}  - stimuli onset times for ith trial type (secs)
% Sess{s}.pst{i}  - peristimulus times  for ith trial type (secs)
% Sess{s}.para{i} - vector of paramters for ith trial type
%_______________________________________________________________________
% %W% Karl Friston %E%


% X and Sess
%-----------------------------------------------------------------------
if nargin == 0
	load(spm_get(1,'.mat','select fMRIDesMtx'))
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

	%-Get Interactive window and delete any previous fMRIDesShow menu
	%---------------------------------------------------------------
	Finter = spm_figure('GetWin','Interactive');
	delete(findobj(get(Finter,'Children'),'flat','Tag','fMRIDesShow'))

	%-Add a scaled design matrix to the design data structure
	%---------------------------------------------------------------
	if ~isfield(X,'nxbX'), X.nxbX = spm_DesMtx('Sca',[X.xX X.bX]); end

	%-Draw menu
	%---------------------------------------------------------------
	hC     = uimenu(Finter,'Label','explore fMRI design',...
		'Separator','on',...
		'Tag','fMRIDesShow',...
		'UserData',struct('X',X,'Sess',{Sess}),...
		'HandleVisibility','on');
	for j = 1:length(Sess)
		h = uimenu(hC,'Label',sprintf('Session %.0f ',j),...
			'HandleVisibility','off');
		for k = 1:length(Sess{j}.name)
			cb = ['tmp = get(get(gcbo,''UserData''),',...
					'''UserData''); ',...
				sprintf(['spm_fMRI_design_show(',...
					'tmp.X,tmp.Sess,%d,%d);'],j,k)];
			uimenu(h,'Label',Sess{j}.name{k},...
	     	   	         'CallBack',cb,...
	     	   	         'UserData',hC,...
	     	   	         'HandleVisibility','off')
		end
	end

% 	uimenu(hC,'Label','Quit','Separator','on',...
% 		'CallBack','delete(get(gcbo,''UserData''))',...
% 		'UserData',hC,...
% 		'HandleVisibility','off')

end


%-Graphics...
%=======================================================================

%-Get Graphics window
%-----------------------------------------------------------------------
Fgraph = spm_figure('GetWin','Graphics');
spm_results_ui('Clear',Fgraph,0)


% Display X
%-----------------------------------------------------------------------
subplot(3,4,1)
if isfield(X,'nxbX')
	image(X.nxbX*32+32)
else
	imagesc(spm_en([X.xX X.bX]))
end
xlabel('effect')
ylabel('scan')
title('Design Matrix','FontSize',16)


% Session subpartition
%-----------------------------------------------------------------------
subplot(3,4,3)
sX   = X.xX(Sess{s}.row,Sess{s}.col);
imagesc(spm_en(sX)')
set(gca,'YTick',[1:size(sX,1)])
set(gca,'YTickLabel',X.Xname(Sess{s}.col)')
title(sprintf('Session %d',s),'FontSize',16)

% Collinearity
%-----------------------------------------------------------------------
subplot(3,4,4)
imagesc(corrcoef(sX))
title('correlations')
axis off, axis square

% Trial-specific regressors
%-----------------------------------------------------------------------
subplot(3,2,3)
plot(Sess{s}.row,sX(:,Sess{s}.ind{i}))
xlabel('scan')
ylabel('regressor[s]')
title(['Regressors for ' Sess{s}.name{i}])
axis tight


% Basis set and peristimulus sampling
%-----------------------------------------------------------------------
subplot(3,2,4)
t    = [1:size(Sess{s}.bf{i},1)]*X.dt;
pst  = Sess{s}.pst{i};
plot(t,Sess{s}.bf{i},pst,0*pst,'.','MarkerSize',16)
xlabel('time (secs}')
title('Basis set and peristimulus sampling')
axis tight
grid on

% onsets and parametric modulation
%-----------------------------------------------------------------------
subplot(3,2,6)
plot(Sess{s}.ons{i},Sess{s}.para{i},'.','MarkerSize',8)
xlabel('time (secs)')
title('trial specific onsets or causes')
grid on

% textual data
%-----------------------------------------------------------------------
subplot(3,2,5)
cla
q     = 0.1;
text(0,1 - 1*q,'Design','FontWeight','bold')
text(0,1 - 2*q,X.DSstr ,'FontSize',10)
text(0,1 - 3*q,'Basis functions','FontWeight','bold')
text(0,1 - 4*q,X.BFstr,'FontSize',10)
text(0,1 - 6*q,sprintf('%0.0f ms time bins',1000*X.dt),'FontSize',10)
text(0,1 - 7*q,sprintf('%0.0f of %0.0f session[s],',s,length(Sess)),'FontSize',10)
text(0,1 - 8*q,sprintf('%0.0f of %0.0f trial[s],',i,length(Sess{s}.name)),'FontSize',10)
text(0,1 - 10*q,Sess{s}.name{i},'FontWeight','bold')
axis off


%-Pop up Graphics figure window
%-----------------------------------------------------------------------
figure(Fgraph);

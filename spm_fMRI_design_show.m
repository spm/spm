function spm_fMRI_design_show(X,Sess,s,i)
% Interactive review of fMRI design matrix
% FORMAT spm_fMRI_design_show(X,Sess)
%
% X.dt    - time bin {secs}
% X.DesN  - X.DesN{1} study type, X.DesN{2} basis function description
% X.xX    - regressors
% X.bX    - session effects
% X.Xname - names of subpartiton columns {1xn}
% X.Bname - names of subpartiton columns {1xn}
%
% Sess{s}.row     - scan indices        for session s
% Sess{s}.col     - Des. matrix offset  for session s
% Sess{s}.name{i} - of ith trial type   for session s
% Sess{s}.ind{i}  - Des. matrix indices for ith trial type
% Sess{s}.bf{i}   - basis functions     for ith trial type
% Sess{s}.ons{i}  - stimuli onset times for ith trial type (secs)
% Sess{s}.pst{i}  - peristimulus times  for ith trial type (secs)
% Sess{s}.para{i} - vector of paramters for ith trial type
%___________________________________________________________________________
% %W% Karl Friston %E%

% display Design matrix {X}
%===========================================================================


% Initialize variables and menu
%---------------------------------------------------------------------------
Finter = spm_figure('FindWin','Interactive');
Fgraph = spm_figure('FindWin','Graphics');
figure(Fgraph);


% X and Sess
%---------------------------------------------------------------------------
if nargin == 0
	load(spm_get(1,'.mat','select fMRIDesMtx'))
	set(gcf,'UserData',{X Sess})
end

% defaults
%---------------------------------------------------------------------------
if nargin < 4
	s = 1;
	i = 1;


	hC    = uimenu(Fgraph,'Label','Display',...
		'HandleVisibility','CallBack','Separator','on');
	for j = 1:length(Sess)
		h = uimenu(hC,'Label',sprintf('Session %.0f ',j),...
	            'HandleVisibility','CallBack');
		for k = 1:length(Sess{j}.name)
			str = sprintf('q = get(gcf,''UserData'');spm_fMRI_design_show(q{1},q{2},%d,%d)',j,k);
			uimenu(h,'Label',Sess{j}.name{k},...
	     	   	'HandleVisibility','CallBack','CallBack',str)
		end
	end

	% Quit
	%-------------------------------------------------------------------
	str = 'spm_clf; h = get(gcf,''Children''); delete(h(length(h)));';
	uimenu(hC,'Label','Quit',...
	'HandleVisibility','CallBack','CallBack',str,'Separator','on');

	% Display X
	%-------------------------------------------------------------------
	subplot(3,4,1)
	imagesc(spm_en([X.xX X.bX]))
	xlabel('effect')
	ylabel('scan')
	title('Design Matrix','FontSize',16)

end

% Session subpartition
%---------------------------------------------------------------------------
subplot(3,4,3)
sX   = X.xX(Sess{s}.row,Sess{s}.col);
imagesc(spm_en(sX)')
set(gca,'YTick',[1:size(sX,1)])
set(gca,'YTickLabel',X.Xname(Sess{s}.col)')
title(sprintf('Session %d',s),'FontSize',16)

% Collinearity
%---------------------------------------------------------------------------
subplot(3,4,4)
imagesc(corrcoef(sX))
title('correlations')
axis off, axis square

% Trial-specific regressors
%---------------------------------------------------------------------------
subplot(3,2,3)
plot(Sess{s}.row,sX(:,Sess{s}.ind{i}))
xlabel('scan')
ylabel('regressor[s]')
title(['Regressors for ' Sess{s}.name{i}])
axis tight


% Basis set and peristimulus sampling
%---------------------------------------------------------------------------
subplot(3,2,4)
t    = [1:size(Sess{s}.bf{i},1)]*X.dt;
pst  = Sess{s}.pst{i};
plot(t,Sess{s}.bf{i},pst,0*pst,'.','MarkerSize',16)
xlabel('time (secs}')
title('Basis set and peristimulus sampling')
axis tight
grid on

% onsets and parametric modulation
%---------------------------------------------------------------------------
subplot(3,2,6)
plot(Sess{s}.ons{i},Sess{s}.para{i},'.','MarkerSize',8)
xlabel('time (secs)')
title('trial specific onsets or causes')
grid on

% textual data
%---------------------------------------------------------------------------
subplot(3,2,5)
cla
q     = 0.1;
text(0,1 - 1*q,'Design','FontWeight','bold')
text(0,1 - 2*q,X.DesN{1},'FontSize',10)
text(0,1 - 3*q,'Basis functions','FontWeight','bold')
text(0,1 - 4*q,X.DesN{2},'FontSize',10)
text(0,1 - 6*q,sprintf('%0.0f ms time bins',1000*X.dt),'FontSize',10)
text(0,1 - 7*q,sprintf('%0.0f of %0.0f session[s],',s,length(Sess)),'FontSize',10)
text(0,1 - 8*q,sprintf('%0.0f of %0.0f trial[s],',i,length(Sess{s}.name)),'FontSize',10)
text(0,1 - 10*q,Sess{s}.name{i},'FontWeight','bold')
axis off

return


% User interface for spm_svd.m - eigenimage analysis
% FORMAT spm_svd_ui
%___________________________________________________________________________
%
% spm_svd_ui prompts for the selection of XA.mat and then computes
% eigenimages and component scores using SVD.  spm_svd.m displays the 
% selected eigenimage in the results window
%
% The eigenimages can be based on the group data or on any single subject.
% The voxels entered into the analysis are those with an F value with
% p < 0.05 (uncorrected)
%
%__________________________________________________________________________
% %W% Karl Friston %E%

% find and clear Interactive window
%---------------------------------------------------------------------------
Finter = spm_figure('FindWin','Interactive');
if isempty(Finter), Finter=spm('CreateIntWin'); end
spm_clf(Finter)
spm_clf('Graphics')
set(Finter,'Name','Eigenimages')

% get adjusted data to analyse
%---------------------------------------------------------------------------
tmp = spm_get(1,'.mat','select adjusted data you wish to analyse','XA');
global CWD
CWD = strrep(tmp,'/XA.mat','');
load([CWD,'/SPM'])
load([CWD,'/XA'])
load([CWD,'/XYZ'])

% create data matrix {M}
%---------------------------------------------------------------------------
M      = 0;
if size(B,2)>1

	% balanced design
	%-------------------------------------------------------------------
	if ~any(diff(sum(B > 0)))
		if spm_input('all subjects {blocks} ?',1,'yes|no',[1 0])
			for i = 1:size(B,2)
				M = M + XA((B(:,i) > 0),:);
			end
		else
			d = spm_input('which subject[s] ?',2);
			for i = 1:length(d)
				M = M + XA((B(:,d(i)) > 0),:);
			end
		end
	else
		d = B(:,spm_input('which subject ?',2)) >  0;
		M = XA(d,:);
	end
else
	M   = XA;
end
clear XA

set(Finter,'Pointer','Watch')


% Singlar Value Decompostion
%----------------------------------------------------------------------------
M       = M - ones(size(M,1),1)*mean(M);
[u s v] = svd(M',0);
s       = diag(s).^2;
s       = length(s)*s/sum(s);

% display
%----------------------------------------------------------------------------
e    = 1;
spm_svd

% setup buttons to move between eigenimages
%----------------------------------------------------------------------------
c    = ['set(Finter,''Pointer'',''watch''); e = e + 1;; e = min([length(s) e]);',...
	'spm_svd; set(Finter,''Pointer'',''arrow'');'];
uicontrol(Finter,'Style','Pushbutton','Position',[80,200,100,30],...
	'String','next','Callback',c);

c    = ['set(Finter,''Pointer'',''watch''); e = e - 1; e = max([1 e]);'...
	'spm_svd; set(Finter,''Pointer'',''arrow'');'];
uicontrol(Finter,'Style','Pushbutton','Position',[220,200,100,30],...
	'String','previous','Callback',c);

% bring Interactive window to the front and set it's pointer
%----------------------------------------------------------------------------
set(Finter,'Pointer','Arrow')
figure(Finter)

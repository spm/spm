function spm_estimate_ui
% Callback function for estimation with spm_spm
% FORMAT spm_estimate_ui
%___________________________________________________________________________
% %W% Karl Friston %E%

% check for previous results
%---------------------------------------------------------------------------
if exist(fullfile('.','SPM.mat'),'file') == 2
	tmp    = ...
	spm_input({'Current directory contains existing SPMstats files:',...
	'        SPMstats results files (inc. SPM.mat)',...
	['(pwd = ',pwd,')'],' ',...
	'Continuing will overwrite existing results!'},...
	1,'bd','stop|continue',[0,1],1);
else
	tmp = 1;
end

% OLS estimation
%---------------------------------------------------------------------------
if tmp
	tmp = load(spm_get(1,'SPMcfg.mat','Select SPMcfg.mat...'));
	if isfield(tmp,'Sess') & ~isempty(tmp.Sess)
		Sess  = tmp.Sess; 	
		xsDes = tmp.xsDes;
		spm_spm(tmp.VY,tmp.xX,tmp.xM,tmp.F_iX0,Sess,xsDes);

	elseif isfield(tmp,'xC')
		xC    = tmp.xC;
		xsDes = tmp.xsDes;
		spm_spm(tmp.VY,tmp.xX,tmp.xM,tmp.F_iX0,xC,xsDes);
	end
end

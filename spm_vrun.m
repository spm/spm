function spm_vrun(Vcode,SPMarg)
% Utility function to run various SPM versions (ML4 & ML5 compatible)
% FORMAT spm_vrun(Vcode,SPMarg)
%
% Vcode   - version identifier code
% SPMarg  - argument to be passed to spm.m
%_______________________________________________________________________
%
% spm_vrun is a standalone MatLab 4/5 compatible function which allows
% the user to choose between multiple versions of SPM.
%
% The path to the specified SPM version is prepended to the current
% path, and `spm` called. Note that since the SPM path is prepended,
% standard SPM functions will take precedence over functions in other
% directories on the path.
% 
% If no particular version of SPM is specified, a default is chosen.
% 
% Systems administrators should code the appropriate definitions into
% the Parameters section in the main body of the code.
%
% See also: spm_choose.m, a GUI for choosing SPM versions, but note that
%		spm_choose is Matlab version dependent. 
%
%_______________________________________________________________________
% %W% Andrew Holmes %E%

%-Parameters - up to 4 versions may be defined
%=======================================================================
Vcodes = str2mat(	'spm_devel',...			%-1
			'spm99',...			%-2
			'spm99r',...			%-3
			'spm97',...			%-4
			'spm96',...			%-5
			'spm95'	);			%-6

Vpaths = str2mat(	'/local/spm/spm_devel',...	%-1
			'/local/spm/spm99_fil',...	%-2
			'/local/spm/spm99_rel',...	%-3
			'/local/spm/spm97_fil',...	%-4
			'/local/spm/spm96_fil',...	%-5
			'/local/spm/spm95_fil'	);	%-6

VMLver = [		5,...				%-1
			5,...				%-2
			5,...				%-3
			4,...				%-4
			4,...				%-5
			4	];			%-6

DefV = [0,0,0,4,1];
		

%-Sort out version to use...
%=======================================================================
MLver = version;
MLver = eval(MLver(1));

if nargin<1
	V = DefV(MLver(1));
else
	V = 0; v = 0;
	while V==0 & v<size(Vcodes,1)
		v=v+1;
		if strcmp(lower(Vcode),deblank(Vcodes(v,:))), V=v; end
	end
	if V==0, error('Unknown version!'), end
end

if exist([deblank(Vpaths(V,:)),'/spm.m'])~=2
	error(sprintf('Invalid path for %s\n\t(File %s/spm.m not found)',...
		upper(deblank(Vcodes(V,:))),deblank(Vpaths(V,:))))
end

if VMLver(V)~=MLver(1)
	error(sprintf('Wrong MatLab for %s, use version %d',...
		upper(deblank(Vcodes(V,:))),VMLver(V)))
end


%-Set path and run...
%=======================================================================
fprintf('spm_vrun: Prepending %s to path...\n\n',deblank(Vpaths(V,:)))

path(deblank(Vpaths(V,:)),path)

if nargin>1
	eval('spm(SPMarg)')
else
	eval('spm')
end

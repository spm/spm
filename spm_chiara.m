function spm_chiara
% Contrast to apply to parameter estimates in nonorthogonal designs
% FORMAT spm_chiara
%___________________________________________________________________________
% %W% Karl Friston %E%
SCCSid = '%I%';

% get parameters
%---------------------------------------------------------------------------
[Finter,Fgraph,CmdLine] = spm('FnUIsetup','Chiara',1);
PWD  = pwd;
CWD  = spm_str_manip(spm_get(1,'SPM.mat','Select SPM.mat for analysis'),'H');
load([CWD,'/SPM'])
cd(CWD)

% get contrasts
%---------------------------------------------------------------------------
m      = size(xX.X,2);
str    = sprintf('contrast[s] n x 1 - %i',m);
C      = spm_input(str,1);
[u v]  = size(C);


% get orthogonlization order
%---------------------------------------------------------------------------
if spm_input('orthogonlization scheme',1,'b','no|yes',[0 1])

	h     = spm_input('orthogonlization order e.g. 1 2 3');

	% orthonalization matrix
	%-------------------------------------------------------------------
	X     = xX.X;
	[m n] = size(X);
	d     = [1:n];
	c     = zeros(u,n);
	c(:,[1:v]) = C;
	C     = c;
	for i = 1:length(h)
		d         = d(d ~= h(i));
		d         = [h(i) d];
	end
	h     = d;
	H     = eye(n);
	for i = 1:length(h)
		d         = [1:n];
		d         = d(d ~= h(i));
		D         = eye(n);
		D(d,h(i)) = -pinv(X(:,d))*X(:,h(i));
		X         = X*D;
		H         = H*D;
	end
	W     = C*inv(H);
else
	W     = C;
end


% compute compound
%---------------------------------------------------------------------------
spm('Pointer','Watch')

%-Create handle template for output images as 16bit Int16's
%---------------------------------------------------------------------------
Vo = struct(	'fname',	'',...
		'dim',		[Vbeta(1).dim(1:3),4],...
		'mat',		Vbeta(1).mat,...
		'pinfo'	,	[1,0,0]',...
		'descrip',	'');

%-Loop over contrasts
%---------------------------------------------------------------------------
for i = 1:u

	%-Implement weighted sum by weighting scalefactors in image handles
	%-------------------------------------------------------------------
	w  = W(i,:);
	Q  = find(abs(w)> 1e-8);
	w  = w(Q); wV = Vbeta(Q);
	for j = 1:length(Q), wV(j).pinfo(1,:) = wV(j).pinfo(1,:)*w(j); end

	%-Write header
	%-------------------------------------------------------------------
	Vo.fname   = sprintf('contrast%d.img',i);
	Vo.descrip = sprintf('contrast');
	Vo.pinfo   = [1,0,0]';
	spm_create_image(Vo);

	%-Compute & rewrite header scalefactor
	%-------------------------------------------------------------------
	Vo.pinfo(1)  = spm_add(wV,Vo,'m');
	spm_create_image(Vo);
end


%-Finished
%---------------------------------------------------------------------------
spm('Pointer','Arrow')
spm_clf(Finter)
cd(PWD)

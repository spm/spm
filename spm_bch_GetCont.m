function status = spm_bch_GetCont(batch_mat,n_con);
%
% create an xCon.mat file from batch contrast, append to 
% previous xCon if it exists
%
% %W% Jean-Baptiste Poline & Stephanie Rouquette %E%

%--------------------------------------------------
%- get xCon if exist; We are in Wdir.
%- WDir = char(spm_input('batch',batch_mat,{},'WDir'));
%--------------------------------------------------


%- initialise status
status.str = '';
status.ok = 1;

swd = pwd;


if exist(fullfile('.','xCon.mat'),'file'), 
	load('xCon.mat'), 
	lxCon = length(xCon);
else, 
	xCon = spm_FcUtil('FconFields');
	lxCon = 0;
end



if exist(fullfile('.','SPM.mat'),'file'), 
	try 
		load(fullfile('.','SPM.mat'),'xX');	
	catch 
		str = ['cannot open ' fullfile('.','SPM.mat') ...
				' file in spm_bch_GetCont ' swd];
		warning(str);
		status.str = str;
		status.ok = 0;
		return;
	end
else 
	str = ['cannot find ' fullfile('.','SPM.mat') ...
				' file in spm_bch_GetCont ' swd];
	warning(str);
	status.str = str;
	status.ok = 0;
	return;
end

%- get contrast to create
%--------------------------------------------------

names = spm_input('batch',batch_mat,{'contrastes',n_con},'names');
types = spm_input('batch',batch_mat,{'contrastes',n_con},'types');
values = spm_input('batch',batch_mat,{'contrastes',n_con},'values');

%- check that the lengths are identical ? NO, should be done in spm_bch_bchmat

len = [length(names) length(types) length(values)];
sX = xX.xKXs;



for n=1:min(len)
   contrast = spm_FcUtil('Set',names{n}, types{n}, 'c', values{n}', sX);
   iFc2 = spm_FcUtil('In', contrast, sX, xCon);
   if ~iFc2, 
      xCon(length(xCon)+1) = contrast;
   else 
      %- 
      fprintf('\ncontrast %s (type %s) already in xCon', names{n}, types{n});
   end
end

try
	save('xCon.mat','xCon')
catch
	str = ['Can''t write xCon.mat to the results directory: ' swd];
	warning(str);
	status.str = str;
	status.ok = 0;
end

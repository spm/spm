function status = spm_bch_GetCont;
%
% The BCH gobal variable is used for spm_input in batch mode : 
%    BCH.bch_mat 
%    BCH.index0  = {'contrasts',index_of_Analysis};
%
% Create an xCon.mat file from contrasts described in the mfile, 
% append to previous xCon if xCon already exists.
% Working directory where the xCon structure will be saved
% is specified in the top level m-file. 
%
% %W%  Jean-Baptiste Poline & Stephanie Rouquette  %E%
%---------------------------------------------------------------

%- initialise status
status.str = '';
status.err = 0;

swd = pwd;


%---------------------------------------------------------------
if exist(fullfile('.','xCon.mat'),'file'), 
	load('xCon.mat'), 
	lxCon = length(xCon);
else, 
	xCon = spm_FcUtil('FconFields');
	lxCon = 0;
end

%---------------------------------------------------------------
if exist(fullfile('.','SPM.mat'),'file'), 
	try 
	   load(fullfile('.','SPM.mat'),'xX');	
	catch 
	   str = ['cannot open ' fullfile('.','SPM.mat') ...
                  ' file in spm_bch_GetCont ' swd];
	   warning(str);
	   status.str = str;
	   status.err = 1;
	   return;
	end
else 
	str = ['cannot find ' fullfile('.','SPM.mat') ...
               ' file in spm_bch_GetCont ' swd];
	warning(str);
	status.str = str;
	status.err = 2;
	return;
end

%- get contrast to create from mat file (in global BCH) 
%---------------------------------------------------------------

names  = spm_input('batch',{},'names');
types  = spm_input('batch',{},'types');
values = spm_input('batch',{},'values');

%- check that the lengths are identical ? 
%- NO, this should be done in spm_bch_bchmat

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
	status.err = 3;
end

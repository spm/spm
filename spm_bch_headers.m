function status = spm_bch_headers
%
% Use the BCH gobal variable for spm_input : 
%    BCH.bch_mat 
%    BCH.index0  = {'headers',index_of_Analysis};
%
% This function will modify the headers 
%
% if ORIGIN is [0 0 0] in the header => set to  (DIM(1:3)+1)/2;
% if ORIGIN is [1 1 1] in the header => set to  (DIM(even)/2)+1;
% if ORIGIN is [1 1 1] in the header => set to  (DIM(odd)+1)/2;
%
% do_mat    : create a .mat file 
% origoff   : add origoff to ORIGIN
%
% Programmers Guide :
% cheching on the arguments read in batch_mat should be done
% in spm_bch_bchmat.m
% %W%  Jean-Baptiste Poline & Stephanie Rouquette  %E%
%---------------------------------------------------------------


%- initialise status
status.str = '';
status.err = 0;

files = spm_input('batch',{},'files');

%------- no files ------- 
if ~size(files,1), return, end;

DIM 		= spm_input('batch',{},'DIM');
VOX 		= spm_input('batch',{},'VOX');
SCALE 		= spm_input('batch',{},'SCALE');
TYPE 		= spm_input('batch',{},'TYPE');
OFFSET 		= spm_input('batch',{},'OFFSET');
ORIGIN 		= spm_input('batch',{},'ORIGIN');
DESCRIP 	= spm_input('batch',{},'DESCRIP');
origoff 	= spm_input('batch',{},'origoff');
do_mat 		= spm_input('batch',{},'do_mat');

str = '[';
if isempty(DIM), str = [str 'DIM ']; else str = [str 'dim ']; end
if isempty(VOX), str = [str 'VOX ']; else str = [str 'vox ']; end
if isempty(SCALE), str = [str 'SCALE ']; else str = [str 'scale ']; end
if isempty(TYPE), str = [str 'TYPE ']; else str = [str 'type ']; end
if isempty(OFFSET), str = [str 'OFFSET ']; else str = [str 'offset ']; end
if isempty(ORIGIN), str = [str 'ORIGIN ']; mptyORIG = 1; 
else str = [str 'origin ']; mptyORIG = 0; end
if isempty(DESCRIP), str = [str 'DESCRIP ']; else str = [str 'descrip ']; end
str = [str ']'];
	
for k=1:size(files,1)

	P = files(k,:);
	
	%----------  read the empty parameters from header
	% [str ' = spm_hread('''  deblank(P)  ''');']
	try 
	   eval([str ' = spm_hread('''  deblank(P)  ''');'])
	catch
		status.str = [status.str 'pb spm_hread ' P];
		status.err = 1; 
	end
	%----------  set default ORIGIN and VOX if necessary
	if all(ORIGIN == 0) & mptyORIG , ORIGIN = (DIM(1:3)+1)/2; end;
	if all(ORIGIN == 1) & mptyORIG, 
	   odd = rem(DIM,2);
	   even = ~odd;
	   ORIGIN(find(even)) = (DIM(find(even))/2)+1;
	   ORIGIN(find(odd)) = (DIM(find(odd))+1)/2;
	end;
	if all(VOX == 0), VOX = [1 1 1]; end;
	
	%----------  add origoff to ORIGIN
	if ~isempty(origoff), ORIGIN = ORIGIN + origoff; end
	
	if isstr(TYPE), TYPE = spm_type(TYPE); end
	try
		spm_hwrite(P,DIM,VOX,SCALE,TYPE,OFFSET,ORIGIN,DESCRIP);
	catch
		status.str = [status.str 'pb in spm_hwrite ' P];
		status.err = 2;
	end

	%----------  if asked to write .mat file 
	if do_mat
	   matname = [spm_str_manip(P,'sd') '.mat'];
	
	   % if (exist(matname) == 2), 
	   % the file is overwritten if permissions ok
	
	   offs = -VOX.*ORIGIN;
	   M   = [VOX(1) 0 0 offs(1) ; ...
		  0 VOX(2) 0 offs(2) ; ...
		  0 0 VOX(3) offs(3) ; 0 0 0 1];
	   try
		eval(['save ' matname ' M -v4']);
	   catch
		status.str = [status.str 'pb in saving ' matname];
		status.err = 3;
	   end

	end % if do_mat

end %- for k=1:size(files,1)



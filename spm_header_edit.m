function sts = spm_header_edit(arg1,arg2,arg3)
% Edit/Create image headers.
% FORMAT spm_header_edit
%
% OR sts = spm_header_edit('set',P,arg3)
% P    - filename of image
% arg3 - parameters to change.
%        - this is a string which is evaluated.
%          eg. 'DIM=[128 128 43]; DESCRIP=''an image'';'
% sts  - 0 means sucess / 1 means failure
%_______________________________________________________________________
%
% The "set" options are for setting values in the headers.
% The "unset" options are for undoing the changes made by "set".
% The "APPLY" option prompts the user for the images which are to be
% modified, before proceeding to modify those headers.
%
% Elements in the headers are preserved, except for the ones which
% are specified by the "set" options.
%_______________________________________________________________________
% %W% John Ashburner %E%

if nargin == 0

	% Act as user interface
	%=======================================================================

	opts = ['  set Image Dimensions - pixels in x, y & z|'...
		'  set Voxel Dimensions - mm in x, y & z|'...
		'  set Scalefactor|'...
		'  set Data Type|'...
		'  set Offset into file - bytes|'...
		'  set Origin - x, y & z coordinate of voxel|'...
		'  set Image description|'...
		'  APPLY to images|'...
		'  QUIT'];

	callbacks = str2mat(...
		['a2 = spm_input(''Image Dimensions [x y z]'',2,''s'');'...
		 'str1 = [''DIM =['' a2 ''];''];'],...
		['a2 = spm_input(''Voxel Dimensions [x y z]'',2,''s'');'...
		 'str2 = [''VOX =['' a2 ''];''];'],...
		['a2 = spm_input(''Scalefactor'',2,''s'');'...
		 'str3 = [''SCALE =['' a2 ''];''];'],...
		['a2 = spm_input(''Data Type'',2,''m'','...
		 ' ''datatype - unsigned char|datatype - signed short|datatype - signed int|datatype - float|datatype - double'','...
		 ' [2 4 8 16 64]);'...
		 'str4 = a2;'],...
		['a2 = spm_input(''Offset into file'',2,''s'');'...
		 'str5 = [''OFFSET =['' a2 ''];''];'],...
		['a2 = spm_input(''Origin [x y z]'',2,''s'');'...
		 'str6 = [''ORIGIN =['' a2 ''];''];'],...
		['a2 = spm_input(''Description'',2,''s'');'...
		 'str7 = [''DESCRIP =['''''' a2 ''''''];''];']);

	posns = [1 (find(opts == '|') + 1)];

	for i=1:(length(posns)-2)
		eval(['str' num2str(i) '='''';']);
	end

	promptstr = 'Options';

	while(1)
		a1 = spm_input(promptstr,1,'m', opts);

		if a1 == length(posns)
			% Quit
			%=======================================================================
			break;

		elseif a1 == length(posns)-1
			% Make requested modifications
			%=======================================================================

			% combine strings together
			%-----------------------------------------------------------------------
			str = '';
			for i=1:(length(posns)-2)
				eval(['str = [str str' num2str(i) '];']);
			end

			% specify images
			%-----------------------------------------------------------------------
			P = spm_get(Inf, 'img','Images to apply changes to');


			% Do the changes
			%-----------------------------------------------------------------------
			spm_progress_bar('Init',size(P,1),'Editing Headers',...
				'Headers Complete');
			for i=1:size(P,1)
				s = spm_header_edit('set',deblank(P(i,:)),str);
				if (s == 1)
					promptstr = 'It failed - try again';
					break;
				end
				spm_progress_bar('Set',i);
			end
			spm_progress_bar('Clear');


		elseif opts(posns(a1)) == 'u'

			% Clear requested modifications
			%-----------------------------------------------------------------------
			eval(['str' num2str(a1) '='''';']);

			opts(posns(a1)  ) = ' ';
			opts(posns(a1)+1) = ' ';

			promptstr = 'Options';
		else
			% Get requested modifications
			%-----------------------------------------------------------------------
			opts(posns(a1)  ) = 'u';
			opts(posns(a1)+1) = 'n';

			eval(callbacks(a1,:));
			promptstr = 'Options';
		end
	end

elseif strcmp(arg1,'set')

	% Read  header
	%-----------------------------------------------------------------------
	[DIM VOX SCALE TYPE OFFSET ORIGIN DESCRIP] = spm_hread(arg2);

	% Make required modifications
	%-----------------------------------------------------------------------
	eval(arg3, ['disp([''Cant evaluate "'' arg3 ''".'']);sts=1;return']);
	if (prod(size(DIM))~=3 | prod(size(VOX))~=3 ...
		| prod(size(SCALE ))~=1 | prod(size(TYPE  ))~=1 ...
		| prod(size(OFFSET))~=1 | prod(size(ORIGIN))~=3)
		disp(['Problem with parameters (' arg3 ')' ]);
		sts = 1;
		return;
	end

	% Write header with modifications
	%-----------------------------------------------------------------------
	eval('spm_hwrite(arg2,DIM,VOX,SCALE,TYPE,OFFSET,ORIGIN,DESCRIP);','sts=1;');
	sts = 0;
else
	error('Incorrect usage');
end

function [Files,Dirs] = spm_list_files_expand(WDir,Filter,opts)
% List files, directories and volumes within analyze files
% FORMAT [Files,Dirs] = spm_list_files_expand(WDir,Filter,opts)
% Dir    - directory to list
% Filter - e.g. '*.img'
%        - alternatively, FILTER can take the value of 'IMAGE', in
%          order to select files with extension '.img', '.mnc', '.v'.
% opts   - options (an optional cell array)
% Files  - filenames
% Dirs   - directories
%
% If the global variable defaults.analyze.multivol is set, then
% spm_list_files_expand will list volumes within each Analyze image.
% This can be disabled by setting opts to {'noexpand'}.
%_______________________________________________________________________
% @(#)spm_list_files_expand.m	1.6 John Ashburner 03/12/10

global defaults

if nargin<2, Filter = '*'; end;
if nargin<1, Wdir   = '.'; end;

Filter = deblank(Filter);

if strcmp(Filter,'IMAGE') | (length(Filter)>5 & strcmp(Filter((end-4):end),'IMAGE')),
	if length(Filter)>5,
		Filter = [Filter(1:(end-5)) '*'];
	else,
		Filter = '*';
	end
	[Files,Dirs] = spm_list_files(WDir,Filter);
	t=cell(size(Files,1),1);
	for i=1:size(Files,1),
		[pth,nam,ext] = fileparts(deblank(Files(i,:)));
		switch ext,
		case {'.img','.mnc','.v'},
			t{i} = Files(i,:);
		end;
	end;
	Files = strvcat(t{:});
else,
	[Files,Dirs] = spm_list_files(WDir,Filter);
end;

if  isempty(defaults) | ~isfield(defaults,'analyze') |...
	 ~isfield(defaults.analyze,'multivol') | ~defaults.analyze.multivol |...
	(nargin>=3 & ~isempty(strmatch('noexpand',lower(opts)))),
	return;
end;

t=cell(size(Files,1),1);
for i=1:size(Files,1),
	[pth,nam,ext] = fileparts(deblank(Files(i,:)));
	pth           = WDir;
	if strcmp(ext,'.img'),
		hname = fullfile(pth,[nam '.hdr']);
		fp    = fopen(hname,'r');
		if fp ~= -1,
			fseek(fp,40,'bof');
			dim = fread(fp,5,'short');
			fclose(fp);
			if dim(1)<0 | dim(1)>15, % Appears to be other-endian
				if spm_platform('bigend'), fp = fopen(hname,'r','ieee-le');
				else,                      fp = fopen(hname,'r','ieee-be'); end;
				fseek(fp,40,'bof');
				dim = fread(fp,5,'short');
				fclose(fp);
			end;
			if dim(1)<0 | dim(1)>15 | dim(5)<=1,
				t{i} = Files(i,:);
			else,
				t{i} = [repmat([deblank(Files(i,:)) ','],dim(5),1) num2str((1:dim(5))')];
			end;
		else,
			t{i} = Files(i,:);
		end;
	else,
		t{i} = Files(i,:);
	end;
end;
Files = strvcat(t{:});

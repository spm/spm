function spm_append(data,fname,options)
% appends columns to a matrix file stored on disk
% FORMAT spm_append(data,fname,options)
% data       - the data to be appended.
% fname	     - name of data file.
%              default file name is ``dat.mad''.
% options(1) - datatype used to represent the data
%              default datatype is 2 (uint8)
%              - see spm_type
%
% The columns in the file are represented by the specified datatype, but
% with additional scalefactors and dc-offsets.  The scalefactos and offsets
% apply to all elements in a column, but take different values for different
% columns.
% 
% See also spm_extract.
%___________________________________________________________________________
% %W% John Ashburner %E%

if nargin<3,
	if nargin<1, return; end;
	if nargin<2, fname = 'dat.mad'; end;
	options = spm_type('uint8');
end;
if isempty(data), return; end;
if length(size(data))~=2, error('Data must have no more than two dimensions.'); end;

data = real(data);

typ = options(1);
hdrlen = 256;
MGC = 141098;

fp = my_fopen(fname,'r+','silent');
if fp.ptr==-1,
	% File doesn't exist, so create it.

	if strcmp(spm_type(typ),'unknown'),
		error(['Unknown data type (' num2str(typ) ').']);
	end;
	if spm_type(typ,'swapped'),
		error(['Cannot work with byteswapped datatypes (' num2str(typ) ').']);
	end;

	fp = my_fopen(fname,'w+','error');
	my_fwrite(fp,MGC,'uint32');                  % Magic number identifying file type
	my_fwrite(fp,size(data), 'uint32');          % Matrix dimensions
	my_fwrite(fp,typ, 'int32');                  % Datatype
	my_fwrite(fp,zeros(hdrlen-3*4,1), 'uchar');  % Padding
	dim = [size(data,1) 0];                      % Existing file dimensions
else,
	% Read header information from existing file.

	mgc = my_fread(fp,1,'uint32');
	if mgc ~= MGC,
		my_fclose(fp);
		error(['"' fname '" appears to be the wrong kind of file.']);
	end;
	dim = my_fread(fp,2,'uint32');
	if dim(1) ~= size(data,1),
		my_fclose(fp);
		error('Incompatible data size.');
	end;
	typ = my_fread(fp,1,'int32');
end;


bits   = spm_type(typ,'bits');
len    = ceil(dim(1)*bits/64)/8*64; % Align to next highest 8 bit boundary
rmndr  = len-dim(1)*bits/8;         % Used so that data is aligned to 8 bit boundaries.
typstr = spm_type(typ);

% Determine scalefactors etc where necessary
range  = [spm_type(typ,'minval') spm_type(typ,'maxval')];
if all(finite(range)),
	mx     = max(data,[],1);
	mn     = min(data,[],1);
	scale  = (mx-mn)/(range(2)-range(1));
	scale(find(scale==0)) = 1;
	off    = mn - scale*range(1);
	rescale = 1;
else,
	scale = ones(1,size(data,2));
	off   = zeros(1,size(data,2));
	rescale = 0;
end;

for i=1:size(data,2),
	% write a column.
	ofset = hdrlen + (len + 2*8)*(dim(2)+i-1);
	my_fseek(fp,ofset,'bof');
	my_fwrite(fp,[scale(i) off(i)], 'float64');
	if rescale, dt = round((data(:,i)-off(i))/scale(i));
	else,       dt = data(:,i); end;
	my_fwrite(fp,dt, typstr);
	if rmndr ~= 0,
		my_fwrite(fp,zeros(rmndr,1), 'uint8');
	end;
end;

% Update matrix dimensions
my_fseek(fp,4*2,'bof');
my_fwrite(fp,dim(2)+size(data,2), 'uint32');
my_fclose(fp);
return;


function fp = my_fopen(fname,perm,flag)
if nargin < 2,
	flag = 'error';
end;
fp = struct('ptr',-1,'fname',fname,'perm',perm);
fp.ptr = fopen(fp.fname,fp.perm,'ieee-be');  % Make everything big-end.
if (fp.ptr == -1) & (strcmp(flag,'error')),
	error(sprintf('Can''t open "%s" (perm "%s")\n%s\n', fname, perm,...
		'Check that you have permission.'));
end;
return;

function my_fwrite(fp,a,prec)
count = fwrite(fp.ptr,a,prec);
if count ~= prod(size(a)),
	er=ferror(fp.ptr);
	my_fclose(fp);
	error(sprintf(...
		'Problems writing to file "%s"\nThe error was:\n"%s"\n%s\n',...
		fp.fname,er, 'Check your disk space.'));
end;
return;

function a = my_fread(fp,count,prec)
a = fread(fp.ptr,count,prec);
if count ~= prod(size(a)),
	er=ferror(fp.ptr);
	my_fclose(fp);
	error(sprintf(...
		'Problems reading from file "%s"\nThe error was:\n"%s"\n%s\n',...
		fp.fname,er, 'Check what has already gone wrong.'));
end;
return;

function my_fseek(fp,offset,origin)
fseek(fp.ptr,0,-1);
sts=fseek(fp.ptr,offset,origin);
if sts == -1,
	er=ferror(fp.ptr);
	my_fclose(fp);
	error(sprintf(...
		'Problems seeking through "%s".  The error was:\n"%s"\n%s\n',...
		fp.fname,er, 'Check what has already gone wrong.'));
end;
return;

function my_fclose(fp)
fclose(fp.ptr);
return;

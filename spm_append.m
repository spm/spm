function spm_append(fname,data)
% appends columns to a matrix file stored on disk
% FORMAT spm_append(fname,data)
% fname	- name of data file.
% data  - the data to be appended.
%
% See also spm_extract.
%___________________________________________________________________________
% %W% John Ashburner %E%

if isempty(data), return; end;

hdrlen = 256;
MGC = 270198;

fp = my_fopen(fname,'r+','silent');
if fp.ptr==-1,
	fp = my_fopen(fname,'w+','error');
	
	my_fwrite(fp,MGC,'uint32');
	my_fwrite(fp,size(data), 'uint32');
	my_fwrite(fp,zeros(hdrlen-3*4,1), 'uchar');
	dim = [size(data,1) 0];
else,
	mgc = my_fread(fp,1,'uint32');
	if mgc ~= MGC,
		my_fclose(fp);
		error(['" fname "' appears to be the wrong kind of file.']);
	end;
	dim = my_fread(fp,2,'uint32');
	if dim(1) ~= size(data,1),
		my_fclose(fp);
		error('Incompatible data size.');
	end;
end;

mx = max(data,[],1);
mn = min(data,[],1);
scale = (mx-mn)/(2^16-1);
off = mn;

len    = 8*ceil(dim(1)*2/8);
rmndr  = rem(dim(1),4);
for i=1:size(data,2),
	ofset = hdrlen + (len + 2*8)*(dim(2)+i-1);
	my_fseek(fp,ofset,'bof');
	my_fwrite(fp,[scale(i) off(i)], 'float64');
	dt = round((data(:,i)-off(i))/scale(i));
	my_fwrite(fp,dt, 'uint16');
	if rmndr ~= 0,
		my_fwrite(fp,zeros(rmndr,1), 'uint16');
	end;
end;

my_fseek(fp,4*2,'bof');
my_fwrite(fp,dim(2)+size(data,2), 'integer*4');
my_fclose(fp);
return;


function fp = my_fopen(fname,perm,flag)
if nargin < 2,
	flag = 'error';
end;
fp = struct('ptr',-1,'fname',fname,'perm',perm);
fp.ptr = fopen(fp.fname,fp.perm); %,'ieee-be');
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

function data = spm_extract(fname,index)
% extracts columns from a matrix file stored on disk
% FORMAT data = spm_extract(fname,index)
% fname - name of data file.
% index - the indexes of the required columns.
%         - if index is Inf or NaN, then the whole
%           matrix is extracted.
% data  - the extracted columns of the matrix.
%
% OR     dim = spm_extract(fname)
% fname - name of data file.
% dim   - the dimensions of the stored matrix.
%
% See also spm_append
%___________________________________________________________________________
% %W% John Ashburner %E%

hdrlen = 256;
MGC = 141098;

fp = my_fopen(fname,'r','error');
mgc = my_fread(fp,1,'uint32');
if mgc ~= MGC,
	my_fclose(fp);
	error(['"' fname '" appears to be the wrong kind of file.']);
end;

dim    = my_fread(fp,2,'uint32')';
typ    = my_fread(fp,1,'int32')';
typstr = spm_type(typ);
bits   = spm_type(typ,'bits');

if nargin<2,
	% Return dimensions of the matrix.
	data = dim;
else,
	% Return the data itself.
	len    = ceil(dim(1)*bits/64)/8*64;
	if prod(size(index))==1 & ~finite(index(1)),
		data = zeros(dim(1),dim(2));
		for i=1:dim(2),
			my_fseek(fp,hdrlen + (len + 2*8)*(i-1),'bof');
			sca = my_fread(fp,2, 'float64');
			tmp = my_fread(fp,dim(1), typstr);
			data(:,i) = tmp*sca(1)+sca(2);
		end;
	else,
		if any((index>dim(2)) | (index<1)),
			my_fclose(fp);
			error(['Trying to read columns from "' fname '" that don''t exist.']);
		end;

		data = zeros(dim(1),prod(size(index)));
		for i=1:prod(size(index)),
			ind=index(i);
			my_fseek(fp,hdrlen + (len + 2*8)*(ind-1),'bof');
			sca = my_fread(fp,2, 'float64');
			tmp = my_fread(fp,dim(1), typstr);
			data(:,i) = tmp*sca(1)+sca(2);
		end;
	end;
end;
my_fclose(fp);
return;


function fp = my_fopen(fname,perm,flag)
if nargin < 2,
	flag = 'error';
end;
fp = struct('ptr',-1,'fname',fname,'perm',perm);
fp.ptr = fopen(fp.fname,fp.perm,'ieee-be');
if (fp.ptr == -1) & (strcmp(flag,'error')),
	error(sprintf('Can''t open "%s" (perm "%s")\n%s\n', fname, perm,...
		'Check that you have permission.'));
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

function xa = spm_readXA(element, filename)

% Read columns from an XA.mat file
% FORMAT xa = spm_readXA(element, filename)
% element  - columns to read.
% filename - name of the XA.mat file
% xa       - the data read from XA.mat
%
% Note: that if the file is not an XA.mat file - then columns from the
% first variable stored in the file `filename' will be read.
% The data must be stored as real, full, etc...
%_______________________________________________________________________
% %W% John Ashburner %E%

%-Default XA.mat file
%-----------------------------------------------------------------------
if (nargin == 1)
	global CWD;
	filename = [CWD '/XA.mat'];
end

%-Open file
%-----------------------------------------------------------------------
fp = fopen(deblank(filename),'r');
if (fp == -1)
	error(['Can''t open ' deblank(filename)]);
end

%-Read dimensions
%-----------------------------------------------------------------------
fseek(fp,4,'bof');
dim = fread(fp,2,'uint32');

%-Range check
%-----------------------------------------------------------------------
element = round(element);
if any(element>dim(2) | element<1)
	error('Elements out of range');
end

%-Need to skip over the variable name
%-----------------------------------------------------------------------
fseek(fp,16,'bof');
nlen=fread(fp,1,'int32');

%-Read the data
%-----------------------------------------------------------------------
xa = zeros(dim(1),length(element));
for i=1:length(element)
	fseek(fp,20+nlen+(element(i)-1)*8*dim(1),'bof');
	xa(:,i) = fread(fp,dim(1),'double');
end

fclose(fp);



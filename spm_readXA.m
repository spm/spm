function xa = spm_readXA(element, filename)

% Read columns from an XA.mat file
% FORMAT xa = spm_readXA(element, filename)
% element  - columns to read.
% filename - name of the XA.mat file
% xa       - the data read from XA.mat
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

%-Check that it is (probably) in the right format
%-----------------------------------------------------------------------
fseek(fp,20,'bof');
nam = fread(fp,3,'char');
if ~all(nam == [88 65 0]')
	error('This is not a proper XA.mat');
end

%-Read the data
%-----------------------------------------------------------------------
xa = zeros(dim(1),length(element));
for i=1:length(element)
	fseek(fp,23+(element(i)-1)*8*dim(1),'bof');
	xa(:,i) = fread(fp,dim(1),'double');
end

fclose(fp);



function bend = spm_bigend
% Checks to see if the computer is big-endian.
% FORMAT bend = spm_bigend
% bend   - 1 for yes, 0 for no, NaN for don't know and Inf for non IEEE
%          format floats.
%
% I don't know about some of the architectures - so it may be worth
% checking the code.
%_______________________________________________________________________
% %W% John Ashburner %E%

computers = str2mat('PCWIN','MAC2','SUN4','SOL2','HP700','SGI',...
	'SGI64','IBM_RS','ALPHA','AXP_VMSG','AXP_VMSIEEE','LNX86',...
	'VAX_VMSG','VAX_VMSD');
endians = [0 1 1 1 1 1 1 1 0 Inf 0 0 Inf Inf];
c=computer;
bend = NaN;
for i=1:size(computers,1),
	if strcmp(c,deblank(computers(i,:))),
		bend = endians(i);
		break;
	end;
end;
if ~finite(bend),
	if isnan(bend),  
		error(['I don''t know if "' c '" is big-endian.']);
	else,
		error(['I don''t think that "' c '" uses IEEE floating point ops.']);   
	end;
end;
return;

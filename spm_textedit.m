% facility to browse ASCII files in SPM
%----------------------------------------------------------------------------
%
%__________________________________________________________________________
% %W% %E%


[s,P] = unix('ls spm_*.man');
u     = findstr(P,'spm_')
for i = 1:length(u)
	d   = P(u(i):length(P));
	d   = d(1:(min(find(d == 10)) - 1));
	unix(['textedit ' d]);
end

[s,P] = unix('ls spm_*.m');
u     = findstr(P,'spm_')
for i = 1:length(u)
	d   = P(u(i):length(P));
	d   = d(1:(min(find(d == 10)) - 1));
	unix(['textedit ' d]);
end

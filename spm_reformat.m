% reformat IOP data
%---------------------------------------------------------------------------
%
%__________________________________________________________________________
% %W% %E%

set(2,'Name','fMRI reformat for IOP data')
P      = spm_get(Inf,'.img','select files to reformat');

[DIM VOX SCALE TYPE OFFSET ORIGIN DESCRIP] = spm_hread(P(1,:));


N      = DIM(3);                              % number in time-series
DIM    = [DIM(1:2) size(P,1)];                % dimnesions of data
OFFSET = 0;
SCALE  = 1;

CWD    = P(1,:);
CWD    = CWD(1:max(find(CWD == '/')));

% set parameters
%---------------------------------------------------------------------------
p      = DIM(1)*DIM(2);


% filenames - output
%---------------------------------------------------------------------------
Q     = [];
for i = 1:N
	Q = [Q; [CWD sprintf('fMRI%3.0i.img',i)]]; end
[m n] = size(Q);
Q     = Q(:);
Q(Q == ' ') = '0'*ones(1,sum(Q == ' '));
Q     = setstr(reshape(Q,m,n));


% create files and headers
%---------------------------------------------------------------------------
for i = 1:size(Q,1)
	fo = fopen(Q(i,:),'w');
	fclose(fo);
	spm_hwrite(Q(i,:),DIM,VOX,SCALE,TYPE,OFFSET);
end
 
% write files
%---------------------------------------------------------------------------
for i = 1:size(Q,1)
	for j = 1:DIM(3)
		fi = fopen(P(j,:),'r');
		fo = fopen(Q(i,:),'a+');
		fseek(fi,(i - 1)*p*2,-1);
		fwrite(fo,fread(fi,p,spm_type(TYPE)),spm_type(TYPE));
		fclose(fo);
		fclose(fi);
	end
end

%---------------------------------------------------------------------------
set(2,'Name',' ')

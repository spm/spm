function spm_image_funks(P,Q,func)
% Perform algebraic functions on images.
% FORMAT spm_image_funks(P,Q,func) OR spm_wi
% P    - matrix of input image filenames.
% Q    - name of output image.
% func - the expression to be evaluated.
%
% Without any arguments, spm_wi acts as a spm_wi_ui.
%_______________________________________________________________________
% The images specified in P, are referred to as i1, i2, i3... in the 
% expression to be evaluated.
% If images are different sizes and orientations, then the size and
% orientation of the first is used for the output image.
% The image Q is written to the same directory as the first input image.
%_______________________________________________________________________
% %W% John Ashburner FIL %E%

if (nargin==0)
	spm_figure('Clear','Interactive');
	P = spm_get(Inf,'.img','Images to work on');
	Q = spm_input('Output filename',1,'s');
	func = spm_input('Evaluated Function',2,'s');
	set(spm_figure('FindWin','Interactive'),'Name','Computing..','Pointer','Watch');
	drawnow;
	spm_image_funks(P,Q,func);
	spm_figure('Clear','Interactive');
	return;
end



Hold = 1;
m = size(P,1);

q = max([find(Q == '/') 0]);
Q = [spm_str_manip(Q((q+1):length(Q)),'sd') '.img'];
p = P(1,:);
p = p(p ~= ' ');
q = max([find(p == '/') 0]);
Q = [p(1:q) Q];

if (nargin < 3) func = 'nothing'; end

p = spm_str_manip(P(1,:), 'd');
M = spm_get_space(p);
[DIM VOX SCALE TYPE OFFSET ORIGIN] = spm_hread(p);

% Get properties of all the images
%----------------------------------------------------------------------------
Matrixes = zeros(16,m);
V        = zeros(12,m);
for i = 1:m,
	p  = deblank(P(i,:));
	C2 = spm_get_space(p);
	M1  = inv(M)*C2;
	Matrixes(:,i) = M1(:);
	V(:,i) = spm_map(p);
end


Output = zeros(prod(DIM(1:2)),DIM(3));

X = kron(ones(DIM(2),1),[1:DIM(1)]');
Y = kron([1:DIM(2)]',ones(DIM(1),1));

% Start progress plot
%----------------------------------------------------------------------------
spm_progress_bar('Init',DIM(3),func,'planes completed');

for j = 1:DIM(3),
	B     = spm_matrix([0 0 -j 0 0 0 1 1 1]);

	for i = 1:m,
		M0  = reshape(Matrixes(:,i),4,4);
		M1 = inv(B*M0);

		d  = spm_slice_vol(V(:,i),M1,DIM(1:2),Hold);
		eval(['i' num2str(i) '=d;']);
	end

	eval(['d1 = ' func ';'], ['error([''Cant evaluate "'' func ''".'']);']);
	if (prod(size(d)) ~= prod(size(d1))) error(['"' func '" produced incompatible image.']); end
	Output(:,j) = d1(:);

	spm_progress_bar('Set',j);
end

% Write output image (16 bit signed)
%------------------------------------------------------------------
mx = max(max(Output));
SCALE  = mx/32767;
if SCALE==0, SCALE=1; end; %For images of all zeros
fp = fopen(Q,'w');
for j = 1:DIM(3),
	d = round(Output(:,j)/SCALE);
	tmp = find(d > 32767);
	d(tmp) = zeros(size(tmp))+32767;
	tmp = find(d < -32768);
	d(tmp) = zeros(size(tmp))-32768;
	fwrite(fp,d,spm_type(4));
end
fclose(fp);
spm_hwrite(Q,DIM,VOX,SCALE,4,0,ORIGIN,'spm - algebra');
spm_get_space(Q,M);

for i = 1:m, spm_unmap(V(:,i)); end
spm_progress_bar('Clear');

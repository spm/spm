function P = spm_imatrix(M)
% returns the parameters for creating an affine transformation
% FORMAT P = spm_imatrix(M)
% M      - Affine transformation matrix
% P      - Parameters (see spm_matrix for definitions)
%___________________________________________________________________________
% %W% John Ashburner & Stefan Kiebel %E%

% Translations and zooms
%-----------------------------------------------------------------------
R         = M(1:3,1:3);
C         = chol(R'*R);
Pt        = [M(1:3,4)' 0 0 0  diag(C)'  0 0 0];
if det(R)<0, Pt(7)=-Pt(7);end % Fix for -ve determinants

% Shears
%-----------------------------------------------------------------------
C         = diag(diag(C))\C;
Pt(10:12) = C([4 7 8]);
R0        = spm_matrix([0 0 0  0 0 0 Pt(7:12)]);
R0        = R0(1:3,1:3);

% This just leaves us with rotations in matrix R1
%-----------------------------------------------------------------------
R1        = R/R0;

% We aren't sure of `x' if we know `sin(x)', so try different
% combinations.
%-----------------------------------------------------------------------
perm = [
 1  1  1
 1  1  0
 1  0  1
 1  0  0
 0  1  1
 0  1  0
 0  0  1
 0  0  0]';

ss = Inf;
for p=perm

	t = R1(1,3);
	t = min(max(t, -1), 1); % Fix for occasional rounding errors
	Pt(5) = asin(t);
	if p(1)==1, Pt(5) = pi-Pt(5); end;

	t = R1(1,2)/cos(Pt(5));
	t = min(max(t, -1), 1); % Fix for occasional rounding errors
	Pt(6) = asin(t);
	if p(2)==1, Pt(6) = pi-Pt(6); end;

	t = R1(2,3)/cos(Pt(5));
	t = min(max(t, -1), 1); % Fix for occasional rounding errors
	Pt(4) = asin(t);
	if p(3)==1, Pt(4) = pi-Pt(4); end;

	% See which version works best
	%-----------------------------------------------------------------------
	M1 =	[1 0 0; 0 cos(Pt(4)) sin(Pt(4)); 0 -sin(Pt(4)) cos(Pt(4))] * ...
		[cos(Pt(5)) 0 sin(Pt(5)); 0 1 0; -sin(Pt(5)) 0 cos(Pt(5))] * ...
		[cos(Pt(6)) sin(Pt(6)) 0; -sin(Pt(6)) cos(Pt(6)) 0; 0 0 1];
	s = sum((M1(:)-R1(:)).^2);
	if (s<ss), ss = s; P = Pt; end
end

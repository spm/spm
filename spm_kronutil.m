% A mex routine for doing things with Kronecker Tensor Products
% FORMAT [alpha,beta]=spm_kronutil(img1,img2,b1,b2)
%
% Performs:
% 
% 	[m1,m2]=size(img1);
% 	[m1,m2]=size(img2);
% 	[m1,n1]=size(b1);
% 	[m2,n2]=size(b2);
% 
% 	alpha = zeros(n1*n2,n1*n2);
% 	beta  = zeros(n1*n2,1);
% 	for i=1:m2
% 		tmp   = kron(img1(:,i),ones(1,n1)).*b1;
% 		alpha = alpha + kron(b2(i,:)'*b1(i,:),  tmp'*tmp);
% 		beta  = beta  + kron(b2(i,:)', tmp'*img2(:,i));
% 	end
% 
% which is equivalent to, but a lot faster than:
% 
% 	B     = kron(b2,b1);
% 	A     = diag(img1(:))*B;
% 	b     = img2(:);
% 	alpha = A'*A;
% 	beta  = A'*b;
%_______________________________________________________________________
% %W% John Ashburner %E%

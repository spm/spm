
function [X, pca_order] = pinv_plus(A,varargin)
% PINV_PLUS   Modified Pseudoinverse which allows a specified
% dimensionality (order) to be passed in.
%
%   pinv_plus(A) is the same as pinv(A)
%
%   [X, pca_order]=pinv_plus(A,pca_order) uses at most only the top pca_order PCs to compute inverse
% 
%   pinv_plus(A,pca_order,1) uses precisely the top pca_order PCs to compute inverse
%
%   If pca_order=-1 then calls spm_pca_order to compute order to use
%
%   Class support for input A: 
%      float: double, single
%
%   Mark Woolrich

tol=[];
r=size(A,1);

if nargin >= 2,
   pca_order = varargin{1};
   if(pca_order==-1),
      pca_order=spm_pca_order(A);
   end;
else
   pca_order = size(A,1);
end;

max_r=pca_order;

if nargin >= 3,
   fix_to_max_r=varargin{2};
else
   fix_to_max_r=0;
end;

if isempty(A)     % quick return
  X = zeros(size(A'),class(A));  
  return  
end

% [U,S,V] = svd(A,0); s = diag(S);
[m,n] = size(A);

if n > m
   X = pinv(A',varargin{:})';
else
   if(max_r<size(A,1))
    %[U,S,V] = svds(A,max_r);
    [U,S,V] = svd(A,0);
    max_r=min(length(diag(S)),max_r);
   else
    [U,S,V] = svd(A,0);
   end;
    

   if m > 1, s = diag(S);
      elseif m == 1, s = S(1);
      else s = 0;
   end
   
   tol=[];
   
   if max_r < 0        

     diffs=abs(diff(s)./s(1:end-1));

     for i=2:length(diffs),
         if(diffs(i) > 6*mean(diffs(1:i))), 
             break;
         end;
     end; 
     
     max_r = i;

     %figure;plot(diffs);ho;plot(r,diffs(r),'*');

   end;

%   max_r
%   fix_to_max_r
   
   if(fix_to_max_r==1),      
       r=max_r;
   else,
       tol=max(m,n) * eps(max(s));

       %tol=max(size(A)) * norm(A) * eps(class(A));

       r = sum(s > tol);                 
       
       r = min(r,max_r);
   end;            
   
   %figure;plot(s);
   
   s2=s;
   
   if (r == 0)
      X = zeros(size(A'),class(A));
   else
      s = diag(ones(r,1)./s(1:r));
      X = V(:,1:r)*s*U(:,1:r)';
   end
   
   pca_order=r;
   
end

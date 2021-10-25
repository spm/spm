classdef ndSparse
%ndSparse - A class of N-dimensional sparse arrays.
%
% by Matt Jacobson
%
% Copyright, Xoran Technologies, Inc. 2010
%
%
% USAGE: 
%   
%   S=ndSparse(X) where X is an ordinary MATLAB sparse matrix converts X into 
%     an ndSparse object. S can be reshaped into an N-dimensional sparse array using
%     its RESHAPE method, for arbitrary N. 
%   
%   S=ndSparse(X,[M,N,P,...]) is equivalent to reshape(ndSparse(X),[M,N,P,...]). 
%   
%The class also has a variety of static methods that can be used to construct instances
%of the class. For example,
%
%        S=ndSparse.build(Coordinates,Values,[m,n,p,...],nzmax) 
%
%lets you generate an N-dimensional sparse array from a table of explicit entries. 
%This is a generalization to N dimensions of S=sparse(i,j,s,m,n,nzmax).
%
%Other such methods include:
%      ndSparse.accumarray 
%      ndSparse.sprand 
%      ndSparse.sprandn 
%      ndSparse.spalloc 
%   
%EXAMPLES: 
%   
%   >> A=ndSparse.build( [1 1 1; 2 1 1;2 2 2] , [50,60 70]) %Builds a 2x2x2 sparse  
%                                                           %array from table
%
%              A(:,:,1) =
%
%                 (1,1) 50 
%                 (2,1) 60 
% 
%              A(:,:,2) =
%
%                 (2,2) 70
%
%Many of the same manipulations common to ordinary multidimensional MATLAB full arrays
%can be performed on the sparse 3D array A generated above. It can be permuted, summed, 
%concatentated, and so forth e.g.,
%   
%   >> B=sum( permute([A,A+10], [3,2,1]) ,2)
%
%              B(:,:,1) =
%
%                (1,1) 120 
%                (2,1) 20 
% 
%              B(:,:,2) =
%
%                (1,1) 140 
%                (2,1) 160
%
%Other overloaded methods include BSXFUN, REPMAT, CIRCSHIFT, CONVN, FLIPDIMS, SQUEEZE,
%SHIFTDIM and many more. Type "methods ndSparse" for a full list and use 
%"help ndSparse.methodname"  to get details of usage. 
%
%When browsing the list of methods, note that certain common operations have 
%different implementations, optimized for different situations. Specifically,
%SUM, ANY,ALL, MIN, MAX... have alternative implementations SUMML, ANYML, ALLML,
%MINML, MAXML which are optimized for "low-dimensional" ndSparse objects OBJ.
%Here, low-dimensional means that a normal N-column MATLAB sparse matrix won't
%consume too much memory on your platform for N=MAX(NUMEL(OBJ)./SIZE(OBJ)).
%
%
%Another feature of the class is that bi-operand operations are allowed between ndSparse
%objects and MATLAB objects of any numeric type (single, uint16, etc...). This is not true
%of ordinary MATLAB sparse matrices, as of R2010b. 
%   
%   >> C=eye(2,'single')*B(:,:,2)
%
%      C =
%
%        (1,1) 140 
%        (2,1) 160 
%
%   >> whos A B C
%
%   Name Size Bytes Class Attributes
%
%   A 2x2x2 136 ndSparse 
%   B 2x1x2 140 ndSparse 
%   C 2x1 104 ndSparse 
%   
%To convert back to an ordinary n-D full array, use the class' overloaded FULL method.
%To convert to a normal 2D sparse matrix, use the methods SPARSE or SPARSE2D. For example,
%SPARSE2D will convert an MxNxPx...xQ ndSparse array to the two dimensional 
%(M*N*P*...)xQ sparse matrix in native MATLAB form. 


 properties (Constant)
    
     oldstyle=false; %Set this to true so that sum, mean, any, all etc...
                     %use memory-liberal algorithms
     
 end

 properties (Access=private)
        data;
        ndShape;       
 end
 
 
 
 methods
    function obj=ndSparse(data,ndShape)
    %ndSparse class constructor
    %
    % OBJ=ndSparse(A) where A is an ordinary MATLAB sparse matrix converts A into 
    % an ndSparse object. OBJ can be reshaped into an N-dimensional sparse array 
    % using its RESHAPE method, for arbitrary N. 
    %   
    % OBJ=ndSparse(A,[M,N,P,...]) is equivalent to reshape(ndSparse(A),[M,N,P,...]) 


         if nargin<2, 
             ndShape=size(data); 
         elseif prod(ndShape)~=numel(data)
             error 'The number of DATA elements is not consistent with the given ndShape.'
         end
         
         [mm,nn,ndShape]=ndShape2MN(ndShape);
         
         if isa(data,'ndSparse')
            obj=data; 
            if nargin>1 
              obj.data=reshape(obj.data,mm,nn);   
              obj.ndShape=ndShape;
            end
            return 
            
         end
        
         
         
         if ~issparse(data), 
             data=reshape(data,mm,nn);
             data=sparse(mkCompat(data));
         end    
        
         data=reshape(data,mm,nn);
         
         obj.data=data;
         obj.ndShape=ndShape;

     end
   
     function obj=set.ndShape(obj,ndShape)
         
         obj.ndShape=untrail1s(ndShape(:).');
         
     end
     
     function varargout=size(obj,varargin)
     %size - same syntax as for full arrays    
     
         varargout=parseSize(obj.ndShape,nargout,varargin{:});
         
     end

     function L=length(obj)
     %length - same behavior as for full arrays    
     
       if isempty(obj)
          L=0;
       else
          L=max(obj.ndShape);
       end
         
     end
     
     function out=end(obj,K,N)
     %end - sameindexing rules as for full arrays
     
         sz= obj;
         if N==1, 
             out=numel(sz); return
         end
         
         [c{1:N}]=size(obj);
         c=[c{:}];
         out=c(K);
         
     end
        
     
     function N=ndims(obj)
     %ndims - same behavior as for full arrays
     
         N=length(size(obj));
         
     end
     
     function N=nnz(obj)
     %nnz - same behavior as for full arrays    
        N=nnz(obj.data); 
     end
  
     function N=nzmax(obj)
     %nzmax - same behavior as for ordinary 2D sparse arrays    
        N=nzmax(obj.data); 
     end
     
     function s=nonzeros(obj)
     %nnonzers - same behavior as for ordinary 2D sparse arrays    
        s=nonzeros(obj.data); 
     end     
     
     
     function N=numel(obj)
     %numel - same behavior as for full arrays    
        N=numel(obj.data); 
     end
     
     function bool=isempty(obj)
      %isempty - same behavior as for full arrays    
        bool=isempty(obj.data); 
     end
     
     function bool=islogical(obj)
      %islogical - returns true of data is stored as type logical   
        bool=islogical(obj.data); 
     end
     
     function bool=issparse(obj)
      %issparse - returns true for ndSparse objects
        bool=true; 
     end
     
     function bool=isnumeric(obj)
      %isnumeric - returns true for non-logical ndSparse objects
        bool=~islogical(obj); 
     end
     
     function bool=isfloat(obj)
     %isfloat - returns true for non-logical ndSparse objects
        bool=~islogical(obj); 
     end 
     
     function bool=isreal(obj)
     %isreal - works as for normal arrays
        bool=isreal(obj.data); 
     end 
     
     
     function bool=isequal(varargin)
     %isequal - works as for normal arrays
        
      N=length(varargin);
      sizeCell=varargin;
      bool=false;
      
           for ii=1:N%for-loop instead of cellfun is deliberate

             if ~isnumeric(varargin{ii})  return; end
             sizeCell{ii}=size(varargin{ii});

           end

           if ~isequal(sizeCell{:})  return; end

           sz=sizeCell{1};
           
           for ii=1:N%for-loop instead of cellfun is deliberate

             varargin{ii}=mkCompat(varargin{ii},sz);

           end
 
       bool=isequal(varargin{:});
       
     end 
 
     function bool=isequalwithequalnans(varargin)
     %isequalwithequalnans - works as for normal arrays
        
      N=length(varargin);
      sizeCell=varargin;
      bool=false;
      
           for ii=1:N%for-loop instead of cellfun is deliberate

             if ~isnumeric(varargin{ii})  return; end
             sizeCell{ii}=size(varargin{ii});

           end

           if ~isequal(sizeCell{:})  return; end

           sz=sizeCell{1};
           
           for ii=1:N%for-loop instead of cellfun is deliberate

             varargin{ii}=mkCompat(varargin{ii},sz);

           end
           
       bool=isequalwithequalnans(varargin{:});
       
     end      
     
     
     function obj=isnan(obj)
     %isnan - works as for normal arrays
        obj.data=isnan(obj.data); 
     end    
     
     function obj=isinf(obj)
     %isinf - works as for normal arrays
        obj.data=isinf(obj.data); 
     end     
     
     function out=isfinite(obj)
     %isinf - works as for normal arrays
        out=reshape( isfinite(obj.data) , obj.ndShape); 
     end         
     
     
     
     
     function obj=real(obj)
     %real - works as for normal arrays
        obj.data=real(obj.data); 
     end    
     
     function obj=imag(obj)
     %imag - works as for normal arrays
        obj.data=imag(obj.data); 
     end        
     
     function obj=conj(obj)
     %conj - works as for normal arrays
        obj.data=conj(obj.data); 
     end        
     
     function obj=abs(obj)
     %abs - works as for normal arrays
        obj.data=abs(obj.data); 
     end      
     
     function obj=sqrt(obj)
     %sqrt - works as for normal arrays
        obj.data=sqrt(obj.data); 
     end        
     
     function c=underlyingClass(obj)
     %underlyingClass(obj) will return the storage type of the data (double/logical)
     
         c=class(obj.data);
     end
     
     function out=logical(obj)
     %logical - converts storage type to double    
         out=logical(obj.data);
     end
     
     function varargout=find(obj,varargin)
     %find - same syntax as for full arrays    
         
        [varargout{1:nargout}]=find(obj.data,varargin{:}); 
        
        if nargout>1,
           II=varargout{1}; 
           JJ=varargout{2};
           
           ndSubs=IJ2ndCoords(II,JJ,obj.ndShape);
           
           II=ndSubs{1};
           
           if length(ndSubs)>2
            JJ=sub2ind(obj.ndShape(2:end),ndSubs{2:end});
           else
            JJ=ndSubs{2};
           end
           
           varargout(1:2)={II,JJ};
           
        end
        
     end
     
     function out=double(obj)
     %double - converts storage type to double
     
         out=double(obj.data);
     end
     
     function out=sparse(obj)
     %sparse(obj) - converts ndSparse object of dimension MxNxPxQx....
     %into an ordinary 2D sparse matrix of dimensions Mx(N*P*Q...)   
         ndShape=obj.ndShape;
         out=reshape(obj.data,ndShape(1),prod(ndShape(2:end)));
     end
     
     function data=sparse2d(obj)
     %sparse2d(obj) - converts ndSparse object of dimension MxNxPxQx...YxZ
     %into an ordinary 2D sparse matrix of dimensions (M*N*P*Q*...*Y)xZ
     
         data=sparse(obj.data);
     end
     

     
     
     function data=full(obj)
     %full - same syntax as for 2D sparse matrices
         
        data=reshape( full(obj.data), obj.ndShape);
        
     end
     
     function obj=spfun(fun,obj)
     %spfun - same syntax as for 2D sparse matrices
     
        obj.data=spfun(fun,obj.data);
        
     end
     
     function obj=spones(obj)
     %spones - same syntax as for 2D sparse matrices
     
        obj.data=spones(obj.data);
        
     end
     
     function obj=reshape(obj,varargin)
     %reshape - same syntax as for full arrays
     
         map=cellfun('isempty',varargin);
         switch nnz(map)
            
             case 0
                 
             case 1
                 missing=prod(obj.ndShape)/prod([varargin{:}]);
                 if mod(missing,1), error 'Bad reshape args'; end
                 varargin(map)={missing};
             otherwise
                 
                 error 'Too many empty args'
             
         end
         
         newshape=[varargin{:}];
         
         obj=ndSparse(obj,newshape);
  
         
     end
     
    function objnew=subsref(obj,S)
    %subsref - same indexing rules as for full arrays
    
        if isempty(S.subs)
            objnew=obj; return
        end
    

        
        [E,targetshape,~,other]=nd2matrixIndex(S.subs,obj.ndShape,'subsref');

        
        if  other.boolConsol 
              
             obj=reshape(obj,other.quasiShape);  
             
        end
        
        
        
          data=builtin('subsref',obj.data,E);
              
          objnew=ndSparse(data,targetshape);
  
    end

   
    
    function objnew=subsasgn(obj,S,rhs)
    %subsasgn - same indexing rules as for full arrays
      
     
        if isempty(S.subs)
            error ' An indexing expression on the left side of an assignment must have at least one subscript.'
        end
        
    

        
        
        if isequal(rhs,[])%NULL ASSIGNMENT
            
          sz=size(obj);
          nn=length(sz);            
          N=length(S.subs);        

            dim=find(~strcmp(':',S.subs));
            qq=length(dim);
            if qq>1, 
                error  'A null assignment can have only one non-colon index.';
            elseif qq==0
                dim=1;
                idx=':';
            else
                idx=S.subs{dim};
            end
           
            
            

              
        
               
                if 1<N && N<nn  %Consolidated indexing of trailing dims
                 
                  outcell=parseSize(size(obj),N);  
                  sz=[outcell{:}];
                  obj=reshape(obj,sz);  
                  nn=N; 
                  
                end
             
   
              if dim>nn
                 T0=1; 
              else
                  T0=1:sz(dim);
              end
              
            try  
            idx=T0(idx); %force idx to be numeric nonlogical
            catch ME
                
                strpres=@(t,p) ~isempty(strfind(t,p));
                if strpres(ME.identifier, 'badsubscript')
                
                   if strpres(ME.message,'index must be a positive integer or logical')
                    error 'Subscript indices must either be real positive integers or logicals.'
                   elseif strpres(ME.message,'index out of bounds')
                       error 'Index exceeds matrix dimensions.'
                   end
                else
                    
                    error 'Bad subscript'
                    
                end
            end
             
            if dim>nn,
                
                sz=trail1s(sz,dim);
                sz(end)=0;
                objnew=ndSparse.build(sz);
                return; 
                
            end
            
            
           if ndSparse.oldstyle%memory liberal style  
          
                data=ExtractDim1Reshaping(obj.data,dim,sz);

                 data(idx,:)=[];
                 sz(dim)=size(data,1);

                data=invExtractDim1Reshaping(data,dim,sz);

                objnew=ndSparse(data,sz);
            
           else
              
                if isequal(idx,':') || ( sz(dim)==1 && ~isempty(idx) )
               
                    sz(dim)=0;
                    objnew=ndSparse.build(sz);
                    
                else
                    
                   
                  
                  
                  
                  [ndSubs,vals]=getEntryTable(obj,length(sz));
                  
                   
                  
                  todelete=ismember(ndSubs{dim},idx);
                  
                  ndSubs=[ndSubs{:}];
                  ndSubs(todelete,:)=[];
                  vals(todelete)=[];
                  
                  T=T0;
                   T(idx)=[];
                                
                  
                  T0(T)=1:length(T);
                   ndSubs(:,dim)=T0(ndSubs(:,dim));
                  

                  
                  sz(dim)=length(T);
                  objnew=ndSparse.build(ndSubs , vals, sz);  
                    
                 
                end
               
           end
         
        elseif isempty(rhs)   
          
            error 'Improper null assignment. Right hand side must be [].'
            
        else%NON-NULL,
            
         [E,targetshape,boolGrow,other]=nd2matrixIndex(S.subs,obj.ndShape,'subsasgn');
        
         
         
          if boolGrow
              
              N = length(targetshape); 
              
              [ndSubs,vals]=getEntryTable(obj,N);
              
              objnew=ndSparse.build(ndSubs,vals,targetshape,nzmax(obj));
 
          elseif  other.boolConsol  
              
            
             objnew=reshape(obj,other.quasiShape);  
             
             
          else
              
             objnew=obj;
              
          end
          
          
          
          
          data=objnew.data;

          if ~isscalar(rhs)  
              
               if prod(other.subshapeND)~=numel(rhs)
                 errorFlag=true;
               else
                  errorFlag=false;
               end
             
              
            if other.LinearIndexing                
                if errorFlag, error ' In an assignment  A(:) = B, the number of elements in A and B must be the same'; end
            else
               if errorFlag
                   error 'Subscripted assignment dimension mismatch.'
               end
               rhs=reshape(rhs,other.subshape);
            end
          end
            
          data=builtin('subsasgn',data,E,rhs);
          objnew=ndSparse(data , targetshape);
          
          
          
        end
        
              
    end

    function idx=subsindex(obj)
    %subsindex - same indexing rules as for full arrays
        
       switch class(obj.data);
           
           case 'double'
             idx=obj.data-1;
           case 'logical'
             idx=find(obj.data)-1; 
       end
    end 
    

    function objnew=sum(obj,varargin)
    %sum - same syntax as for full arrays    
        
      objnew=sumEngine(0,obj,varargin{:});
        
    end
     
    function objnew=summl(obj,varargin)
    %summl - An alternative implementation of SUM. It
    %makes more liberal use of memory allocation, but is
    %sometimes faster. Best suited for low dimensional arrays.
        
      objnew=sumEngine(1,obj,varargin{:});
        
    end       
       
              
    function objnew=any(obj,varargin)
    %any - same syntax as for full arrays    
        
      objnew=anyEngine(0,obj,varargin{:});
        
    end
     
    function objnew=anyml(obj,varargin)
    %anyml - An alternative implementation of ANY. It
    %makes more liberal use of memory allocation, but is
    %sometimes faster. Best suited for low dimensional arrays.
        
      objnew=anyEngine(1,obj,varargin{:});
        
    end       
              
    function objnew=all(obj,varargin)
    %all - same syntax as for full arrays    
        
      objnew=allEngine(0,obj,varargin{:});
        
    end
     
    function objnew=allml(obj,varargin)
    %allml - An alternative implementation of ALL. It
    %makes more liberal use of memory allocation, but is
    %sometimes faster. Best suited for low dimensional arrays. 
        
      objnew=allEngine(1,obj,varargin{:});
        
    end       
    
    
    function objnew=mean(obj,varargin)
    %mean - same syntax as for full arrays    
        
      objnew=meanEngine(0,obj,varargin{:});
        
    end
     
    function objnew=meanml(obj,varargin)
    %meanml - An alternative implementation of MEAN. It
    %makes more liberal use of memory allocation, but is
    %sometimes faster. Best suited for low dimensional arrays.    
        
      objnew=meanEngine(1,obj,varargin{:});
        
    end  
       
    function objnew=circshift(obj,varargin)
    %circshift - same syntax as for full arrays    
        
      objnew=circshiftEngine(0,obj,varargin{:});
        
    end
     
    function objnew=circshiftml(obj,varargin)
    %circshiftml - An alternative implementation of CIRCSHIFT. It
    %makes more liberal use of memory allocation, but is
    %sometimes faster. Best suited for low dimensional arrays.   
        
      objnew=circshiftEngine(1,obj,varargin{:});
        
    end      
       
    function objnew=cat(obj,varargin)
    %cat - same syntax as for full arrays    
        
      objnew=catEngine(0,obj,varargin{:});
        
    end
     
    function objnew=catml(obj,varargin)
    %catml - An alternative implementation of CAT. It
    %makes more liberal use of memory allocation, but is
    %sometimes faster. Best suited for low dimensional arrays.   
        
      objnew=catEngine(1,obj,varargin{:});
        
    end          
       

    function varargout=max(obj,varargin) %max
    %max - same syntax as for full arrays    
        
      if numel(varargin)==1
          varargout{1}=bsxfun(@max,obj,varargin{1});
      else
    
      [varargout{1:nargout}]=maxEngine(0,obj,varargin{:});
        
      end
      
      
    end
     
    function varargout=maxml(obj,varargin)%maxml
    %maxml - An alternative implementation of MAX. It
    %makes more liberal use of memory allocation, but is
    %sometimes faster. Best suited for low dimensional arrays.   
        
      [varargout{1:nargout}]=maxEngine(1,obj,varargin{:});
        
    end       
       
    function varargout=min(obj,varargin)%min
    %min - same syntax as for full arrays    
        
      if numel(varargin)==1
          varargout{1}=bsxfun(@min,obj,varargin{1});
      else
    
      [varargout{1:nargout}]=minEngine(0,obj,varargin{:});
        
      end
    end
     
    function varargout=minml(obj,varargin)%minml
    %minml - An alternative implementation of MIN. It
    %makes more liberal use of memory allocation, but is
    %sometimes faster. Best suited for low dimensional arrays.   
        
      [varargout{1:nargout}]=minEngine(1,obj,varargin{:});
        
    end       
       
       function obj=horzcat(varargin)
       %horzcat - same syntax as for full arrays    
           
           obj=cat(2,varargin{:});
       end
       
       
       function obj=vertcat(varargin)
       %vertcat - same syntax as for full arrays
       
           obj=cat(1,varargin{:});
       end
       

       

         
         function objnew=permute(obj,order)
         %permute - same syntax as for full arrays

           
           ndShape=obj.ndShape;
           N=length(ndShape);
           
           if length(order)<N, 
               error 'ORDER must have at least N elements for an N-D array'
           else
               [ndShape,N]=trail1s(ndShape,max(order));
           end
           
                        
            if issorted(order), objnew=obj; return; end
             
           [ndSubs,vals]=getEntryTable(obj,N);
     
           ndSubs=ndSubs(:,order);
           newshape=ndShape(order);
           

           objnew=ndSparse.build(ndSubs,vals,newshape,nzmax(obj));
           
           
       end
       
  
       function objnew=ipermute(obj,order) 
       %ipermute  - same syntax as for full arrays    
           
           [null,order]=sort(order);
           objnew=permute(obj,order);
           
       end
 
       
        %%Unary ops    
        function  obj=uplus(obj)
        %uplus - same behavior as for full arrays
            
        end

        function  obj=uminus(obj)
        %uminus - same behavior as for full arrays
        
           obj.data=-obj.data;
            
        end
        
        
        function  obj=not(obj)
        %not - same behavior as for full arrays  
            obj.data = ~obj.data ;
            
        end
        
        %%Binary ops
        
        function  obj=plus(L,R)
        %plus - same behavior as for full arrays  
            
        
              obj=bsxfun(@plus,L,R);
        
%            [s,L,R]=getshape(L,R);  
%           obj = finalObject( mkCompat(L,s) + mkCompat(R,s),s);
            
        end

        function  obj=minus(L,R)
        %minus - same behavior as for full arrays
        
        
        
               obj=bsxfun(@minus,L,R);       
        
%            [s,L,R]=getshape(L,R);
%            obj = finalObject( mkCompat(L,s) - mkCompat(R,s),s);
            
        end        
        
        
        
       function  obj=times(L,R)
       %times - same behavior as for full arrays
       
       
              obj=bsxfun(@times,L,R);     
       
%             [s,L,R]=getshape(L,R);
%             obj = finalObject( mkCompat(L,s) .* mkCompat(R,s),s);
            
       end

        
       function  obj=rdivide(L,R)
       %rdivide - same behavior as for full arrays
       
       
               obj=bsxfun(@rdivide,L,R);   
       
%                [s,L,R]=getshape(L,R);
%                obj = finalObject( mkCompat(L,s) ./ mkCompat(R,s), s);
            
        end

        function  obj=ldivide(L,R)
        %ldivide - same behavior as for full arrays
        
        
              obj=bsxfun(@ldivide,L,R); 
        
%              [s,L,R]=getshape(L,R);
%              obj = finalObject( mkCompat(L,s) .\ mkCompat(R,s),s);
            
        end        
        



        function  obj=power(L,R)
        %power - same behavior as for full arrays
        
        
        
           obj=bsxfun(@power,L,R); 
        
%              [s,L,R]=getshape(L,R);
%              
%              obj = finalObject( mkCompat(L,s) .^ mkCompat(R,s),s);
            
        end     

        

        function  obj=lt(L,R)
        %lt - same behavior as for full arrays

        
             obj=bsxfun(@lt,L,R);    
        
%             [s,L,R]=getshape(L,R);
%             obj = finalObject( mkCompat(L,s) < mkCompat(R,s),s);

        end   
        
        function  obj=le(L,R)
        %le - same behavior as for full arrays    

        
              obj=bsxfun(@le,L,R);               
        
%              [s,L,R]=getshape(L,R);
%              obj = finalObject( mkCompat(L,s) <= mkCompat(R,s),s);

        end         
        
        function  obj=gt(L,R)
        %gt - same behavior as for full arrays    

                obj=bsxfun(@gt,L,R);         
        
%              [s,L,R]=getshape(L,R);
%              obj = finalObject( mkCompat(L,s) > mkCompat(R,s),s);

        end      
        
  
        function  obj=ge(L,R)
        %ge - same behavior as for full arrays    
 
        
                obj=bsxfun(@ge,L,R);         

%              [s,L,R]=getshape(L,R);
%              obj = finalObject( mkCompat(L,s) >= mkCompat(R,s),s);

        end   
     
        function  obj=eq(L,R)
        %eq - same behavior as for full arrays    
            
        
           obj=bsxfun(@eq,L,R);             
        
        
%            [s,L,R]=getshape(L,R);
%            obj = finalObject( mkCompat(L,s) == mkCompat(R,s),s);

        end     
        
        function  obj=ne(L,R)
        %ne - same behavior as for full arrays    

              obj=bsxfun(@ne,L,R);    
        
%            [s,L,R]=getshape(L,R);
%            obj = finalObject( mkCompat(L,s) ~= mkCompat(R,s),s);

        end    

        
        
        function  obj=and(L,R)
        %and - same behavior as for full arrays    

        
             obj=bsxfun(@and,L,R);           
        
%            [s,L,R]=getshape(L,R);
%            obj = finalObject( mkCompat(L,s) & mkCompat(R,s),s);

        end         
        
        function  obj=or(L,R)
        %or - same behavior as for full arrays
        
        
             obj=bsxfun(@or,L,R);                   
        
%              [s,L,R]=getshape(L,R);
%              obj = finalObject( mkCompat(L,s) | mkCompat(R,s),s);

        end    
        
        function obj=xor(L,R)
        %xor - same behavior as for full arrays    
            
          obj=bsxfun(@xor,L,R);
        
        end
        
        
        function obj=rem(L,R)
        %rem - same behavior as for full arrays    
            
          obj=bsxfun(@rem,L,R);
        
        end
        
        function obj=mod(L,R)
        %mod - same behavior as for full arrays    
            
          obj=bsxfun(@mod,L,R);
        
        end       
        
        function obj=atan2(L,R)
        %atan2 - same behavior as for full arrays    
            
          obj=bsxfun(@atan2,L,R);
        
        end    
        
        function obj=atan2d(L,R)
        %atan2d - same behavior as for full arrays    
            
          obj=bsxfun(@atan2d,L,R);
        
        end         
        
        function obj=hypot(L,R)
        %hypot - same behavior as for full arrays    
            
          obj=bsxfun(@hypot,L,R);
        
        end          
        
        
        
       %%Methods for 2D Matrices 
       function  obj=inv(obj)
       %inv - same behavior as for full arrays
       
        if ndims(obj)>2, error 'Operation transpose defined for 2D arrays only'; end 
        
          obj = ndSparse(inv(sparse2d(obj)));
            
       end
       
       
       function  obj=triu(obj,varargin)
       %triu
       
        if ndims(obj)>2, error 'Operation transpose defined for 2D arrays only'; end 
        
          obj = ndSparse(  triu(  sparse2d(obj) , varargin{:} )  );
            
       end
       
        function  obj=tril(obj,varargin)
        %tril - same behavior as for full arrays
       
        if ndims(obj)>2, error 'Operation transpose defined for 2D arrays only'; end 
        
          obj = ndSparse( tril(  sparse2d(obj) , varargin{:} ));
            
       end      
       
       
       function  obj=transpose(obj)
       %transpose - same behavior as for full arrays
       
          if ndims(obj)>2, error 'Operation transpose defined for 2D arrays only'; end  %#ok<*ISMAT>
           
          obj = ndSparse(sparse2d(obj).');
            
       end
       
       function  obj=ctranspose(obj)
       %ctranspose - same behavior as for full arrays
           
          if ndims(obj)>2, error 'Operation ctranspose defined for 2D arrays only'; end  
           
          obj = ndSparse(sparse2d(obj)');
            
       end   
       
       function  obj=mtimes(L,R)
       %mtimes - same behavior as for full arrays
       
          if ( isscalar(R) ||isscalar(L))
              obj=L.*R; return 
          end
           
          if ~all([ndims(L),ndims(R)]==2) , 
              error 'Operation mtimes defined for 2D arrays only'; end  
           
          obj = finalObject( mkCompat(L)*mkCompat(R));
            
        end
       
       function  obj=mrdivide(L,R)
       %mrdivide - same behavior as for full arrays
       
           if ( isscalar(R))
              obj=L./R; return 
          end          
           
          if ~all([ndims(L),ndims(R)]==2), 
              error 'Operation mrdivide defined for 2D arrays only'; end 
          
          obj = finalObject( mkCompat(L)/mkCompat(R));
            
        end

        function  obj=mldivide(L,R)
        %mldivide - same behavior as for full arrays   
            
          if ( isscalar(L))
              obj=L.\R; return 
          end           
            
           if ~all([ndims(L),ndims(R)]==2) && ~isscalar(L), 
               error 'Operation mldivide defined for 2D arrays only'; end 
          
          obj = finalObject( mkCompat(L)\mkCompat(R));
            
        end        
   
       function  obj=mpower(L,R)
       %mpower - same behavior as for full arrays
       
         if ~all([ndims(L),ndims(R)]==2), 
             error 'Operation mpower defined for 2D arrays only'; end   
        
         obj = finalObject( mkCompat(L)^mkCompat(R));
            
       end   
       
    function display(obj)
    %display           
        
        
            l=inputname(1);
            nn=ndims(obj);
            
            simplecase=false;
            if ~prod(obj.ndShape) %display empty array
                
               T=evalc('full(obj)'); %only way to get at the builtin display method
               
               simplecase=true;
               
            elseif nn<=2
                
               T=evalc('sparse2d(obj)'); %only way to get at the builtin display method        
            
               simplecase=true;
                              
            elseif ~nnz(obj)
                
               dimchar=repmat({'-by-'},2,ndims(obj));
               dimchar(1,:)=arrayfun(@(c) num2str(c), size(obj), 'Unif',0);
               dimchar=[dimchar{1:end-1}];
               T=evalc( ' [''   All zero sparse: '' dimchar] ' );
               
               simplecase=true;
               
            end
            
            if simplecase
                
              T=strrep(T,'ans =',[l ' =']);

              jj=find(T~=sprintf('\n'),1,'last');
            
              T=T(1:jj);
           
              disp(T), disp ' '
              return
              
            end
            
            [ndSubs,vals,ndShape]=getEntryTable(obj,ndims(obj));
            
            trailingdims=ndShape(3:end);
            
            data=sparse2d(obj);
            
            
            sz=obj.ndShape(1:2);
            M=prod(sz);
            
            sheets=unique([ndSubs{:,3:end}],'rows');
            sheetCorners=num2cell([ones(size(sheets,1),2), sheets], 1);
            cornerRows=sub2ind(ndShape(1:end-1),sheetCorners{1:end-1});
            cornerCols=sheetCorners{end};
            
          
            rangecell=num2cell(sheets);
            numSheets=size(rangecell,1);

            commas=rangecell(1,:); commas(:)={','};
            rangecell_aschar=cellfun(@(x)num2str(x),rangecell,'uni',0);

              
          
        loopctr=0;    
            
        rg=0:M-1;
            
          for ii=1:numSheets

            loopctr=loopctr+1;  
              
            bookNum=cornerCols(ii);
            rr=cornerRows(ii)+rg;
            
            T=evalc('reshape(data(rr,bookNum),sz)'); %only way to get at the builtin display method

             
            
            if ~isempty(l),
                
                lbl=l;                
                %if numBooks>1
                 cc=[rangecell_aschar(loopctr,:); commas];    
                 lbl=[l '(:,:,' cc{1:end-1} ')'];
                %end
                 T=strrep(T,'ans =',[lbl ' =']);
                 
            end

           
            jj=find(T~=sprintf('\n'),1,'last');
            
            T=T(1:jj);
           
            disp(T), disp ' '
        
            
          end 
     
    end

   
    function objnew=repmat(obj,varargin)
    %repmat - same syntax as for full arrays
     
        repDims=[varargin{:}];
        ndShape=size(obj);
        N=max(length(repDims),length(ndShape));
        
        repDims=trail1s(repDims,N);
        ndShape=trail1s(ndShape,N);

        newshape=ndShape.*repDims;  
        
        if ~all(round(repDims)==repDims)
           error 'Array repititions must be integer-valued' 
        end
        
        if ~all(repDims)%some repDims are 0
           objnew=ndSparse([],newshape); return
        end  

         if ~any(repDims~=1),%No actual repeated data 
             objnew=obj;return
         end

        
         [ndSubs,vals,ndShape]=getEntryTable(obj,N);
         
         if isempty(vals)
             objnew=ndSparse.build(newshape); return 
         end

         entrytable=[ndSubs{:},vals(:)];
         
                for ii=1:N %loop over dims

                   thisRep=repDims(ii);

                   if thisRep==1, continue; end    
                    
                   repCoords=bsxfun(@plus, entrytable(:,ii), (0:thisRep-1)*ndShape(ii) );
                   
                   entrytable=repmat(entrytable,thisRep,1); 
                   entrytable(:,ii)=repCoords(:);

                end
           
           objnew=ndSparse.build(entrytable(:,1:end-1),...
                                 entrytable(:,end),...
                                 newshape);     
                
    end
    
    
    function  out=bsxfun(fun,L,R)%bsxfun
    %bsxfun - same syntax as for full arrays
    
        isL=isa(L,'ndSparse');
        isR=isa(R,'ndSparse');
        isfun=isa(fun,'function_handle');
        
        if ~isfun, error 'First Argument must be a function handle.'; end
        
        
        N=max(ndims(L),ndims(R));
        
        [Lshape{1:N}]=size(L); Lshape=[Lshape{:}];
        [Rshape{1:N}]=size(R); Rshape=[Rshape{:}];          
        
        if any(Lshape~=Rshape & Lshape~=1 & Rshape~=1)
          error 'Non-singleton dimensions of the two input arrays must match each other.'
        end
        
        newshape=max([Lshape;Rshape], [],1);
         newshape=newshape.*logical(Lshape.*Rshape); %mask out the zero-dimensions  
        
        
        Lreps=ones(size(newshape));
        Lreps(Lshape==1)=newshape(Lshape==1);
        
        Rreps=ones(size(newshape));
        Rreps(Rshape==1)=newshape(Rshape==1);
        

        
        if isL && isR
            
         L=repmat(L,Lreps);
         R=repmat(R,Rreps);
         
         out=ndSparse(bsxfun(fun,L.data,R.data) , newshape);
         return
         
          
        elseif isL
            
             Lsamp=L.data;  Rsamp=R;
             if isempty(Lsamp), Lsamp=reshape(Lsamp,0,1); else Lsamp=Lsamp(1); end
             if isempty(Rsamp), Rsamp=reshape(Rsamp,0,1); else Rsamp=Rsamp(1); end
   
              isOutFull=~issparse( bsxfun(fun, Lsamp, mkCompat(Rsamp) )  );

   
              if isOutFull%cheapest to do everything as full array operations

                  
                 out=bsxfun(fun,full(L),full(R)); return        

   
              else%sparse output - permutation approach is probably cheaper    
                 
                L=repmat(L,Lreps);
                
                idx=(Rshape~=1);
              
                order=1:N;
                  order=[order(idx), order(~idx)];   
                  
                
                L=permute(L,order);  
                R=permute(R,order);  
                
                data=bsxfun(fun,reshape(L.data,numel(R),[]),mkCompat(R(:)));
                
                out=ndSparse( data , newshape(order));
                
                out=ipermute(out,order); 
                
                return
                
              end
          
              

        else%isR
            
          
             Lsamp=L;  Rsamp=R.data;
             if isempty(Lsamp), Lsamp=reshape(Lsamp,0,1); else Lsamp=Lsamp(1); end
             if isempty(Rsamp), Rsamp=reshape(Rsamp,0,1); else Rsamp=Rsamp(1); end
             
              isOutFull=~issparse( bsxfun(fun, mkCompat(Lsamp), Rsamp )  );

              

              

              if isOutFull%cheapest to do everything as full array operations

                  
                 out=bsxfun(fun,full(L),full(R)); return        

   
              else%sparse output - permutation approach is probably cheaper    
                 
                R=repmat(R,Rreps);
                
                idx=(Lshape~=1);
              
                order=1:N;
                  order=[order(idx), order(~idx)];   
                  
                
                L=permute(L,order);  
                R=permute(R,order);  
                
                
                data=bsxfun(fun,mkCompat(L(:)),  reshape(R.data,numel(L),[]));
                
                out=ndSparse( data , newshape(order));
                
                out=ipermute(out,order); 
                
                return
                
              end
              
          
        end
        
        
        
    end
    
    
    function obj=squeeze(obj)
    %squeeze - same syntax as for full arrays
    
        ndShape=obj.ndShape;
        ndShape(ndShape==1)=[];
        obj=reshape(obj,ndShape);
        
    end
    
 
    function [obj,N]=shiftdim(obj,N)
    %shiftdim - same syntax as for full arrays
    
       
       if nargin>1 
           
             order=1:ndims(obj);

             if N>=0
              order=circshift(order,[0,-N]);
              obj=permute(obj,order);
             else
              obj=reshape(obj, [ones(1,-N) ,obj.ndShape]);
             end
  
         
       else
           
           sz=obj.ndShape;
           N=find(sz~=1,1,'first');
           obj=reshape(obj,sz(N:end));
           N=N-1;
           
       end
        
       
    end
    
    

    
    
    function objnew=flipdim(obj,dim)
    %flipdim - same behavior as for full arrays
    
       [ndSubs,vals,ndShape]=getEntryTable(obj,max(ndims(obj), dim));
       
       ndSubs{dim}=ndShape(dim)+1 - ndSubs{dim};
       nzm=nzmax(obj);
       
       objnew=ndSparse.build([ndSubs{:}], vals, ndShape,nzm);
        
    end
  
    function objnew=flipud(obj)
    %flipud - same behavior as for full arrays, except that it works
    %for 3D and higher-dimensional arrays as well. That is,
    %flipud(obj)=flipdim(obj,1)
    
      objnew=flipdim(obj,1);
        
    end    
    
    function objnew=fliplr(obj)
    %fliplr - same behavior as for full arrays, except that it works
    %for 3D and higher-dimensional arrays as well. That is,
    %fliplr(obj)=flipdim(obj,2)
    
      objnew=flipdim(obj,2);
        
    end    
    
    function objnew=rot90(obj,K)
    %rot90 - same behavior as for full arrays, except that it works
    %for 3D and higher-dimensional arrays as well. The rotation will
    %be applied to the first two dimensions.
    
       if nargin<2, K=1; else K=mod(K,4); end
    
     

       switch K
       
           case 0
               
               objnew=obj;
               
           case 1
               
               order=1:ndims(obj);
                   order([1 2])=order([2 1]);
               
               objnew=permute(obj,order);
               objnew=flipud(objnew);
               
           case 2
               
               objnew=flipud(fliplr(obj));
               
               
           case 3
               
               order=1:ndims(obj);
                   order([1 2])=order([2 1]);
               
               objnew=permute(obj,order);
               objnew=fliplr(objnew);
       
       end
       
       
       
        
    end
    
    
 
    function obj=convn(varargin)
    %convn  method for ndSparse
    %
    %SYNTAXES:
    %
    %
    %(1) OBJ=CONVN(A,B,SHAPE) performs the N-dimensional convolution of
    %    arrays A and B. By MATLAB dispatching rules, at least one of which
    %    A and B will be ndSparse.
    %
    %(2) OBJ=CONVN(K_1,K_2,...K_N,A,SHAPE) performs separable convolution of
    %    array A with kernels K_1,...,K_N. Each K_i must be a vector and will
    %    be applied along dimension i of A. Whether K_i is a row or column
    %    vector has no impact here on the effect of the convolution. 
    %    For this syntax to be triggered, it is required that N>2 and that
    %    N=ndims(A)-isvector(A).
    %
    %The SHAPE input argument is optional and works like MATLAB's usual
    %convolution functions:
    %
    %
    %       'full'   - (default) returns the full N-D convolution
    %       'same'   - returns the central part of the convolution that
    %                  is the same size as A.
    %       'valid'  - returns only the part of the result that can be
    %                  computed without assuming zero-padded arrays.
    %
    %
    %CAUTION: The output of this routine will always be ndSparse and is intended
    %for situations where the output data is expected to be sparse
    %intrinsically, e.g., when sparse data is convolved with small or
    %sparse kernels. If there is reason to think the output will not be sparse,
    %it would be advisable to pre-cast to full arrays and use the native
    %MATLAB convn method.
    
    
    
        if ischar(varargin{end})
           truncType=varargin{end}; varargin(end)=[];
        else
           truncType='full';  
        end
    
        B=varargin{end}; 
          varargin(end)=[];
        B=ndSparse(B);
        sizeB=size(B); 
        ndimsB=length(sizeB);
        
        A=varargin{1};
        
        nv=length(varargin);
        sepFlag=(nv>1);%separable input kernel

        if sepFlag
            

            nvReq=ndimsB-isvector(B);
            
            if nv~=nvReq
              error('The number N of separable kernels K_i in CONVN(K_1,K_2,...,K_N,A) must equal N=ndims(A)-isvector(A)');      
            end
            
            if ~all(cellfun(@(c) isvector(c) , varargin) );
             error 'Separable kernel data appears to have been input. All the 1D kernels need to be vectors'
            end
            
            varargin=cellfun(@(c) c(:), varargin,'uni',0);
       
    
            for ii=2:nv
                A=kron(varargin{ii},A);
            end
            
            
            Lengths=cellfun(@(c) length(c),varargin);
            A=reshape(A,Lengths);
            
            A=ndSparse(A);
            [A,B]=deal(B,A);
       
        else

           A=ndSparse(A);
        
        end
 
        ndimsA=ndims(A);
        
        N=max(ndimsA,ndimsB);
       

        [Acoords,Avals,Ashape]=getEntryTable(A,N); 
        
        [Bcoords,Bvals,Bshape]=getEntryTable(B,N);       
        
         convShape=Ashape+Bshape-1;
         
         switch truncType
            
            case 'full'
                %Nothing to do
                
               lb=[]; ub=[];
                
            case 'same'
               
                lb=floor(Bshape/2)+1;
                ub=Ashape+lb-1;
                
            case 'valid'
            
                lb=Bshape;
                ub=convShape+1-Bshape;
                
                if any(lb>ub) %result will be empty
                    convShape=ub-lb+1;
                    convShape(convShape<0)=0;
                    obj=ndSparse.build([],[],convShape );
                    return
                end
                
            otherwise
                
                error 'SHAPE must be ''full'', ''same'', or ''valid''.'
         end
         
           Bcoords=[Bcoords{:}];
        
           Corner=cellfun(@(c) min(c), Acoords,'uni',0);
           Acoords=[Acoords{:}];
           
           S.type='()'; S.subs=Corner;
           cornerval=full(subsref(A,S));
           
           if ~isequal(cornerval,0)
             Acoords=[[Corner{:}];Acoords];
             Avals=[0;Avals];
           end
               


        
        convVals=Avals(:)*Bvals(:).';

        convCoords=bsxfun(@plus, reshape(Bcoords.',N,1,[]) ,Acoords.'-1);
        
         convCoords=reshape( convCoords,N,[]).';
        
     
        
       if ~isempty(lb)%Only necessary for 'same' and 'valid'
          
           idx = all( bsxfun(@ge,convCoords,lb) &  bsxfun(@le,convCoords,ub) ,2);

           convCoords=bsxfun(@minus, convCoords(idx,:)  , lb-1);
           convVals=convVals(idx);
           convShape=ub-lb+1;
           
       end
       
       obj=ndSparse.build(convCoords,convVals,convShape);
       
    end
    
    
    
    
 end%NONSTATIC NON HIDDEN METHODS
   
 
   methods (Hidden)
      function data=getData(obj) 
          data=obj.data;
      end

      function ndShape=getndShape(obj) 
          ndShape=obj.ndShape;
      end
      
       function objnew=sumEngine(memliberal, obj,varargin)
       %sumEngine
       
         memliberal=memliberal | ndSparse.oldstyle;  
           
         if ndims(obj)<=2 %Matrix case
             objnew=ndSparse( sum(obj.data, varargin{:}) );
             return;
         end
         
            if  isempty(varargin) || ischar(varargin{1}) 
               dim=find(obj.ndShape,1,'first');
            else
               dim=varargin{1};
            end
        
            sz=size(obj);
            nn=length(sz);
            if dim>nn || sz(dim)==1,objnew=obj; return; end
            
            if memliberal

                data=ExtractDim1Reshaping(obj.data, dim, sz);

                 data=sum(data,1);
                 sz(dim)=1;

                data=invExtractDim1Reshaping(data,dim,sz);
                objnew=ndSparse(data,sz);

            else

                [ndSubs,vals]=getEntryTable(obj,nn);

                ndSubs{dim}(:)=1;
                sz(dim)=1;
                
               
                 objnew=ndSparse.accumarray([ndSubs{:}],vals,sz);           

                
            end
            
            
       end  
      
      
       function objnew=anyEngine(memliberal,obj,varargin)
       %anyEngine 
       
       
         memliberal=memliberal | ndSparse.oldstyle;  
       
         if ndims(obj)<=2 %Matrix case
             objnew=ndSparse( any(obj.data, varargin{:}) );
             return;
         end
         
            if  isempty(varargin) 
               dim=find(obj.ndShape,1,'first');
            else
               dim=varargin{1};
            end
        
            sz=size(obj);
            nn=length(sz);
            if dim>nn || sz(dim)==1,
                objnew=obj;
                if ~islogical(objnew.data)
                 objnew.data=logical(objnew.data);
                 end
                return; 
            end
         
            if memliberal
            
                data=ExtractDim1Reshaping(obj.data,dim,sz);

                 data=any(data,1);
                 sz(dim)=1;

                data=invExtractDim1Reshaping(data,dim,sz);

                objnew=ndSparse(data,sz);

            else
            
              objnew =  sumEngine(0, obj~=0 ,dim)>0; 
            
            end

       end  

       function objnew=allEngine(memliberal,obj,varargin)
       %allEngine
       
         
         memliberal=memliberal | ndSparse.oldstyle;  
           
         if ndims(obj)<=2 %Matrix case
             objnew=ndSparse( all(obj.data, varargin{:}) );
             return;
         end
         
            if  isempty(varargin) 
               dim=find(obj.ndShape,1,'first');
            else
               dim=varargin{1};
            end
        
            sz=size(obj);
            nn=length(sz);
            if dim>nn || sz(dim)==1,
                objnew=obj;
                if ~islogical(objnew.data)
                 objnew.data=logical(objnew.data);
                 end
                return; 
            end
  
            if memliberal
            
                data=ExtractDim1Reshaping(obj.data,dim,sz);

                 data=all(data,1);
                 sz(dim)=1;

                data=invExtractDim1Reshaping(data,dim,sz);

                objnew=ndSparse(data,sz);
            
            else
            
                objnew=( sumEngine(0, obj~=0 , dim)  == sz(dim)  ); 
            
            end
       end          
       
      
       function objnew=meanEngine(memliberal, obj,varargin)
       %meanEngine  
           
         memliberal=memliberal | ndSparse.oldstyle;  
       
         if ndims(obj)<=2 %Matrix case
             objnew=ndSparse( mean(obj.data, varargin{:}) );
             return;
         end
         
            if  isempty(varargin) 
               dim=find(obj.ndShape,1,'first');
            else
               dim=varargin{1};
            end
        
            sz=size(obj);
            nn=length(sz);
            if dim>nn || sz(dim)==1,objnew=obj; return; end
            
            if memliberal
            
             data=ExtractDim1Reshaping(obj.data,dim,sz);
            
              data=mean(data,1);
              sz(dim)=1;
            
             data=invExtractDim1Reshaping(data,dim,sz);
        
            
             objnew=ndSparse(data,sz);
 
            else
            
              objnew=sumEngine(0,obj,dim)/sz(dim); 
            
            end  
              
              
       end 
    
       
    function objnew=circshiftEngine(memliberal,obj,circDims)
    %circshiftEngine 
      
    
        memliberal=memliberal | ndSparse.oldstyle;  

        ndShape=size(obj);
        nn=length(ndShape);
        
        if nn==2
         obj.data=circshift(obj.data,circDims);
         objnew=obj;
         return;
        end
            
        N=max(length(circDims),nn);
        circDims=[circDims,ones(1,N-length(circDims))];%BUG? should be zeros()?
    
        
        sz=trail1s(ndShape,N);
        initialshape=sz;
            
        map=(circDims~=0);
        firstShift=find(map,1,'first');
         if isempty(firstShift) || firstShift>nn,
             objnew=obj;return
         end
      
       if memliberal 
           
           %reformulate
           sz=circshift(sz,[0,1-firstShift]);
           circDims=circDims(firstShift:end);
  
         
              data=ExtractDim1Reshaping(obj.data,firstShift,initialshape);
              mm=length(circDims);

            for ii=1:mm

               data=reshape(data,sz(ii),[]);

               data=data.';

               thisShift=circDims(ii);

               if thisShift
                data=circshift(data,[0,thisShift]);
               end           

            end


            objnew=ndSparse(data, initialshape);
        
        
       else%~memliberal
          
           
         [ndSubs,vals]=getEntryTable(obj,nn);
         
         if isempty(vals); objnew=obj; return; end %all zero array
         
         circDims=circDims(1:nn);
         
         for ii=1:nn
          if circDims(ii)==0, continue; end   
          ndSubs{ii}=mod(ndSubs{ii}+circDims(ii)-1,sz(ii))+1;
         end
         
         objnew=ndSparse.build(ndSubs,vals,sz,nzmax(obj));
           
       end
        
        
    end      
       
    
       function objnew=catEngine(memliberal,dim,varargin)
       %catEngine - same syntax as for full arrays
       
       
         memliberal=memliberal | ndSparse.oldstyle;  
       
       
          Nargs=length(varargin);
          
          map=cellfun('isclass',varargin,'ndSparse');
          first=find(map,1); 
          if first==1
             L=varargin{first};
             if Nargs==1, objnew=L; return; end
             R=varargin{first+1};
             LastProcessed=2;
          else
             L=cat(dim,varargin{1:first-1});
             R=varargin{first};
             LastProcessed=first;
          end
          
          [szL,ndimsL]=trail1s(size(L),dim);
          [szR,ndimsR]=trail1s(size(R),dim);
       
          idx=true(1,length(szL));
               idx(dim)=0;
          if ndimsL~=ndimsR || any(szL(idx)~=szR(idx))
              error 'CAT arguments dimensions are not consistent.'   
          else
            newshape=szL;
            newshape(dim)=newshape(dim)+szR(dim);
          end
         
          
          
         if dim>=ndimsL || (ndimsL<=2) %equivalent to 2D cat

             
            L=mkCompat(L,szL);
            R=mkCompat(R,szR);
              
            objnew=ndSparse( cat(min(dim,2),L,R),newshape);
             
         else %Fully Multi-dimensional case
         
           if memliberal
          
                         
            L=mkCompat(L,szL);
            R=mkCompat(R,szR);   
               
               
            L=ExtractDim1Reshaping(L,dim,szL);
            R=ExtractDim1Reshaping(R,dim,szR);
             

             L=[L.',R.'].';

             L=invExtractDim1Reshaping(L,dim,newshape);
             
             objnew=ndSparse(L,newshape);
            
           else
              
               if ~map(1), L=ndSparse(L); end
               if ~map(2), R=ndSparse(R); end              
                            
               [LSubs,Lvals]=getEntryTable(L,ndimsL);
               [RSubs,Rvals]=getEntryTable(R,ndimsR);
                
               RSubs{dim}=RSubs{dim}+szL(dim);
                              
               objnew=ndSparse.build( [[LSubs{:}];[RSubs{:}]] , [Lvals;Rvals] , newshape);
               
               
           end
             
         end
         
         %Use recursion to handle further CAT args
         if LastProcessed<Nargs
            objnew=catEngine(memliberal,dim,objnew,varargin{LastProcessed+1:Nargs}); 
         end
         
       end

    
       function varargout=minEngine(memliberal,L,varargin)
       %minEngine
          
       memliberal=memliberal | ndSparse.oldstyle;  
       
         if nargin>2,
            R=varargin{1}; 
         else
            R=[];
         end        

         
         if nargin>3
             dim=varargin{2};
         elseif isempty(R); 
             dim=find(L.ndShape,1,'first');     
         end
         
         if isempty(R)
             
            sz=size(L);
            nn=length(sz);
            if dim>nn || sz(dim)<=1,
                varargout{1}=L;
                if nargout>1,
                   varargout{2}=ones(size(L)); 
                end
                return; 
            end
            
            
            
          if memliberal  
            
                data=ExtractDim1Reshaping(L.data,dim,sz);

                [varargout{1:nargout}]=min(data,[],1);
                sz(dim)=1;

                varargout{1}=invExtractDim1Reshaping(varargout{1},dim,sz);
                  varargout{1}=ndSparse(varargout{1},sz); 

                if nargout>1  
                varargout{2}=invExtractDim1Reshaping(varargout{2},dim,sz);
                  varargout{2}=reshape(varargout{2},sz);
                end
            
          else
              
                data=ExtractDimNReshaping(L.data,dim,sz);

                [varargout{1:nargout}]=min(data,[],2);
                sz(dim)=1;

                varargout{1}=invExtractDimNReshaping(varargout{1},dim,sz);
               
                if nargout>1  
                 varargout{2}=invExtractDimNReshaping(varargout{2},dim,sz);
                end  

              
          end
            
            
            
            
            
         elseif length(varargin)>1 %R and dim given
             
          error('MIN with two matrices to compare and a working dimension is not supported.')
          
         elseif nargout>1%R given, dim not given, but 2nd argout requested
           error('MIN with two matrices to compare and two output arguments is not supported.')
         
         else%bi-operand min 
         
             
          %varargout{1}=bsxfun(@min,L,R);   
             
           [s,L,R]=getshape(L,R);
           varargout{1} = finalObject( min(mkCompat(L,s),mkCompat(R,s)),s); 
             
         end
        
      end
       
      function varargout=maxEngine(memliberal,L,varargin)
      %maxEngine
          
          
          memliberal=memliberal | ndSparse.oldstyle;
          
          if nargin>2,
              R=varargin{1};
          else
              R=[];
          end
          
          
          
          if nargin>3
              dim=varargin{2};
          elseif isempty(R);
              dim=find(L.ndShape,1,'first');
          end
          
          if isempty(R)
              
              sz=size(L);
              nn=length(sz);
              if dim>nn || sz(dim)<=1,
                  varargout{1}=L;
                  if nargout>1,
                      varargout{2}=ones(size(L));
                  end
                  return;
              end
              
              if memliberal
                  
                  data=ExtractDim1Reshaping(L.data,dim,sz);
                  
                  [varargout{1:nargout}]=max(data,[],1);
                  sz(dim)=1;
                  
                  varargout{1}=invExtractDim1Reshaping(varargout{1},dim,sz);
                  varargout{1}=ndSparse(varargout{1},sz);
                  
                  if nargout>1
                      varargout{2}=invExtractDim1Reshaping(varargout{2},dim,sz);
                      varargout{2}=reshape(varargout{2},sz);
                  end
                  
              else
                  
                  data=ExtractDimNReshaping(L.data,dim,sz);
                  
                  [varargout{1:nargout}]=max(data,[],2);
                  sz(dim)=1;
                  
                  varargout{1}=invExtractDimNReshaping(varargout{1},dim,sz);
                  
                  if nargout>1
                      varargout{2}=invExtractDimNReshaping(varargout{2},dim,sz);
                  end
                  
                  
              end
              
              
              
          elseif length(varargin)>1 %R and dim given
              
              error('MAX with two matrices to compare and a working dimension is not supported.')
              
          elseif nargout>1%R given, dim not given, but 2nd argout requested
              error('MAX with two matrices to compare and two output arguments is not supported.')
              
          else%bi-operand max
              
              %varargout{1}=bsxfun(@max,L,R);
              
              [s,L,R]=getshape(L,R);
              varargout{1} = finalObject( max(mkCompat(L,s),mkCompat(R,s)),s);
              
          end
          
      end
      
      
      
   end%hidden methods

  methods (Static)
     
     
     function obj=loadobj(B)
     
         obj=ndSparse(B.data,B.ndShape);
         
     end
         
     function obj=build(varargin)
     %ndSparse.build - This is a generalization of sparse(i,j,s,m,n,nzmax)    
     %to N dimensions, allowing one to build an ndSparse array from a
     %user-supplied table of explicit array entries.
     %
     %SYNTAX:
     %
     %  S = ndSparse.build(explicitCoords,s,ndShape, nzmax)
     %
     %
     %out:
     %
     % S: the N-dimensional sparse output array as an ndSparse object.
     %
     %in:
     %
     % explicitCoords: a PxN matrix whos rows are the N-dimensionsal
     %                 coordinates of the explicit entries in the desired
     %                 output array S.
     %
     %          s:     a P-vector of array values at the locations 
     %                  supplied in explicitCoords. Can also be a scalar.
     %
     %    ndShape: The desired dimensions of S. Defaults to max(explicitCoords,[],1).
     %
     %      nzmax:  The output S will have memory allocated for nzmax
     %              non-zeros. Default = nnz(S).
     %
     %
     %To create an ndSparse array of zeros, one can do
     %
     %  S=ndSparse.build(ndShape)
     %
     %which is a simplifcation of S=ndSparse.build([],[],ndShape,0).
     %
     %See also ndSparse.spalloc, ndSparse.sprand, ndSparse.sprandn

         nzmaxARG={};
         switch nargin
             
             case 1  %ndSparse.build(ndShape)
             
                 ndShape=varargin{1};
                 
                 [M,N,ndShape]=ndShape2MN(ndShape);
                 data=sparse(M,N);
                 obj=ndSparse(data,ndShape);
                 
                 return;
                 
             case 2 %ndSparse.build(explicitCoords,s) 
                 
                 explicitCoords=varargin{1};
                 s=varargin{2};
                 
                 if ~iscell(explicitCoords)
                  ndShape=max(explicitCoords,[],1);
                  %nzmax=size(explicitCoords,1);                
                 else
                   ndShape=cellfun(@max,explicitCoords);
                   %nzmax=length(explicitCoords{1});
                 end
                  
             case 3 %ndSparse.build(explicitCoords,s,ndShape) 
                 
                 [explicitCoords,s,ndShape]=deal(varargin{:});
                 
                    %  if ~iscell(explicitCoords)
                    %    nzmax=size(explicitCoords,1);                
                    %  else
                    %    nzmax=length(explicitCoords{1});
                    %  end
 
             case 4 %ndSparse.build(explicitCoords,s,ndShape,nzmax)   
                 
                [explicitCoords,s,ndShape,nzmaxARG{1}]=deal(varargin{:});
               
                 
             otherwise
                 
                 error 'Too many inputs'
                 
         end
         
         if iscell(explicitCoords)
            emptyCoords=any(cellfun('isempty',explicitCoords));
         else
            emptyCoords=isempty(explicitCoords);
         end
         
         
          if emptyCoords && isempty(s)
              
                    [M,N]=ndShape2MN(ndShape);
                    data=sparse([],[],[], M, N, nzmaxARG{:});
                    obj=ndSparse(data,ndShape); 
                    return 
          end
         
          
          [I,J,dims]=ndCoords2IJ(explicitCoords,ndShape);
 
          S=s;
          M=dims(1);
          N=dims(2); 
          
          data=sparse(I,J,S,M,N,nzmaxARG{:});  
          
          
          obj=ndSparse(data,ndShape);
          
     end
     
     function obj=spalloc(ndShape,nzmax)
     %spalloc method (Static) for ndSparse
     %
     %SYNTAX:
     %
     %  obj=ndSparse.spalloc(ndShape,nzmax)
     %
     %in:
     %
     % ndShape: desired dimensions of the output ndSparse array obj
     % nzmax : desired capacity for nonzeros
     %
     %out: output all-zero ndSparse object of dimensions ndShape 
    
         obj=ndSparse.build([],[],ndShape,nzmax);
         
     end
     
     
     function obj=sprand(ndShape,density)
     %sprand method (Static) for ndSparse
     %
     %SYNTAX:
     %
     %  obj=ndSparse.sprand(ndShape,density)
     %
     %in:
     %
     % ndShape: desired dimensions of the output ndSparse array obj
     % density : desired density of nonzeros
     %
     %out: output ndSparse object of dimensions ndShape 
     
         [mm,nn,ndShape]=ndShape2MN(ndShape);
         
         obj=ndSparse( sprand(mm,nn,density) , ndShape);
         
     end
     
     
     function obj=sprandn(ndShape,density)
     %sprandn method (Static) for ndSparse
     %
     %SYNTAX:
     %
     %  obj=ndSparse.sprand(ndShape,density)
     %
     %in:
     %
     % ndShape: desired dimensions of the output ndSparse array obj
     % density : desired density of nonzeros
     %
     %out: output ndSparse object of dimensions ndShape 
     
         [mm,nn,ndShape]=ndShape2MN(ndShape);
         
         obj=ndSparse( sprandn(mm,nn,density) , ndShape);
         
     end
     
     
     function obj=accumarray(subs,val,sz,fun)
     %ndSparse.accumarray 
     %
     %SYNTAX:
     %
     %  obj = ndSparse.accumarray(subs,val,sz,fun)
     %
     %This works pretty much like MATLAB's built-in accumarray(subs,val,sz,fun)    
     %except that the output obj is always ndSparse. Defaults for sz and fun
     %are the same as for the native accumarray as well.
     %
     %Note that the normal built-in accumarray also supports a 5th and 6th argument
     %FILLVAL and an ISSPARSE, which are not supported here. They make no sense
     %for an ndSparse output obj, since it will always be sparse, and hence also the only
     %sensible FILLVAL is zero. 

         
         if nargin==4, do_fun=true; else do_fun=false; end
         
         if iscell(subs)
            
             msg='Cells of first input SUBS must contain real, full, numeric vectors of equal length.';
             
             subs=cellfun(@(c) c(:), subs, 'uni', 0);
             try
              subs=[subs{:}];
             catch
               error(msg)
             end
             
             if ~isreal(subs) || issparse(subs) || ~isnumeric(subs)
               error(msg) 
             end 
             
         end
 
         
          N=size(subs,2);
          
          if ~isreal(subs) || issparse(subs) || ~isnumeric(subs) || N<1
           error 'First input SUBS must be a real, full, numeric matrix with at least one column.'
          end
          
          minsz=max(subs,[],1); 
          if isempty(minsz), minsz=zeros(1,N); end
          if isscalar(minsz), minsz=[minsz,1]; end
          
          if nargin<3 || isempty(sz), 
              sz=minsz;
          end
          
          if N>1 && length(sz)~=N
            error 'Third input SZ must be a full row vector with one element for each column of SUBS.'
          elseif N==1 && length(sz)<2
              error 'When input SUBS is a column vector, input SZ must be of the form [M,1]'
          end

          if ~all(sz>=minsz)
              error 'First input SUBS and third input SZ must satisfy ALL(MAX(SUBS)<=SZ).'
          end
             

         if isempty(subs)
            obj=ndSparse.build(sz); return 
         end
          
         desiredShape=sz;
         
         [ii,jj,dims]=ndCoords2IJ(subs,desiredShape);

         
         if ~isa(val,'double')
           val=double(val);    
         end
         
         if do_fun
            data=accumarray([ii,jj],val,dims,fun,0,true);
         else
            data=accumarray([ii,jj],val,dims,[],0,true);
         end
 
         obj=ndSparse(data,desiredShape);
         
     end
     
 end%STATIC METHODS
 
end






 function [E,targetshape,boolGrow,other]=nd2matrixIndex(ndSubs,ndShape,operation)
 %nd2matrixIndex 
 %
 % ndShape is expected to be untrailed by ones
 
     N=length(ndSubs);
     [ndShape, numdims]=untrail1s(ndShape);
    
     
     [objIs1D,obj1Ddim]=proc1D(ndShape);
     
     LinearIndexing=(N==1);
     
     boolConsol = 1<N && N<numdims;   %consolidating indexing of trailing dims.
     
     
     if ~LinearIndexing%subscript indexing
         
      ShapeCell=parseSize(ndShape,N);
      quasiShape=[ShapeCell{:}];
      
     else%linear indexing
         
       idx1=ndSubs{1};  
       quasiShape=prod(ndShape);
    
     end
  
      ndSubs0=ndSubs;%initial substruct
      ndSubs=cellfun(@mkCompat, ndSubs,'uni',0);
      targetshape=zeros(1,N);
      accessedshape=zeros(1,N);
      subshape=zeros(1,N);
      boolGrow=false;
               
          for ii=1:N %loop over indices
               

           idx=ndSubs{ii};

               if islogical(idx) %logical index
                   
                   idx=find(idx);
                   maxAccessed=max(idx);
                   if isempty(maxAccessed), maxAccessed=quasiShape(ii); end
                   
                   subshape(ii)=numel(idx);
                   
               elseif ischar(idx) && ~LinearIndexing %colon index
                   
                   idx=(1:quasiShape(ii)).';
                   maxAccessed=quasiShape(ii);
                   
                   subshape(ii)=maxAccessed;
                   
               elseif ischar(idx) && LinearIndexing %colon index                  
                   
                   idx=':';
                   maxAccessed=quasiShape(ii);                
                   
                   subshape(ii)=maxAccessed;
                   
               else %numeric index
                   
                   maxAccessed=max(idx);
                   
                   subshape(ii)=numel(idx);
               end
               
              
                  
                    switch operation

                       case 'subsref'

                         if maxAccessed>quasiShape(ii)
                             error  'Index exceeds matrix dimensions.'
                         end


                         targetshape(ii)=subshape(ii);
                         accessedshape(ii)=quasiShape(ii);

                       case 'subsasgn'

                         accessedshape(ii)=max(maxAccessed,quasiShape(ii));
                         targetshape(ii)=accessedshape(ii);
                         
                         
                         if maxAccessed>quasiShape(ii)
                           boolGrow=true;

                         end
                         
                        otherwise
                           
                         error 'Unknown operation'

                    end

            ndSubs{ii}=idx;
 
          end
      
       if ~LinearIndexing,   %Clean up trailing dimensions
          
        [accessedshape,N]=untrail1s(accessedshape);
        targetshape=targetshape(1:N);     
        subshape=subshape(1:N);
        ndSubs=ndSubs(1:N);
        
       end

          
          if boolGrow
                  
            if LinearIndexing && ~objIs1D %linear indexing of non-vector

                if numdims==2
                  error 'In an assignment  A(I) = B, a matrix A cannot be resized.'
                elseif ndims>2
                  error 'Attempt to grow array along ambiguous dimension.'
                end
                
            elseif boolConsol  
            
                error 'Attempt to grow array along ambiguous dimension.'
                
            end
            
          elseif isequal(operation,'subsasgn')
              
              targetshape=ndShape;

          end    
          
        
                  
          if LinearIndexing && isequal(operation,'subsref') %linear indexing in subsref
              
                   idx1=ndSubs{1};
                   if ~ischar(idx1)
                    idxShape=size(idx1);
                   else%single colon
                     idxShape=[quasiShape(1),1];   
                   end
                   idxIs1D=proc1D(idxShape);
                   
                   
                   
                if  idxIs1D && objIs1D
                    
                    tmp=ndShape;
                    tmp(obj1Ddim)=targetshape;
                    targetshape=tmp;

                else  % shape of idx1 dictates resulting shape

                    targetshape=idxShape;
                    
                end

          end
        

           
           
       
      %%Everything up to here was to convert indices to numeric form
      %%and to compute targetshape and accessedshape - now prepare
      %%subscripts
            
 
      E.type='()';
      E.subs{1}=ndSubs{1};   

      multicolons=false; 
      switch N
          
             
          case 1 %linear indexing, use index as is.




          case 2
       
               E.subs{2}=ndSubs{2};       
               
          otherwise%N>2
              

            E.subs{2}=ndSubs{N};
            
            %=numdims; 
            nnn=N;

            if cellfun('isclass', ndSubs0(1:nnn-1),'char') 
                 %Colons in the first numdims-1 subscripts
              
              multicolons=true;   
              E.subs{1}=':';  
                
            else     
                
              E.subs{1}=sub2allind(accessedshape(1:N-1),ndSubs{1:N-1});
              
            end
            
            
            
      end             
              
   other.subshapeND=subshape;
   other.subshape=cellfun('length',E.subs);
    if LinearIndexing, 
        other.subshape=subshape;
    elseif multicolons
        other.subshape(1)=prod(ndShape(1:N-1));
    end
   other.boolConsol=boolConsol;
   other.quasiShape=quasiShape;
   other.LinearIndexing=LinearIndexing;
   
   
   
   
   
 end 
  
 
 
 function outcell=parseSize(sz,numargsout,dim)
 %parseSize


     idx=find(sz~=1,1,'last');
     if isempty(idx)
         sz=[1,1];
     elseif idx==1
         sz=[sz(1), 1];
     else
         sz=sz(1:idx);
     end
         
       
    dimSpecified=logical(exist('dim','var'));


    if dimSpecified,

      if ~isscalar(dim), error('DIM argument must be scalar.'); end

      if dim>length(sz), sz(end+1:dim)=1; end

      sz=sz(dim);

    end


    if numargsout<=1,

        outcell{1}=sz;
        return;

    elseif dimSpecified,

      error('Too many output arguments.')

    else


      if length(sz)>=numargsout

        sz(numargsout)=prod(sz(numargsout:end));
        sz(numargsout+1:end)=[];

      else

        sz(end+1:numargsout)=1;  

      end

      outcell=num2cell(sz);

    end

 end
 
 
 function X=mkCompat(X,s)
 %mkCompat
 
     if isa(X,'ndSparse')
        X=X.data;
     elseif ~(isa(X,'double') || islogical(X) || ischar(X))
          %char is for handling the case where X is the index ':'
          
         X=double(X);
     end
     
     if nargin>1 && ~isscalar(X)
        X=reshape(X,[],s(end)); 
     end
     
 end

 
 function out=finalObject(obj, ndShape) 
 %finalObject 
 
      if nargin<2,
          out=ndSparse(obj);
      elseif ~issparse(obj)%a full object type is produced (not ndSparse or native sparse)
          out=reshape(obj,ndShape);
      else
          out=ndSparse(obj,ndShape); 
      end
 
 end
 
function [shape,L,R]=getshape(L,R)
%getshape

   szL=size(L);
   szR=size(R);
   
   shape=szL; %Initial assignment
   
   %Overrides
   if isscalar(L)
     shape=szR;
   elseif isscalar(R)
       shape=szL;
   elseif ~isequal(szL,szR)
       
     N=max(numel(szL),numel(szR));
     
     szL=trail1s(szL,N); szR=trail1s(szR,N);
     
     
     inExpandable=any(szL~=szR & min(szL,szR)~=1);  
       
     if inExpandable  
      error 'Array dimensions must match for binary array op.' 
     else %"implicit" expansion
         
        shape=max(szL,szR);

        L=repmat(L,shape./szL);
        R=repmat(R,shape./szR);

     end
      
      
   end
       
end

function data=ExtractDim1Reshaping(data,dim,ndShape)
%ExtractDim1Reshaping

      %ndShape=obj.ndShape;   
      %data=obj.data;
      
      nd=length(ndShape); 
      
     if dim==1
         
         data=reshape(data,ndShape(1),[]);
   
     elseif dim==nd
         
         data=data.';
      
     elseif dim<nd %it is an implied condition here that dim>1
         
          data=reshape(data, prod(ndShape(1:dim-1)),[]); 
          data=reshape(data.',ndShape(dim),[]);
           
     else% dim>nd
         
        %ndShape=trail1s(ndShape,dim);  
         error 'DIM>nD. This situation should be handled by separate processing'
      
     end

    
  
end

function data=invExtractDim1Reshaping(data,dim,newshape)
%invExtractDim1Reshaping - not a precise inverse, only up to a reshaping


     nd=length(newshape);

     if dim==1
         
        
   
     elseif dim==nd
         
         data=data.';
      
     elseif dim<nd %it is an implied condition here that dim>1
         

         data=reshape( data, [], prod( newshape(1:dim-1) ) ).' ;
          
          
     else% dim>nd
         
         error 'DIM>nD. This situation should be handled by separate processing'
      
     end

  
    %obj=ndSparse(data,newshape);
  
  
end


function data=ExtractDimNReshaping(data,dim,ndShape)
%ExtractDimNReshaping


      
      nd=length(ndShape); 
      
     if dim==nd
         
         %do nothing
   
      
     elseif dim<nd 
         
         
       order=[1:dim-1, dim+1:nd, dim];
       obj=ndSparse(data,ndShape);
       obj=permute(obj,order);
       
       
       data=obj.data; 
       if ndShape(dim)==1
         data=data(:);  
       end
        
     else% dim>nd
         
         error 'DIM>nD. This situation should be handled by separate processing'
      
     end

    
    
  
end

function obj=invExtractDimNReshaping(data,dim,newshape)
%invExtractDimNReshaping
 

      nd=length(newshape); 
      
     if dim<=nd
         

         
       order=[1:dim-1, dim+1:nd, dim];
       obj=ndSparse(data,newshape(order));
       obj=ipermute(obj,order);
           
     else% dim>nd
         
         error 'DIM>nD. This situation should be handled by separate processing'
      
     end

    
  
end




function [ndShape,N]=trail1s(ndShape,dim)
%trail1s

  N=length(ndShape);
  if N==1, 
      ndShape=[ndShape,1]; 
      N=2; 
  end
  if dim>N, %add trailing ones
     ndShape=[ndShape, ones(1,dim-N)]; 
     N=dim;
  end 

end

function [ndShape,N]=untrail1s(ndShape)
%untrail1s

  idx=find(ndShape~=1,1,'last');

   if isempty(idx)
       ndShape=[1,1];
   elseif idx==1
       ndShape=[ndShape(1),1];
   else
       ndShape=ndShape(1:idx);
   end
   
   N=length(ndShape);

end

function [MM,NN,ndShape]=ndShape2MN(ndShape)
%ndShape2MN
    
  ndShape=untrail1s(ndShape);

  MM=prod(ndShape(1:end-1));
  NN=ndShape(end);
  
     
end


function  [ndSubs,vals,ndShape]=getEntryTable(obj,maxdim)
%getEntryTable
%
%[ndSubs,vals,ndShape]=getEntryTable(obj,maxdim)
%
%Produces a cell array ndSubs{1:maxdim} of N-dimensional
%subscripts of explicit entries of obj and a corresponding array, 
%val, of values. The output NDSHAPE will be the dimension of obj with any 
%required trailing 1 dimensions. 

         N=maxdim;
         ndShape=trail1s(obj.ndShape,N);


         [II,JJ,vals]=find(obj.data); 
          
         ndSubs=IJ2ndCoords(II,JJ,ndShape); 

         vals=vals(:);
         
end




function [II,JJ,dims]=ndCoords2IJ(ndCoords,ndShape)
%ndCoords2IJ


        if ~iscell(ndCoords),  ndCoords=num2cell(ndCoords,1); end
         
        if any(cellfun('isempty',ndCoords))
            II=[]; JJ=[];  
            dims=[ prod( ndShape(1:end-1) ) , ndShape(end)];
            return
        end
         
         [ndShape,N]=trail1s(ndShape,length(ndCoords));
         
         minShape=trail1s( cellfun(@max, ndCoords)  , N );

         if any(minShape>ndShape)
            error 'Given array ndShape is not large enough to hold given n-D coordinates.' 
         end
         
         [ndShape,numdims]=untrail1s(ndShape);
          if ndShape(2)==1 && numdims<=2, 
              numdims=1; 
          end
         ndCoords=ndCoords(1:numdims);
         
         
          if numdims==1
              II=ndCoords{1}(:);
              JJ=ones(length(II),1);
              dims=ndShape;
          elseif numdims==2
              II=ndCoords{1}(:);
              JJ=ndCoords{2}(:);
              dims=ndShape;
          else

             II=sub2ind(ndShape(1:end-1),ndCoords{1:end-1});
               II=II(:);
             JJ=ndCoords{end}(:);
             dims=[ prod( ndShape(1:end-1) ) , ndShape(end)];
             
          end
          

end

function [ndCoords,ndShape]=IJ2ndCoords(II,JJ,ndShape)
%IJ2ndCoords

           N0=length(ndShape);
           
           [ndShape,N]=untrail1s(ndShape);

           ndCoords=[ cell( 1, N-1 )  , {JJ} ];
           
           [ndCoords{1:N-1}]=ind2sub(ndShape(1:end-1), II); 
                
           ndCoords=cellfun(@(c) c(:), ndCoords, 'uni',0);
           
           if N0>N
              ndCoords=[ndCoords, num2cell( ones(size(ndCoords{1},1),N0-N) , 1) ];
           end
end



function bool=iscolumn(A)

  sz=size(A);
  bool=length(sz)==2 && sz(2)==1;
  
end

function bool=isrow(A)

  sz=size(A);
  bool=length(sz)==2 && sz(1)==1;
  
end


%%%%%%%%%%%%%%%%%%%%%  EXTERNAL CONTRIBUTIONS  %%%%%%%%%%%%%%%%%%%%%%%%


function idx = sub2allind(sz, varargin)
%   idx = sub2allind(sz, sub1, sub2, sub3, ... )
%
%   Like sub2ind, sub2allind computes the equivalent linear indices for
%   given subscripts of an array with size SZ.
%   Unlike sub2ind, it computes them for all combinations of
%   subscripts. So, instead of calling A( 2:3, 1, 4:11) you might
%   use
%       linIdx = sub2allind( size(A), 2:3, 1, 4:11 );
%
%   and then call A(linIdx).
%
%   Using the colon operator is allowed:
%   linIdx = sub2allind( sz, sub1, :, sub3 );

% Michael Völker, 2011
% This is Matlab FileEx #30096

    % ==================================================================
    % Make sure, "sz" is fine
    %
        if nargin < 2
            error('sub2allind:NoOfInputs', 'At least two inputs, please.')
        end
        if ~exist('sz', 'var') || length(sz(:)) < 2 || ~all(isnumeric(sz(:))) || ~all(isreal(sz(:)))
            error('sub2allind:szNotSane', 'Size vector must be real and have at least 2 elements.');
        end
        sz  = sz(:).';
        Ndims = length(sz);

        if Ndims < nargin-1
            ResDims = nargin-1 - Ndims;
            sz = [sz   ones(1,ResDims)];                    % Adjust for trailing singleton dimensions
        elseif Ndims > nargin-1
            sz = [sz(1:nargin-2)  prod(sz(nargin-1:end))];  % Adjust for linear indexing on last element
        end
    % ==================================================================

    Ndims = length(sz);     % No. of dimensions

    % ==================================================================
    % Make sure all subscripts are fine, including the use of ":"
    %
        for d = 1:Ndims
          flagNum  = isnumeric( varargin{d}(:) );
          flagReal = isreal( varargin{d}(:) );
          flagInt  = all( varargin{d}(:) == floor(varargin{d}(:)) );
          if ~flagNum || ~flagReal || ~flagInt

            if isequal( varargin{d}, ':' )  % We allow exactly one type of non-numeric input
                varargin{d} = 1:sz(d);      % interpret ":" as in usual subscript syntax
            else
                error('sub2allind:SubsNotSane', 'Subscripts must either be valid matrix subscripts or the well known '':''.')
            end

          end
        end
    % ==================================================================


    % ==================================================================
    %  Compute linear indices, e.g. in 3D, with edges L1, L2, L3:
    %
    %       idx = x1 + (x2-1)*L1 + (x3-1)*L1*L2,
    %
    %  for every permutation (x1,x2,x3)
    k = [1 cumprod(sz(1:end-1))];           % k(d) holds the accumulated No. of array elements
                                            % up to the d'th subscript (or dimension)
    idx = 1;                                % smallest possible index
    for d = 1:Ndims

        xd = varargin{d}(:);                        % the d'th subscripts

        if any(xd < 1) || any(xd > sz(d))           % Verify subscripts are within range
            error('sub2allind:IndexOutOfRange', 'Out of range subscript.');
        end

        reshrule      = ones(1,d+1);                % how to reshape the size of xd
        reshrule(end) = length(xd);                 %

        xd = reshape( xd, reshrule );               % prepare for bsxfun's needs

        idx = bsxfun( @plus, idx, (xd-1)*k(d) );    % iteratively calculate the sum pointed out above

    end
    % ==================================================================

    idx = idx(:);       % linearify indices

end

function  [bool,dim] = proc1D(sz)

 dim=find(sz>1);
 ldim=length(dim);
 map=double(sz>0);
 bool= (ldim<=1) & all(map);
  
 if ~bool
    dim=[];
 elseif ldim==0 %means the object is a scalar
     dim=1;
 end
 

end
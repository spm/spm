function varargout=ndstest(TOL)
%Performs numerous tests of ndSparse math operations, 
%
%  ndstest(TOL)
%
%TOL is a tolerance value on the percent error. Execution will pause in debug
%mode for inspection if any one of the tests exhibits an error greater than
%TOL.

if nargin<1
 TOL=inf; %default tolerance value on discrepancies
end

CHECKTYPES=false;

%%function for measuring error

 err=@(x,y) DiscrepancyMeasure(x,y,TOL,CHECKTYPES);
 
  
Pf=srand(3,3,2,4);
Qf=srand(size(Pf));
Af=srand(3,2)*1i;  
Bf=srand(3,2);
Cf=srand(3,2);
Sf=srand(3);
Vrf=1:(numel(Pf)-3);
Vcf=Vrf.';
three=single(3);

%%Representations of the above as ndSparse objects 
P=ndSparse(Pf(:),size(Pf));   % =Pf
Q=ndSparse(Qf(:),size(Qf));   % =Qf
A=ndSparse(Af);   % =Af
B=ndSparse(Bf);   % =Bf
C=ndSparse(Cf);   % =Cf
S=ndSparse(Sf);   % =Sf
Vc=ndSparse(Vcf);   % =Vcf
Vr=ndSparse(Vrf);   % =Vrf

%%%%%%%%%%%%%%%%%%%%%%%%%%%TESTS%%%%%%%%%%%%%%%%%%%%%%%%%%


%%Test of full()
Error(1)=     err( full(P), Pf  );  
Error(end+1)= err( full(Q), Qf  );  
Error(end+1)= err( full(S), Sf );  

%%Test of sparse2d, sparse

  %Error(end+1)= ~isequal( sparse(A), A ); %obsolete test of sparse(ndSparse)

Error(end+1)= ~isequal( sparse2d(A), reshape(sparse(Af),[],size(Af,ndims(Af)) ) );
Error(end+1)= ~isequal( sparse(A), sparse(Af) );

   %%No sense in proceding if the tests so far didn't pass - the error
   %%calculations rely on the functionality of full@ndSparse()
    if max(Error)>TOL,
        Error, %#ok<NOPRT>
        error 'Something wrong with sparse() and full() methods'; 
    end

%Test of logical, double

Error(end+1)=err( class(logical(P)), class(logical(Pf)));
Error(end+1)=err( class(double(P)), class(double(Pf)));


%%Test isnumeric, islogical, isempty, issparse, isfloat, isreal,isinf, 
%      isnan, isfinite, isequal, isequalwithequalnans            
Error(end+1)= err( isnumeric(P),1);
Error(end+1)= err( ~isnumeric(P>0),1);
Error(end+1)= err( isfloat(P),1);
Error(end+1)= err( ~isfloat(P>0),1);
Error(end+1)= err( ~islogical(P),1);
Error(end+1)= err( islogical(P>0),1);           
Error(end+1)= err( ~isempty(P),1);   
Error(end+1)= err( isempty(ndSparse([])),1);
Error(end+1)= err( issparse(P),1);
Error(end+1)= err( issparse(A), 1 );
Error(end+1)= err( ~isreal(A), 1 );
Error(end+1)= err( isreal(B), 1 );

     Z=P; Z(1)=nan; Z(2)=inf;
     Zf=full(Z);
     
Error(end+1)= err( isnan(Z), isnan(Zf) );
Error(end+1)= err( isinf(Z), isinf(Zf) );     
Error(end+1)= err( isfinite(Z), isfinite(Zf) );      
Error(end+1)= err( ~isequal(Z,Zf), 1 ); 
Error(end+1)= err( ~isequal(Zf,Z), 1 ); 
Error(end+1)= err( ~isequal(Z,Zf,Z), 1 ); 
Error(end+1)= err( ~isequal(Zf,Z,Z), 1 ); 
Error(end+1)= err( isequalwithequalnans(Z,Zf), 1 );  %#ok<*DISEQN>
Error(end+1)= err( isequalwithequalnans(Zf,Z), 1 );
Error(end+1)= err( isequalwithequalnans(Z,Zf,Z), 1 ); 
Error(end+1)= err( isequalwithequalnans(Zf,Z,Zf), 1 ); 


%Test of real, imag, conj, abs,sqrt

     Z=ndSparse(Bf+Af); %complex result
     Zf=full(Z);
     
Error(end+1)= err( real(Z), real(Zf) );
Error(end+1)= err( imag(Z), imag(Zf) );     
Error(end+1)= err( conj(Z), conj(Zf) );   
Error(end+1)= err( abs(Z), abs(Zf) );
Error(end+1)= err( sqrt(Z), sqrt(Zf) );
            


%%Test of size

Error(end+1)=err( size(P), size(Pf) );
Error(end+1)=err( size(P,2), size(Pf,2) );

       [mm,nn]  =size(Q);
       [mmm,nnn]=size(Qf);
       
Error(end+1)=err( [mm,nn] , [mmm,nnn] );

%%Test of reshape

Error(end+1)= err( reshape(P,size(P,1),[]) , reshape(Pf,size(Pf,1),[]));

%%Test of permute, ipermute

       ord=randperm(ndims(Pf));
       Z=permute(P,ord);
       Zf=permute(Pf,ord);
       
Error(end+1)= err( Z , Zf );
Error(end+1)= err( ipermute(Z,ord) ,  Pf );


%%Test of transpose, ctranspose
Error(end+1)= err( A.' ,  Af.'   );
Error(end+1)= err( A'  ,  Af'   );


%%Test of uplus, uminus

Error(end+1)= err( +P  , +Pf );
Error(end+1)= err( -Q  , -Qf );

%%Test of plus, minus

Error(end+1)= err( P+Q  , Pf+Qf );
Error(end+1)= err( P-Q ,  Pf-Qf );


%%Test of inv

Error(end+1)= err( inv(S)  , inv(Sf) );

%%Test of find

            [II,JJ,KK]=find(P);
            [IIf,JJf,KKf]=find(Pf);      
            
Error(end+1)= err( II  , IIf );
Error(end+1)= err( JJ  , JJf );
Error(end+1)= err( KK  , KKf );

%%Test of mtimes

Error(end+1)= err( A*three , Af*three  );  %scalar with ndSparse
Error(end+1)= err( three*P , three*Pf  );

Error(end+1)= err( S*A , Sf*Af  ); % 2 ndSparses

   x=Af(end,:).';
   y=Af(:,end).';

Error(end+1)= err( A*[x,x] , Af*[x,x]  ); %pre-mult with columnized data
Error(end+1)= err( [y;y]*A , [y;y]*Af  ); %post-mult with columnized data

Error(end+1)= err( S*Af , Sf*Af  ); %mixed op



%%Test of mldivide
Error(end+1)= err( three\P , three\Pf  );  %scalar with ndSparse
Error(end+1)= err( three\A , three\Af  );
Error(end+1)= err( B\A , Bf\Af  ); % 2 ndSparses
Error(end+1)= err( B\Af , Bf\Af  );


%%Test of mrdivide
  
    Bt=B.'; Btsp=Bf.';

Error(end+1)= err( Bt/three , Btsp/three  );%scalar with ndSparse

Error(end+1)= err( P/three , Pf/three  );%scalar with ndSparse
Error(end+1)= err( A.'/Bt , Af.'/Btsp ); % 2 ndSparses
Error(end+1)= err( Af.'/Bt , Af.'/Btsp );


%%Test of times

Error(end+1)= err( P.*three , Pf.*three  );  %scalar with ndSparse
Error(end+1)= err( three.*P , three.*Pf  );

Error(end+1)= err( P.*Q , Pf.*Qf  ); % 2 ndSparses
Error(end+1)= err( Pf.*Q , Pf.*Qf  ); 

%%Test of rdivide

Error(end+1)= err( Q./three , Qf./three  );  %scalar with ndSparse
Error(end+1)= err( three./Q , three./Qf  );

Error(end+1)= err( P./Q , Pf./Qf  ); %2 ndSparses
Error(end+1)= err( Q./P , Qf./Pf  ); %2 ndSparses
Error(end+1)= err( Qf./P , Qf./Pf  ); %mixed

%%Test of ldivide

Error(end+1)= err( Q.\three , Qf.\three  );  %scalar with ndSparse
Error(end+1)= err( three.\Q , three.\Qf  );

Error(end+1)= err( Q.\P , Qf.\Pf  ); %2 ndSparses
Error(end+1)= err( P.\Q , Pf.\Qf  ); %2 ndSparses
Error(end+1)= err( Pf.\Q , Pf.\Qf  );

%%Test of power, mpower
Error(end+1)= err( P.^three ,  Pf.^three   );
Error(end+1)= err( three.^P ,  three.^Pf   );
Error(end+1)= err( Q.^P ,  Qf.^Pf   );
Error(end+1)= err( Qf.^P ,  Qf.^Pf   );

%%Test of mpower
Error(end+1)= err( S^three  ,  Sf^three   );


%%Test of relops

Error(end+1)= err( P>P ,  Pf>Pf   );
Error(end+1)= err( P>Pf ,  Pf>Pf   );
Error(end+1)= err( P>three/6  , Pf>three/6   );

Error(end+1)= err( P>=P ,  Pf>=Pf   );
Error(end+1)= err( P>=Pf ,  Pf>=Pf   );
Error(end+1)= err( P>=three/6  , Pf>=three/6   );

Error(end+1)= err( P<P ,  Pf<Pf   );
Error(end+1)= err( P<Pf ,  Pf<Pf   );
Error(end+1)= err( P<three/6  , Pf<three/6   );

Error(end+1)= err( P<=P ,  Pf<=Pf   );
Error(end+1)= err( P<=Pf ,  Pf<=Pf   );
Error(end+1)= err( P<=three/6  , Pf<=three/6   );


Error(end+1)= err( P==P ,  Pf==Pf   );
Error(end+1)= err( P==Pf ,  Pf==Pf   );
Error(end+1)= err( P==three/6  , Pf==three/6   );

Error(end+1)= err( P~=P ,  Pf~=Pf   );
Error(end+1)= err( P~=Pf ,  Pf~=Pf   );
Error(end+1)= err( P~=three/6  , Pf~=three/6   );

%%Test of logical ops

Error(end+1)= err( P&P ,  Pf&Pf   );
Error(end+1)= err( P&Pf ,  Pf&Pf   );
Error(end+1)= err( P|three/6  , Pf|three/6   );

Error(end+1)= err( ~P ,   ~Pf   );
Error(end+1)= err( ~Q  ,  ~Qf  );


%%Test of sum

Error(end+1)= err( sum(P,1)  , sum(Pf,1) );
Error(end+1)= err( sum(Q,2)  , sum(Qf,2) );
Error(end+1)= err( sum(P)    , sum(Pf) );
Error(end+1)= err( sum(P,'native')    , sum(Pf,'native') );
Error(end+1)= err( sum(P,3,'double')    , sum(Pf,3,'double') );
Error(end+1)= err( sum(P,4)    , sum(Pf,4) );
Error(end+1)= err( sum(P,5)    , sum(Pf,5) );

%%Test of cat, horzcat, vertcat

Error(end+1)= err([P,P] ,[Pf,Pf]);
Error(end+1)= err([Pf,P] ,[Pf,Pf]);
Error(end+1)= err([P,Pf] ,[Pf,Pf]);

Error(end+1)= err([P;P] ,[Pf;Pf]);
Error(end+1)= err([Pf;P] ,[Pf;Pf]);
Error(end+1)= err([P;Pf] ,[Pf;Pf]);

Error(end+1)= err([P,P] ,[Pf,Pf]);
Error(end+1)= err([Pf,P] ,[Pf,Pf]);
Error(end+1)= err([P,Pf] ,[Pf,Pf]);

Error(end+1)= err(cat(3,P,P,P) ,cat(3,Pf,Pf,Pf));
Error(end+1)= err(cat(3,Pf,P) ,cat(3,Pf,Pf));
Error(end+1)= err(cat(3,P,Pf) ,cat(3,Pf,Pf));

Error(end+1)= err(cat(4,P,Pf) ,cat(4,Pf,Pf));
Error(end+1)= err(cat(5,P,Pf) ,cat(5,Pf,Pf));


%%Test of spfun

 f=@(x) cos(x).^2;
 
Error(end+1)= err( spfun(f,P) ,    reshape(full(spfun(f,Pf(:))) ,size(Pf)) );




%%Test of subsindex

   idx=1:3;
   Zf=rand(3);

Error(end+1)= err( Zf(ndSparse(idx)) , Zf(idx) );   
   

%%Test of subsref


Error(end+1)= err( P(P<.5) ,  Pf(Pf<.5)   ); %logical indexing
Error(end+1)= err( Pf(P<.5) ,  Pf(Pf<.5)   ); 
Error(end+1)= err( P(Pf<.5) ,  Pf(Pf<.5)   );

Error(end+1)= err( P(:) ,  Pf(:)   );  %linear indexing
Error(end+1)= err( P(1) ,  Pf(1)   );  
Error(end+1)= err( P(1:4) ,  Pf(1:4)   );
Error(end+1)= err( P((1:4).') ,  Pf((1:4).')   );
Error(end+1)= err( Vc(1:3) ,  Vcf(1:3)   );
Error(end+1)= err( Vr(1:3) ,  Vrf(1:3)   );

Error(end+1)= err( P(Vc(1:3)) ,  Pf(Vcf(1:3))   ); %indexing vectors test different shaping rules
Error(end+1)= err( Pf(Vc(1:3)) ,  Pf(Vcf(1:3))   );
Error(end+1)= err( P(Vcf(1:3)) ,  Pf(Vcf(1:3))   );
Error(end+1)= err( Pf(Vcf(1:3)) ,  Pf(Vcf(1:3))   );


Error(end+1)= err( P(2,1,2,2) ,  Pf(2,1,2,2)   ); %subscript indexing
Error(end+1)= err( P(2,:,2,2) ,  Pf(2,:,2,2)   );
Error(end+1)= err( P(:,1,:,2) ,  Pf(:,1,:,2)   ); 
Error(end+1)= err( P(:,2,2,2) ,  Pf(:,2,2,2)   );
Error(end+1)= err( P(1,2,2,:) ,  Pf(1,2,2,:)   );
Error(end+1)= err( P(:,:,:,2) ,  Pf(:,:,:,2)   );

Error(end+1)= err( P(2,1,2) ,  Pf(2,1,2)   ); %truncated subscript indexing
Error(end+1)= err( P(2,:,2) ,  Pf(2,:,2)   );
Error(end+1)= err( P(:,1,:) ,  Pf(:,1,:)   ); 
Error(end+1)= err( P(:,2,2) ,  Pf(:,2,2)   );

            lidx=logical([1 0]);
            
Error(end+1)= err( P(lidx,lidx,1:2) ,  Pf(lidx,lidx,1:2)   ); %combine all types of indexing
Error(end+1)= err( P(lidx,lidx,:) ,  Pf(lidx,lidx,:)   ); 
Error(end+1)= err( P(lidx,:,lidx) ,  Pf(lidx,:,lidx)   ); 


%%Test of subsasgn

        Z=P; 
        Zf=Pf;

        Z(Z<.5)=6;   Zf(Zf<.5)=6;
       
        
Error(end+1)= err(  Z, Zf); %logical indexing

        Z(1)=7;   Zf(1)=7;
       
        
Error(end+1)= err(  Z, Zf); %scalar indexing



        Z(Zf<.5)=13; Zf(Zf<.5)=13;
       
        
Error(end+1)= err( Z ,  Zf   );

           Z(:)=999; Zf(:)=999;
           
Error(end+1)= err(  Z, Zf    ); 

           Z(1:4)=88;  Zf(1:4)=88;  

Error(end+1)= err( Z ,  Zf   );

          Z((1:4).')=77;  Zf((1:4).')=77; 

Error(end+1)= err(Z,Zf);
         
          d=rand;
          Z(2,1,2,2)=d;   Zf(2,1,2,2)=d;
          
Error(end+1)= err(Z,Zf); %subscript indexing

           d=rand;         
          Z(2,:,2,2)=d;  Zf(2,:,2,2)=d;

Error(end+1)= err(Z,Zf);

          d=rand;
          Z(:,2,2,2)=d;  Zf(:,2,2,2)=d; 

Error(end+1)= err(Z,Zf);

          d=rand;         
          Z(:,:,:,2)=d;  Zf(:,:,:,2)=d;  
          


Error(end+1)= err(Z,Zf);

           mm=size(Z,1)+3;
           nn=size(Z,2)+3;

          Z(mm,nn,:,:)=5;  Zf(mm,nn,:,:)=5; %matrix expansion test

          
Error(end+1)= err(Z,Zf);

          Z(:,1:2:end,:,:)=[];  Zf(:,1:2:end,:,:)=[]; %null assignment test          
          
Error(end+1)= err(Z,Zf);

           Z=P; Zf=Pf;
           Z(:,mod(1:2:end,2)==1,:,:)=[];  %null assignment test logical indexing
           Zf(:,mod(1:2:end,2)==1,:,:)=[];
           
           
Error(end+1)= err(Z,Zf);          
          
      Z=P; Zf=Pf;
      Z(:,:,:,:,1)=0;  Zf(:,:,:,:,1)=0; %matrix expansion with n-1 colons          
      Z(:,:,:,:,1)=[];  Zf(:,:,:,:,1)=[]; %null assignment with n-1 colons 
      
      
Error(end+1)= err(Z,Zf); 

           Z=P; Zf=Pf;

           d=rand(size(Zf(:,:,1)));
           
           Z(:,:,1)=d;  Zf(:,:,1)=d; %nonscalar assignment
           
Error(end+1)= err( Z , Zf   ); 

            idx={':'};
            
            Z=P; Zf=Pf;
            d=rand(size(Zf(idx{:})));
            
            Z(idx{:})=d;  Zf(idx{:})=d;%nonscalar assignment
            
                     
Error(end+1)= err(Z,Zf);            

            lidx=logical([1 0]);
            
            idx={lidx,lidx,1:2};
            
            Z=P; Zf=Pf;
            d=rand(size(Zf(idx{:})));
            
            Z(idx{:})=d;  Zf(idx{:})=d;%combine all types of indexing
            
            
Error(end+1)= err( Z , Zf   ); 






%%Test of all(), any(), mean(), max/min

      Z=P>.5; Zf=Pf>.5;
      
      Z(:,1,1,1)=0; Zf(:,1,1,1)=0;
     
      Args={ {},{1},{2},{3},{4},{5} };
      for jj=1:length(Args)
          
          args=Args{jj};
          
Error(end+1)= err( all(Z,args{:})  , all(Zf,args{:}) );

Error(end+1)= err( any(Z,args{:})  , any(Zf,args{:}) );

Error(end+1)= err( mean(Z,args{:})  , mean(Zf,args{:}) );

Error(end+1)= err( all(P,args{:})  , all(Pf,args{:}) );

Error(end+1)= err( any(P,args{:})  , any(Pf,args{:}) );

Error(end+1)= err( mean(P,args{:})  , mean(Pf,args{:}) );
      end
 

        Args={ {1},{2},{3},{4},{5} };
 
      for jj=1:length(Args)
          
          args=Args{jj};        
        
          [Z,idx]=max(P,[],args{:});
          [Zf,idxf]=max(Pf,[],args{:});
          
Error(end+1)= err( Z ,Zf );       
Error(end+1)= err( idx ,idxf );    

      end
      
Error(end+1)= err( max(P,Q) , max(Pf,Qf) );      
Error(end+1)= err( max(P,Qf) , max(Pf,Qf) );      
Error(end+1)= err( max(Pf,Q) , max(Pf,Qf) );   

 
      for jj=1:length(Args)
          
          args=Args{jj};        
        
          [Z,idx]=min(P,[],args{:});
          [Zf,idxf]=min(Pf,[],args{:});
          
Error(end+1)= err( Z ,Zf );       
Error(end+1)= err( idx ,idxf );    

      end
      
Error(end+1)= err( min(P,Q) , min(Pf,Qf) );      
Error(end+1)= err( min(P,Qf) , min(Pf,Qf) );      
Error(end+1)= err( min(Pf,Q) , min(Pf,Qf) );   


       %%%%%%2D cases%%%%%

      Z=B>.5; Zf=Bf>.5;
      
      Z(:,1,1,1)=0; Zf(:,1,1,1)=0;
     
      Args={ {},{1},{2},{3},{4},{5} };
      for jj=1:length(Args)
          
          args=Args{jj};
          
Error(end+1)= err( all(Z,args{:})  , all(Zf,args{:}) );

Error(end+1)= err( any(Z,args{:})  , any(Zf,args{:}) );

Error(end+1)= err( mean(Z,args{:})  , mean(Zf,args{:}) );


      end
 

        Args={ {1},{2},{3},{4},{5} };
 
      for jj=1:length(Args)
          
          args=Args{jj};        
        
          [Z,idx]=max(B,[],args{:});
          [Zf,idxf]=max(Bf,[],args{:});
          
Error(end+1)= err( Z ,Zf );       
Error(end+1)= err( idx ,idxf );    

      end
      
Error(end+1)= err( max(B,C) , max(Bf,Cf) );      
Error(end+1)= err( max(B,Cf) , max(Bf,Cf) );      
Error(end+1)= err( max(Bf,C) , max(Bf,Cf) );   

 
      for jj=1:length(Args)
          
          args=Args{jj};        
        
          [Z,idx]=min(B,[],args{:});
          [Zf,idxf]=min(Bf,[],args{:});
          
Error(end+1)= err( Z ,Zf );       
Error(end+1)= err( idx ,idxf );    

      end
      
Error(end+1)= err( min(B,C) , min(Bf,Cf) );      
Error(end+1)= err( min(B,Cf) , min(Bf,Cf) );      
Error(end+1)= err( min(Bf,C) , min(Bf,Cf) );   
       
       
%Test of numel
Error(end+1)= err( numel(P) , numel(Pf) );

%Test of repmat, circshift

       N=ndims(Pf)+1;
       z=ones(1,N);

       
       for kk=0:N

           idx=nchoosek(1:N,kk);
           
           for jj=1:size(idx,1)
               
             arg2=z;
             arg2(idx(jj,:))=2;
  
Error(end+1)= err( repmat(P,arg2) , repmat(Pf,arg2) );
Error(end+1)= err( repmat(P,arg2-1) , repmat(Pf,arg2-1) );
Error(end+1)= err( circshift(P,arg2) , circshift(Pf,arg2) );
Error(end+1)= err( circshift(P,arg2-1) , circshift(Pf,arg2-1) );             
Error(end+1)= err( circshift(A,arg2) , circshift(Af,arg2) );
Error(end+1)= err( circshift(A,arg2-1) , circshift(Af,arg2-1) );             

           end
       end

              Z=P*0; Zf=full(Z);
       
Error(end+1)= err( repmat(Z,arg2) , repmat(Zf,arg2) );
Error(end+1)= err( circshift(Z,arg2) , circshift(Zf,arg2) );
Error(end+1)= err( repmat(Z,arg2-1) , repmat(Zf,arg2-1) );
Error(end+1)= err( circshift(Z,arg2-1) , circshift(Zf,arg2-1) );
       
       
%Test of squeeze, shiftdim      

     N=ndims(Pf);
     z=ones(1,2*N); z(2:2:end)=size(P);
     Z=reshape(P,[1 1 z]); %add some singleton dimensions
     Zf=full(Z);
     
Error(end+1)= err( squeeze(Z) , squeeze(Zf) );     

    [Z,n]=shiftdim(Z); [Zf,nf]=shiftdim(Zf);

Error(end+1)= err( Z , Zf);  
Error(end+1)= err( n , nf);  

Error(end+1)= err( shiftdim(Z,3) , shiftdim(Zf,3));
Error(end+1)= err( shiftdim(Z,-3) , shiftdim(Zf,-3));

%Test of bsxfun


     funcs={@plus,@minus,@times,@rdivide,@ldivide,@power,...
            @max,@min,@rem,@mod,@atan2,@hypot,@(a,b) a.^2-b};

     lfuncs={@eq,@ne,@lt,@le,@gt,@ge,@and,@or,@xor};
      
     HH=P;
     for qq=1:2 
         
     
     
     %test non-logical funcs
     H=HH; Hf=full(H);
     Z=mean(H,1); Z=mean(Z,3);   Zf=full(Z);
     
     
     for kk=1:2
        for ii=1:length(funcs)
         
          fun=funcs{ii};
          
Error(end+1)= err( bsxfun(fun,Z,H) ,     bsxfun(fun,Zf,Hf)  );     %#ok<*AGROW>
Error(end+1)= err( bsxfun(fun,Z,Hf) ,     bsxfun(fun,Zf,Hf)  );           
Error(end+1)= err( bsxfun(fun,Zf,H) ,     bsxfun(fun,Zf,Hf)  );    
Error(end+1)= err( bsxfun(fun,three,H) ,     bsxfun(fun,three,Hf)  );  
Error(end+1)= err( bsxfun(fun,Zf,three) ,     bsxfun(fun,Zf,three)  ); 

       end
      
      Z=mean(H,2); Zf=full(Z);
      
      end
      
      %test logical lfuncs      
      H=(HH<0.5); Hf=full(H);
      Z=mean(H,1); Z=mean(Z,3); 
      Z=(Z>0.5);  Zf=full(Z);
     
     
      onebit=(rand>=0.5);
      
      for kk=1:2
       for ii=1:length(lfuncs)
         
          fun=lfuncs{ii};
          
Error(end+1)= err( bsxfun(fun,Z,H) ,     bsxfun(fun,Zf,Hf)  );    
Error(end+1)= err( bsxfun(fun,Z,Hf) ,     bsxfun(fun,Zf,Hf)  );           
Error(end+1)= err( bsxfun(fun,Zf,H) ,     bsxfun(fun,Zf,Hf)  );    
Error(end+1)= err( bsxfun(fun,onebit,H) ,     bsxfun(fun,onebit,Hf)  );  
Error(end+1)= err( bsxfun(fun,Zf,onebit) ,     bsxfun(fun,Zf,onebit)  ); 

       end
      
      Z=mean(H,2); 
      Z=(Z>0.5); Zf=full(Z);
      
      
      end
     
      HH=repmat(P,[1,1,1,0]);
     
     end
 
%Test of implicit expansion


     funcs={@plus,@minus,@times,@rdivide,@ldivide,@power,...
            @max,@min,@rem,@mod,@atan2,@hypot,@(a,b) a.^2-b};

     lfuncs={@eq,@ne,@lt,@le,@gt,@ge,@and,@or,@xor};
      
     HH=P;
     for qq=1:2 
         
     
     
     %test non-logical funcs
     H=HH; Hf=full(H);
     Z=mean(H,1); Z=mean(Z,3);   Zf=full(Z);
     
     
     for kk=1:2
        for ii=1:length(funcs)
         
          fun=funcs{ii};
          
Error(end+1)= err( feval(fun,Z,H) ,     bsxfun(fun,Zf,Hf)  );     %#ok<*AGROW>
Error(end+1)= err( feval(fun,Z,Hf) ,     bsxfun(fun,Zf,Hf)  );           
Error(end+1)= err( feval(fun,Zf,H) ,     bsxfun(fun,Zf,Hf)  );    
Error(end+1)= err( feval(fun,three,H) ,     bsxfun(fun,three,Hf)  );  
Error(end+1)= err( feval(fun,Zf,three) ,     bsxfun(fun,Zf,three)  ); 

       end
      
      Z=mean(H,2); Zf=full(Z);
      
      end
      
      %test logical lfuncs      
      H=(HH<0.5); Hf=full(H);
      Z=mean(H,1); Z=mean(Z,3); 
      Z=(Z>0.5);  Zf=full(Z);
     
     
      onebit=(rand>=0.5);
      
      for kk=1:2
       for ii=1:length(lfuncs)
         
          fun=lfuncs{ii};
          
Error(end+1)= err( feval(fun,Z,H) ,     bsxfun(fun,Zf,Hf)  );    
Error(end+1)= err( feval(fun,Z,Hf) ,     bsxfun(fun,Zf,Hf)  );           
Error(end+1)= err( feval(fun,Zf,H) ,     bsxfun(fun,Zf,Hf)  );    
Error(end+1)= err( feval(fun,onebit,H) ,     bsxfun(fun,onebit,Hf)  );  
Error(end+1)= err( feval(fun,Zf,onebit) ,     bsxfun(fun,Zf,onebit)  ); 

       end
      
      Z=mean(H,2); 
      Z=(Z>0.5); Zf=full(Z);
      
      
      end
     
      HH=repmat(P,[1,1,1,0]);
     
     end     
     
     
%Test of ndSparse.build, nzmax, nonzeros   
      
      nzm=10;
      Z=ndSparse.build([1 1 1; 2 2 2; 3 3 3; 4 4 4], 1:4,[4 4 6],nzm);
      Zf=zeros(4,4,6); for ii=1:4, Zf(ii,ii,ii)= ii; end

Error(end+1)= err( Z , Zf  );    
Error(end+1)= err( nzmax(Z) , nzm );
Error(end+1)= err( nonzeros(Z) , nonzeros(Zf) );

      nzm=10;
      Z=ndSparse.build([1 1 1; 2 2 2; 3 3 3; 4 4 4], 5,[4 4 6],nzm);
      Zf=zeros(4,4,6); for ii=1:4, Zf(ii,ii,ii)= 5; end

Error(end+1)= err( Z , Zf  );    
Error(end+1)= err( nzmax(Z) , nzm );
Error(end+1)= err( nonzeros(Z) , nonzeros(Zf) );


       Z=ndSparse.build([4 4 6]);
       Zf=zeros(4,4,6); 

Error(end+1)= err( Z , Zf  );   
Error(end+1)= err( ndSparse.build([],[],[4 4 6]) , Zf  );    
Error(end+1)= err( nzmax(Z) , 1 );
Error(end+1)= err( nonzeros(Z) , nonzeros(Zf) );       
       

%Test of ndSparse.spalloc

      nzm=10;
      Z=ndSparse.spalloc([3,5,3],nzm);
      Zf=full(Z);
      
Error(end+1)= err( nzmax(Z) , nzm );  
Error(end+1)= err( Z , Zf );  

%Test of ndSparse.accumarray

         args={[1 1 1; 2 2 2; 3 3 3;2 2 2], [5,4,3,2]};
         Z=ndSparse.accumarray(args{:});
         Zf=accumarray(args{:});
         
Error(end+1)= err( Z , Zf );  

         args={[1 1 1; 2 2 2; 3 3 3;2 2 2], [5,4,3,2],[4 4 4]};
         Z=ndSparse.accumarray(args{:});
         Zf=accumarray(args{:});

Error(end+1)= err( Z , Zf );          

         args={[1;2;3;2], [5,4,3,2]};
         Z=ndSparse.accumarray(args{:});
         Zf=accumarray(args{:});

Error(end+1)= err( Z , Zf );   

         args={[1;2;3;2], [5,4,3,2],[4,1]};
         Z=ndSparse.accumarray(args{:});
         Zf=accumarray(args{:});

Error(end+1)= err( Z , Zf ); 

    
         args={[1;2;3;2], [5,4,3,2],[4,1],@prod};
         Z=ndSparse.accumarray(args{:});
         Zf=accumarray(args{:});

Error(end+1)= err( Z , Zf ); 

%Test of spones

          Z=ndSparse.build([1,1,1], 5,[3,5,3]);
          Zf=full(Z); Zf(~~Zf)=1;
       
Error(end+1)= err( spones(Z) , Zf );  
      
%Test of length

Error(end+1)= err( length(P) , length(Pf) );  

%Test of triu,triul

Error(end+1)= err( triu(A) , triu(Af) );
Error(end+1)= err( triu(A,1) , triu(Af,1) );
Error(end+1)= err( triu(A,-1) , triu(Af,-1) );

Error(end+1)= err( tril(A) , tril(Af) );
Error(end+1)= err( tril(A,1) , tril(Af,1) );
Error(end+1)= err( tril(A,-1) , tril(Af,-1) );

Error(end+1)= err( triu(S) , triu(Sf) );
Error(end+1)= err( triu(S,1) , triu(Sf,1) );
Error(end+1)= err( triu(S,-1) , triu(Sf,-1) );

Error(end+1)= err( tril(S) , tril(Sf) );
Error(end+1)= err( tril(S,1) , tril(Sf,1) );
Error(end+1)= err( tril(S,-1) , tril(Sf,-1) );

%Test of flipdim

Error(end+1)= err( flipdim(P,1) , flipdim(Pf,1) );
Error(end+1)= err( flipdim(P,2) , flipdim(Pf,2) );
Error(end+1)= err( flipdim(P,3) , flipdim(Pf,3) );


%Test of  fliplr

   Z=P(:,:,:,1); Zf=full(Z);   for ii=1:size(Zf,3), Zf(:,:,ii)=fliplr( Zf(:,:,ii) ); end

Error(end+1)= err( fliplr(Z) ,  Zf );

%Test of flipud 

   Z=P(:,:,:,1);  Zf=full(Z);    for ii=1:size(Zf,3), Zf(:,:,ii)=flipud( Zf(:,:,ii) ); end
   
Error(end+1)= err( flipud(Z) ,  Zf );

%Test of rot90

  Z=P(:,:,:,1);  Zf=full(Z);    for ii=1:size(Zf,3), Zf(:,:,ii)=rot90( Zf(:,:,ii) ); end
   
Error(end+1)= err( rot90(Z) ,  Zf );

             for kk=-5:5

   Z=P(:,:,:,1);  Zf=full(Z);    for ii=1:size(Zf,3), Zf(:,:,ii)=rot90( Zf(:,:,ii) ,kk); end
   
Error(end+1)= err( rot90(Z,kk) , Zf);

             end

%Test of convn

            args={{},{'full'},{'same'},{'valid'}};
            Zs=sparse([0;0;1;2]); Zf=full(Zs);
            
            for ii=1:length(args)
                
Error(end+1)= err( convn(P,P,args{ii}{:}) , convn(Pf,Pf,args{ii}{:}) );
Error(end+1)= err( convn(P,A,args{ii}{:}) , convn(Pf,Af,args{ii}{:}) );

Error(end+1)= err( convn(Pf,P,args{ii}{:}) , convn(Pf,Pf,args{ii}{:}) );
Error(end+1)= err( convn(Pf,A,args{ii}{:}) , convn(Pf,Af,args{ii}{:}) );

Error(end+1)= err( convn(P,Pf,args{ii}{:}) , convn(Pf,Pf,args{ii}{:}) );
Error(end+1)= err( convn(P,Af,args{ii}{:}) , convn(Pf,Af,args{ii}{:}) );
             

Error(end+1)= err( convn(P,Zs,args{ii}{:}) , convn(Pf,Zf,args{ii}{:}) );
Error(end+1)= err( convn(Zs,P,args{ii}{:}) , convn(Zf,Pf,args{ii}{:}) );

Error(end+1)= err( convn(P,Zf,args{ii}{:}) , convn(Pf,Zf,args{ii}{:}) );
Error(end+1)= err( convn(Zf,P,args{ii}{:}) , convn(Zf,Pf,args{ii}{:}) );

            end

            Hf={rand(2,1),rand(1,3)};
            H=Hf;
            H{1}=ndSparse(H{1}); 
            
            Zf=diag([0 0 0 0 0 1 1]);
            Zs=sparse(Zf); 
            
           for ii=1:length(args)
                
Error(end+1)= err( convn(H{:},A,args{ii}{:}) , conv2(Hf{:},Af,args{ii}{:}) );

Error(end+1)= err( convn(H{:},Zs,args{ii}{:}) , conv2(Hf{:},Zf,args{ii}{:}) );


Error(end+1)= err( convn(H{:},Zf,args{ii}{:}) , conv2(Hf{:},Zf,args{ii}{:}) );

           end


           
           
%%Test of allml(), anyml(), meanml(), maxml/minml, summl, catml

      Z=P>.5; Zf=Pf>.5;
      
      Z(:,1,1,1)=0; Zf(:,1,1,1)=0;
     
      Args={ {},{1},{2},{3},{4},{5} };
      for jj=1:length(Args)
          
          args=Args{jj};
          
Error(end+1)= err( allml(Z,args{:})  , all(Zf,args{:}) );

Error(end+1)= err( anyml(Z,args{:})  , any(Zf,args{:}) );

Error(end+1)= err( meanml(Z,args{:})  , mean(Zf,args{:}) );

Error(end+1)= err( allml(P,args{:})  , all(Pf,args{:}) );

Error(end+1)= err( anyml(P,args{:})  , any(Pf,args{:}) );

Error(end+1)= err( meanml(P,args{:})  , mean(Pf,args{:}) );
      end
 

        Args={ {1},{2},{3},{4},{5} };
 
      for jj=1:length(Args)
          
          args=Args{jj};        
        
          [Z,idx]=maxml(P,[],args{:});
          [Zf,idxf]=max(Pf,[],args{:});
          
Error(end+1)= err( Z ,Zf );       
Error(end+1)= err( idx ,idxf );    

      end
      
Error(end+1)= err( maxml(P,Q) , max(Pf,Qf) );      
Error(end+1)= err( maxml(P,Qf) , max(Pf,Qf) );      
Error(end+1)= err( maxml(Pf,Q) , max(Pf,Qf) );   

 
      for jj=1:length(Args)
          
          args=Args{jj};        
        
          [Z,idx]=minml(P,[],args{:});
          [Zf,idxf]=min(Pf,[],args{:});
          
Error(end+1)= err( Z ,Zf );       
Error(end+1)= err( idx ,idxf );    

      end
      
Error(end+1)= err( minml(P,Q) , min(Pf,Qf) );      
Error(end+1)= err( minml(P,Qf) , min(Pf,Qf) );      
Error(end+1)= err( minml(Pf,Q) , min(Pf,Qf) );   


       %%%%%%2D cases%%%%%

      Z=B>.5; Zf=Bf>.5;
      
      Z(:,1,1,1)=0; Zf(:,1,1,1)=0;
     
      Args={ {},{1},{2},{3},{4},{5} };
      for jj=1:length(Args)
          
          args=Args{jj};
          
Error(end+1)= err( allml(Z,args{:})  , all(Zf,args{:}) );

Error(end+1)= err( anyml(Z,args{:})  , any(Zf,args{:}) );

Error(end+1)= err( meanml(Z,args{:})  , mean(Zf,args{:}) );


      end
 

        Args={ {1},{2},{3},{4},{5} };
 
      for jj=1:length(Args)
          
          args=Args{jj};        
        
          [Z,idx]=maxml(B,[],args{:});
          [Zf,idxf]=max(Bf,[],args{:});
          
Error(end+1)= err( Z ,Zf );       
Error(end+1)= err( idx ,idxf );    

      end
      
Error(end+1)= err( maxml(B,C) , max(Bf,Cf) );      
Error(end+1)= err( maxml(B,Cf) , max(Bf,Cf) );      
Error(end+1)= err( maxml(Bf,C) , max(Bf,Cf) );   

 
      for jj=1:length(Args)
          
          args=Args{jj};        
        
          [Z,idx]=minml(B,[],args{:});
          [Zf,idxf]=min(Bf,[],args{:});
          
Error(end+1)= err( Z ,Zf );       
Error(end+1)= err( idx ,idxf );    

      end
      
Error(end+1)= err( minml(B,C) , min(Bf,Cf) );      
Error(end+1)= err( minml(B,Cf) , min(Bf,Cf) );      
Error(end+1)= err( minml(Bf,C) , min(Bf,Cf) );  

Error(end+1)= err( summl(P,1)  , sum(Pf,1) );
Error(end+1)= err( summl(Q,2)  , sum(Qf,2) );
Error(end+1)= err( summl(P)    , sum(Pf) );
Error(end+1)= err( summl(P,'native')    , sum(Pf,'native') );
Error(end+1)= err( summl(P,3,'double')    , sum(Pf,3,'double') );
Error(end+1)= err( summl(P,4)    , sum(Pf,4) );
Error(end+1)= err( summl(P,5)    , sum(Pf,5) );

Error(end+1)= err(catml(3,P,P,P) ,cat(3,Pf,Pf,Pf));
Error(end+1)= err(catml(3,Pf,P) ,cat(3,Pf,Pf));
Error(end+1)= err(catml(3,P,Pf) ,cat(3,Pf,Pf));

Error(end+1)= err(catml(4,P,Pf) ,cat(4,Pf,Pf));
Error(end+1)= err(catml(5,P,Pf) ,cat(5,Pf,Pf));



       N=ndims(Pf)+1;
       z=ones(1,N);

       
       for kk=0:N

           idx=nchoosek(1:N,kk);
           
           for jj=1:size(idx,1)
               
             arg2=z;
             arg2(idx(jj,:))=2;
 
Error(end+1)= err( circshiftml(P,arg2) , circshift(Pf,arg2) );
Error(end+1)= err( circshiftml(P,arg2-1) , circshift(Pf,arg2-1) );             
Error(end+1)= err( circshiftml(A,arg2) , circshift(Af,arg2) );
Error(end+1)= err( circshiftml(A,arg2-1) , circshift(Af,arg2-1) );             

           end
       end

              Z=P*0; Zf=full(Z);
       
Error(end+1)= err( circshiftml(Z,arg2) , circshift(Zf,arg2) );
Error(end+1)= err( circshiftml(Z,arg2-1) , circshift(Zf,arg2-1) );


%%%%%%%%%%%%%%%%%%%%%%%END OF TESTS%%%%%%%%%%%%%%%%%%%%%%%%%%


MAX_ERROR=max(Error);

if any(~isfinite(Error)), warning 'There were improper Error values'; keyboard; end

disp(['Maximum observed error was   ' num2str(MAX_ERROR) ' percent.'])

if nargout, varargout{1}=Error; end

function errval=DiscrepancyMeasure(X,Y,TOL,CHECKTYPES)


  if Discrepancy(0,Y)
    errval=Discrepancy(X,Y)/Discrepancy(0,Y)*100; %normalize
  else
    errval=Discrepancy(X,Y); 
  end
  
 
  
  isndSparse=@(c) ~isempty(strfind(class(c),'ndSparse'));
  
  if CHECKTYPES
    if ~isndSparse(X),
       warning(['X is not ndSparse class, but rather class ' class(X)]) 
    end
  end

  if errval>TOL || ~isfinite(errval), 
      disp ' '; disp 'Discrepancy detected'
      errval, 
      x=full(X); y=full(Y);
      keyboard;
  end 
  
  
function errval=Discrepancy(X,Y) 
%Primary error measurement function

  fin=@(a) reshape( a(isfinite(a)),[],1);
  nonfin=@(a) reshape(  a(~isfinite(a))  ,[],1); 

  x=full(X); y=full(Y);

  if ( isequal(x,y));
     errval=0; return
  elseif xor(isempty(x),isempty(y))
      errval=1; return
  end
  
  errval= norm( fin(x-y) , inf)+...   
       ~isequalwithequalnans(nonfin(x),nonfin(y))*...
       ~isempty([nonfin(x);nonfin(y)])+ ~isequal(size(X),size(Y)); 



function out=srand(varargin)

   out=rand(varargin{:})-.5;



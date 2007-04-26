function grad = fif2grad(filename);

% FIF2GRAD constructs a gradiometer definition from a Neuromag *.fif file
% The resulting gradiometer definition can be used by Fieldtrip for forward
% and inverse computations.
%
% Use as
%   grad = fif2grad(filename);
%
% See also READ_FCDC_HEADER

% Copyright (C) 2004, Joachim Gross
%
% $Log: fif2grad.m,v $
% Revision 1.2  2007/03/15 12:44:42  roboos
% added try-catch to deal with fif files that do not contain coil information (code by chrhes)
%
% Revision 1.1  2006/08/31 13:32:11  roboos
% moved from fieldtrip to fileio module
%
% Revision 1.1  2004/09/20 13:53:04  roboos
% moved gradiometer definition for Neuromag fif files into separate function
%

% this try-catch construct ensures that missing gradiometer information is
% handeled in a "graceful" way
grad = [];
try
   megmodel('head',[0 0 0],filename);
   [n,s,t]=chaninfo;
   [TY,NA]=chaninfo('type');
   nCoils=sum(TY+1);      % number of coils
   nSensors=length(TY);   % number of sensors
   grad.pnt=zeros(nCoils,3);
   grad.ori=zeros(nCoils,3);
   grad.tra=zeros(nSensors,nCoils);
   grad.unit='cm';
   % define coils
   kCoil=1;
   for k=1:nSensors,
      if (TY(k)==0),  %magnetometer
         grad.pnt(kCoil,:)=100*(t{k}(1:3,4));
         grad.ori(kCoil,:)=t{k}(1:3,3);
         grad.tra(k,kCoil)=1;
         kCoil=kCoil+1;
         grad.label{k}=deblank(s(k,:));
      elseif (TY(k) == 1), %planar gradiometer
         grad.pnt(kCoil,:)=100*(t{k}(1:3,4)-0.008*t{k}(1:3,1)); % multiply with 100 to get cm
         grad.ori(kCoil,:)=t{k}(1:3,3);
         grad.tra(k,kCoil)=1;
         kCoil=kCoil+1;
         grad.pnt(kCoil,:)=100*(t{k}(1:3,4)+0.008*t{k}(1:3,1));
         grad.ori(kCoil,:)=t{k}(1:3,3);
         grad.tra(k,kCoil)=-1;
         kCoil=kCoil+1;
         grad.label{k}=deblank(s(k,:));
      else
         error('unknown sensor type');
      end
   end
catch
   warning(['gradiometer information could not be extracted from file']);
   warning(['returning an empty grad structure']);
   return;
end

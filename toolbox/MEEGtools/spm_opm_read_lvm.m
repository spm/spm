function [lbv] = spm_opm_read_lvm(S)
% Read LVM file
% FORMAT [lbv] = spm_opm_read_lvm(S)
%   S               - input structure
% Optional fields of S:
%   S.filename           - filepath to LVM file             -Default: no Default                      
%   S.headerlength       - integer specifying how many      -Default: 23
%                          lines of file are header
%   S.timeind            - integer specifying which         -Default: 1
%                          column is time variable
%   S.decimalTriggerInds - Indices of trigger Channels      -Default: 74:81
%   S.binaryTriggerInds  - Indices of trigger Channels      -Default: []
%   S.trigThresh         - Value to threshold triggers at   -Default: Auto
%   
% Output: lbv - output Structure
%  Fields of lbv:
%   lbv.B                  - MEG data
%   lbv.Time               - Time variable
%   lbv.decimalTrigs       - Trigger Channels
%   lbv.binaryTrigs        - Trigger Channels
%   lbv.pinout             - pinout of lbv file(coming soon)
%__________________________________________________________________________
% Copyright (C) 2018 Wellcome Trust Centre for Neuroimaging

% Tim Tierney
% $Id: spm_opm_read_lvm.m 7414 2018-09-07 11:00:29Z spm $


%-Set default values
%--------------------------------------------------------------------------
msg = 'filename needs to be provided.';
if ~isfield(S, 'filename'),            error(msg); end
if ~isfield(S, 'headerlength'),        S.headerlength = 23; end
if ~isfield(S, 'timeind'),             S.timeind = 1; end
if ~isfield(S, 'decimalTriggerInds'),  S.decimalTriggerInds = 74:81; end
if ~isfield(S, 'binaryTriggerInds'),   S.binaryTriggerInds = []; end

%-Check for cells 
%--------------------------------------------------------------------------
if iscell(S.filename)
    S.filename=S.filename{1};
end

%-Check for zipped files and read data
%--------------------------------------------------------------------------
[fold,fi,ext]= fileparts(S.filename);
zipped = strmatch(ext,'.zip');

if(zipped)
    cellFile   = unzip(S.filename,fold);
    S.filename = cellFile{1};
end

data = dlmread(S.filename, '\t',S.headerlength,0);


%-Subset  magnetic fields and triggers
%--------------------------------------------------------------------------

chans = 1:size(data,2);
time = data(:,S.timeind);
Bind = setdiff(chans,S.timeind);
B = data(:,Bind);

%-Binarise triggers
%--------------------------------------------------------------------------
decTrigs = data(:,S.decimalTriggerInds);
binTrigs = data(:,S.binaryTriggerInds);

if(isfield(S,'trigThresh'))
    tTrigsDecimal= decTrigs>S.trigThresh;
    tTrigsBinary= binTrigs>S.trigThresh;
else
    
    thresh= repmat(std(decTrigs)*2,size(decTrigs,1),1);
    tTrigsDecimal = decTrigs>thresh;
    
    thresh= repmat(std(binTrigs)*2,size(binTrigs,1),1);
    tTrigsBinary = binTrigs>thresh;
    
end
m1= 'Triggers are not bitwise consistent. ';
m2= 'Attempting a Correction. ';
m3='Please check triggers before epoching.';
msg= [m1,m2,m3];
dBefore= tTrigsDecimal;
    for i = 1:size(tTrigsDecimal,2)
        for j= 1:size(tTrigsDecimal,2)
            di =tTrigsDecimal(:,i)-tTrigsDecimal(:,j);
            prLoc = [find(di==1);find(di==-1)];
            if(all(abs(diff(prLoc))>1))
                tTrigsDecimal(prLoc,i)=1;
                tTrigsDecimal(prLoc,j)=1;
            end
            
        end
    end
 err= sum(sum(abs(dBefore-tTrigsDecimal)));
 if(err>0)
     warning(msg)
 end
 
%-Convert triggers
%--------------------------------------------------------------------------
S.nbits = length(S.decimalTriggerInds);
cFactor = repmat([2.^(0:S.nbits-1)],size(tTrigsDecimal,1),1);
aTrigs  = sum(cFactor.*tTrigsDecimal,2);
uTrigs  = unique(aTrigs);
nTrigs  = length(uTrigs)-1;

outTrigs = zeros(length(aTrigs),nTrigs);

for i =1:nTrigs
    outTrigs(:,i)= (aTrigs==uTrigs(i+1))*uTrigs(i+1);
end

%-Output Struct
%--------------------------------------------------------------------------
lbv=[];
lbv.B= B;
lbv.time= time;
lbv.decimalTrigs= outTrigs;
lbv.binaryTrigs= tTrigsBinary;

function [pos, name] = spm_eeg_elec_pos3d
% returns positions and names of electrodes in biosemi-setup.
%_______________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% James Kilner
% $Id: spm_eeg_elec_pos3d.m 1748 2008-05-28 17:55:25Z guillaume $

inter_ring_dist=2;
pos=zeros(128,4); %x and y coordinates for each electrode

scalars=pi/2:-(pi/2+pi/4)/11:-pi/4;


ring1={'a2';'d15';'d1';'c1';'b1'};
ring2={'a3';'d16';'d14';'d2';'c23';'c2';'b20';'b2'};
ring3={'a4';'junk';'d17';'d18';'d13';'c24';'c22';'c11';'b32';'b21';'b19';'junk'};
ring4={'a19';'a5';'a6';'d28';'d19';'d12';'d3';'c25';'c21';'c12';'c3';'b31';'b22';'b18';'b3';'a32'};
ring5={'a20';'a18';'a7';'d27';'d20';'d11';'d4';'c26';'c20';'c13';'c4';'b30';'b23';'b17';'b4';'a31'};
ring6={'a21';'a17';'a8';'d29';'d26';'d21';'d10';'d5';'c32';'c27';'c19';'c14';'c10';'c5';'b29';'b24';'b16';'b13';'b5';'a30'};
ring7={'a22';'a16';'a9';'d30';'d25';'d22';'d9';'d6';'c31';'c28';'c18';'c15';'c9';'c6';'b28';'b25';'b15';'b12';'b6';'a29'};
ring8={'a23';'a15';'a10';'d31';'d24';'d23';'d8';'d7';'c30';'c29';'c17';'c16';'c8';'c7';'b27';'b26';'b14';'b11';'b7';'a28'};
ring9={'a24';'a14';'a11';'d32';'b10';'b8';'a27'};
ring10={'a25';'a13';'a12';'b9';'a26'};
allelec={'a1';'a2';'d15';'d1';'c1';'b1';'a3';'d16';'d14';'d2';'c23';'c2';'b20';'b2';'a4';'d17';'d18';'d13';'c24';'c22';'c11';'b32';'b21';'b19'; ...
        'a19';'a5';'a6';'d28';'d19';'d12';'d3';'c25';'c21';'c12';'c3';'b31';'b22';'b18';'b3';'a32';'a20';'a18';'a7';'d27';'d20';'d11';'d4';'c26';'c20';'c13';'c4';'b30';'b23';'b17';'b4';'a31'; ...
        'a21';'a17';'a8';'d29';'d26';'d21';'d10';'d5';'c32';'c27';'c19';'c14';'c10';'c5';'b29';'b24';'b16';'b13';'b5';'a30';'a22';'a16';'a9';'d30';'d25';'d22';'d9';'d6';'c31';'c28';'c18';'c15';'c9';'c6';'b28';'b25';'b15';'b12';'b6';'a29'; ...
        'a23';'a15';'a10';'d31';'d24';'d23';'d8';'d7';'c30';'c29';'c17';'c16';'c8';'c7';'b27';'b26';'b14';'b11';'b7';'a28';'a24';'a14';'a11';'d32';'b10';'b8';'a27';'a25';'a13';'a12';'b9';'a26'};
allelecs={'a1';'a2';'a3';'a4';'a5';'a6';'a7';'a8';'a9';'a10';'a11';'a12';'a13';'a14';'a15';'a16';'a17';'a18';'a19';'a20';'a21';'a22';'a23';'a24';'a25';'a26';'a27';'a28';'a29';'a30';'a31';'a32'; ...
        'b1';'b2';'b3';'b4';'b5';'b6';'b7';'b8';'b9';'b10';'b11';'b12';'b13';'b14';'b15';'b16';'b17';'b18';'b19';'b20';'b21';'b22';'b23';'b24';'b25';'b26';'b27';'b28';'b29';'b30';'b31';'b32'; ...
        'c1';'c2';'c3';'c4';'c5';'c6';'c7';'c8';'c9';'c10';'c11';'c12';'c13';'c14';'c15';'c16';'c17';'c18';'c19';'c20';'c21';'c22';'c23';'c24';'c25';'c26';'c27';'c28';'c29';'c30';'c31';'c32'; ...
        'd1';'d2';'d3';'d4';'d5';'d6';'d7';'d8';'d9';'d10';'d11';'d12';'d13';'d14';'d15';'d16';'d17';'d18';'d19';'d20';'d21';'d22';'d23';'d24';'d25';'d26';'d27';'d28';'d29';'d30';'d31';'d32'};

pos(1,1:3)=[0,0,1];
[x,y,z]=sph2cart((0:-(2*pi)/length(ring1):-(2*pi-0.1))-pi/2,ones(1,size(ring1,1)).*scalars(2),1);
pos(2:6,1)=x';
pos(2:6,2)=y';
pos(2:6,3)=z';

[x,y,z]=sph2cart((0:-(2*pi)/length(ring2):-(2*pi-0.1))-pi/2,ones(1,size(ring2,1)).*scalars(3),1);
pos(7:14,1)=x';
pos(7:14,2)=y';
pos(7:14,3)=z';
[x,y,z]=sph2cart((0:-(2*pi)/length(ring3):-(2*pi-0.1))-pi/2,ones(1,size(ring3,1)).*scalars(4),1);
pos(15:24,1)=x([1,3:end-1])';
pos(15:24,2)=y([1,3:end-1])';
pos(15:24,3)=z([1,3:end-1])';
[x,y,z]=sph2cart((0:-(2*pi)/length(ring4):-(2*pi-0.1))-pi/2,ones(1,size(ring4,1)).*scalars(5),1);
pos(25:40,1)=x';
pos(25:40,2)=y';
pos(25:40,3)=z';
[x,y,z]=sph2cart((0:-(2*pi)/length(ring5):-(2*pi-0.1))-pi/2,ones(1,size(ring5,1)).*scalars(6),1);
pos(41:56,1)=x';
pos(41:56,2)=y';
pos(41:56,3)=z';
[x,y,z]=sph2cart((0:-(2*pi)/length(ring6):-(2*pi-0.1))-pi/2,ones(1,size(ring6,1)).*scalars(7),1);
pos(57:76,1)=x';
pos(57:76,2)=y';
pos(57:76,3)=z';
[x,y,z]=sph2cart((0:-(2*pi)/length(ring7):-(2*pi-0.1))-pi/2,ones(1,size(ring7,1)).*scalars(8),1);
pos(77:96,1)=x';
pos(77:96,2)=y';
pos(77:96,3)=z';
[x,y,z]=sph2cart((0:-(2*pi)/length(ring8):-(2*pi-0.1))-pi/2,ones(1,size(ring8,1)).*scalars(9),1);
pos(97:116,1)=x';
pos(97:116,2)=y';
pos(97:116,3)=z';
[x,y,z]=sph2cart((0:-(2*pi)/20:-(2*pi-0.1))-pi/2,ones(1,20).*scalars(10),1);
pos(117:123,1)=x([1,2,3,4,end-2,end-1,end])';
pos(117:123,2)=y([1,2,3,4,end-2,end-1,end])';
pos(117:123,3)=z([1,2,3,4,end-2,end-1,end])';
[x,y,z]=sph2cart((0:-(2*pi)/20:-(2*pi-0.1))-pi/2,ones(1,20).*scalars(11),1);
pos(124:128,1)=x([1,2,3,end-1,end])';
pos(124:128,2)=y([1,2,3,end-1,end])';
pos(124:128,3)=z([1,2,3,end-1,end])';
for i=1:128
    [pos(i,4)]=strmatch(char(allelec(i)),char(allelecs),'exact');
end
    

[i,j]=sort(pos(:,4));
pos(i,4)=pos(j,4);
pos(i,1)=pos(j,1);
pos(i,2)=pos(j,2);
pos(i,3)=pos(j,3);
name = allelecs;

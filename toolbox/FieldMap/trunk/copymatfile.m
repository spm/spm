Fin=spm_get(1,'*.mat','Select single mat file to copy from');
nFin=1;
Fout=spm_get(Inf,'*.mat','Select mat files to copy to');
nFout=size(Fout,1);


for i=1:nFout
     unix(['\cp ' Fin ' ' Fout(i,:)]);
end

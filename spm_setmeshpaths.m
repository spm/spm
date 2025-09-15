function D=spm_setmeshpaths(D,newpath,val)
%function D=spm_setmeshpaths(D,newpath,val)
%% update paths in spm mesh structure
if nargin<3,
    val=[];
end;

if isempty(val),
    val=1;
end;

[a1,b1,c1]=fileparts(D.inv{val}.mesh.tess_iskull);
D.inv{val}.mesh.tess_iskull=[newpath filesep b1 c1];

[a1,b1,c1]=fileparts(D.inv{val}.mesh.tess_oskull);
D.inv{val}.mesh.tess_oskull=[newpath filesep b1 c1];

[a1,b1,c1]=fileparts(D.inv{val}.mesh.tess_scalp);
D.inv{val}.mesh.tess_scalp=[newpath filesep b1 c1];

[a1,b1,c1]=fileparts(D.inv{val}.mesh.tess_ctx);
D.inv{val}.mesh.tess_ctx=[newpath filesep b1 c1];

[a1,b1,c1]=fileparts(D.inv{val}.mesh.sMRI);
D.inv{val}.mesh.sMRI=[newpath filesep b1 c1];

D.save;





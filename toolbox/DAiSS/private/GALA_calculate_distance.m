function distance = GALA_calculate_distance(mesh)

Vertices = mesh.Vertices;
Faces = mesh.Faces;

%first hemisphere
A = spm_mesh_distmtx(struct('vertices',Vertices,'faces',Faces),0);
A = A(1:length(Vertices)/2,1:length(Vertices)/2);
AA = A;
B = A+speye(length(Vertices)/2);
DD1 = full(A);
i=2;
while ~all(all(DD1+eye(length(Vertices)/2)))
    AA = AA*A;
    in = find(AA);
    Bn=zeros(length(Vertices)/2,length(Vertices)/2);
    Bn(in)=1;
    D=full(Bn-B);
    B=Bn;
    DD1 = DD1+i*D;
    i=i+1;
end

%second hemisphere
A = spm_mesh_distmtx(struct('vertices',Vertices,'faces',Faces),0);
A = A(length(Vertices)/2+1:length(Vertices),length(Vertices)/2+1:length(Vertices));
AA = A;
B = A+speye(length(Vertices)/2);
DD2 = full(A);
i=2;
while ~all(all(DD2+eye(length(Vertices)/2)))
    AA = AA*A;
    in = find(AA);
    Bn=zeros(length(Vertices)/2,length(Vertices)/2);
    Bn(in)=1;
    D=full(Bn-B);
    B=Bn;
    DD2 = DD2+i*D;
    i=i+1;
end

distance = Inf*ones(length(Vertices),length(Vertices));
distance(1:length(Vertices)/2,1:length(Vertices)/2)=DD1;
distance(length(Vertices)/2+1:length(Vertices),length(Vertices)/2+1:length(Vertices))=DD2;

return


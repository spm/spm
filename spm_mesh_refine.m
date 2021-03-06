function M = spm_mesh_refine(M)
% Refine a triangle mesh
% FORMAT M = spm_mesh_refine(M)
% M        - a patch structure or gifti object
%__________________________________________________________________________
%
% See also:
%
% R.E. Bank, A.H. Sherman and A. Weiser. Refinement Algorithms and Data 
% Structures for Regular Local Mesh Refinement. Scientific Computing 
% (Applications of Mathematics and Computing to the Physical Sciences)
% (R. S. Stepleman, ed.), North-Holland (1983), 3-17.
% https://ccom.ucsd.edu/~reb/reports/a23.pdf.gz
%__________________________________________________________________________

% Guillaume Flandin
% Copyright (C) 2018-2022 Wellcome Centre for Human Neuroimaging


V  = M.vertices;
F  = M.faces;

Nv = size(V,1);
Nf = size(F,1);

Vo = V;
Fo = zeros(4*Nf,3);
A  = spm_mesh_adjacency(M);
if isfield(M,'cdata')
    cdata = M.cdata;
else
    cdata = [];
end

for f=1:Nf
    T0 = F(f,:);
    T1 = T0([2 3 1]);
    v  = (V(T0,:) + V(T1,:)) / 2;
    if ~isempty(cdata), C = (cdata(T0,:) + cdata(T1,:)) / 2; end
    
    s = 1:3;
    b = [false false false];
    for j=1:3
        if A(T0(j),T1(j)) == 1
            s(j) = Nv + 1;
            A(T0(j),T1(j)) = s(j);
            A(T1(j),T0(j)) = s(j);
            Nv = s(j);
            b(j) = true;
        else
            s(j) = A(T0(j),T1(j));
        end
    end
    Vo(end+1:end+nnz(b),:) = v(b,:);
    if ~isempty(cdata), cdata(end+1:end+nnz(b),:) = C(b,:); end
    T0(4:6) = s;
    
    Fo(4*f+(-3:0),:) = T0([1 4 6;4 2 5;6 5 3;4 5 6]);
end

M.faces = Fo;
M.vertices = Vo;
if ~isempty(cdata), M.cdata = cdata; end

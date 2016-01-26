function M = spm_mesh_polyhedron(name)
% Return one of the Platonic solids with triangle faces
% FORMAT M = spm_mesh_polyhedron(name)
% name     - polyhedron name
%            (one of {'tetrahedron','octahedron','icosahedron'})
% 
% M        - patch structure
%__________________________________________________________________________
%
% See http://en.wikipedia.org/wiki/Platonic_solid
%__________________________________________________________________________
% Copyright (C) 2012 Wellcome Trust Centre for Neuroimaging

% Guillaume Flandin
% $Id: spm_mesh_polyhedron.m 6694 2016-01-26 17:09:11Z guillaume $


switch lower(name)
    
    case 'tetrahedron'
        M.vertices = [  1   0  -1/sqrt(2) ;
                       -1   0  -1/sqrt(2) ;
                        0   1   1/sqrt(2) ;
                        0  -1   1/sqrt(2) ];
        M.faces    = [  1   2   3 ;
                        1   2   4 ;
                        2   3   4 ;
                        3   1   4];

    case 'octahedron'
        M.vertices = [  0   0   1 ;
                        1   0   0 ;
                        0   1   0 ;
                       -1   0   0 ;
                        0  -1   0 ;
                        0   0  -1 ];
        M.faces    = [  1   2   3 ;
                        1   3   4 ;
                        1   4   5 ;
                        1   5   2 ;
                        6   3   2 ;
                        6   4   3 ;
                        6   5   4 ;
                        6   2   5 ];
        
    case 'icosahedron'
        r = (1 + sqrt(5)) / 2;
        M.vertices = [  1   0  -r ;
                        0   r  -1 ;
                        r   1   0 ;
                        r  -1   0 ;
                        0  -r  -1 ;
                       -1   0  -r ;
                       -r   1   0 ;
                        0   r   1 ;
                        1   0   r ;
                        0  -r   1 ;
                       -r  -1   0 ;
                       -1   0   r ];
        M.faces    = [  1   2   3   
                        3   4   1
                        1   4   5
                        5   6   1
                        1   6   2
                        2   6   7
                        7   8   2
                        2   8   3
                        3   8   9
                        9   4   3
                        4   9  10
                        4  10   5
                       10  11   5
                        5  11   6
                        6  11   7
                        7  12   8
                        8  12   9
                        9  12  10
                       10  12  11
                       11  12   7 ];
        
    otherwise
        error('''%s'' is not a Platonic solid with triangle faces.',name);
end

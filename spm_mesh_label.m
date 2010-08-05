function [C, N] = spm_mesh_label(M,opt)
% Label connected components of a surface mesh
% FORMAT C = spm_mesh_label(M)
% M        - a [nx3] faces array or a patch structure
% opt      - return connected components on faces/vertices:
%            {['faces'] ,'vertices'}
%
% C        - a [nx1] vector containing labels for the connected components
%            in M
% N        - number of vertices per connected component
%__________________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% Guillaume Flandin
% $Id: spm_mesh_label.m 4035 2010-08-05 18:54:32Z guillaume $


%-Parse input arguments
%--------------------------------------------------------------------------
if ishandle(M)
    M = struct('faces',get(M,'Faces'), 'vertices',get(M,'Vertices'));
elseif isnumeric(M)
    M = struct('faces',M);
end

if nargin < 2, opt = 'faces'; end

%-Compute the adjacency matrix and perform Dulmage-Mendelsohn decomposition
%--------------------------------------------------------------------------
A = spm_mesh_adjacency(M);
A = A + speye(size(A));
[p,q,r] = dmperm(A);

N = diff(r);

switch lower(opt)
    
    case 'faces'
        %-Label faces
        %------------------------------------------------------------------
        C = zeros(size(M.faces,1),1);
        for i=1:length(r)-1
            C(any(ismember(M.faces,p(r(i):r(i+1)-1)),2)) = i;
        end

    case 'vertices'
        %-Label vertices
        %------------------------------------------------------------------
        C = zeros(size(A,1),1);
        for i=1:length(r)-1
            C(p(r(i):r(i+1)-1)) = i;
        end
        
    otherwise
        error('Unknown option.');
end

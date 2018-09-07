function [B,C] = spm_mesh_borders(M)
% Return borders of a triangle mesh
% FORMAT [B,C] = spm_mesh_borders(M)
% M            - a [nx3] faces array or a patch handle/structure
%
% B            - a [mx1] vector of indices of border vertices
% C            - a cell array of indices of contiguous border vertices
%__________________________________________________________________________
% Copyright (C) 2018 Wellcome Trust Centre for Neuroimaging

% Guillaume Flandin
% $Id: spm_mesh_borders.m 7413 2018-09-07 09:53:34Z guillaume $


%-Parse input arguments
%--------------------------------------------------------------------------
if ishandle(M)
    M = get(M,'Faces');
elseif ~isnumeric(M)
    M = M.faces;
end

%-Compute borders
%--------------------------------------------------------------------------
M  = sort(M,2);
M  = [M(:,[1 2]);M(:,[2 3]);M(:,[1 3])];
[E,IA,IC] = unique(M,'rows'); % M = E(IC,:)
IC = sort(IC);
i  = diff(IC) > 0;
I  = IC([true;i] & [i;true]);
E  = E(I,:);
B  = unique(E);

%-Detect contiguous vertices
%--------------------------------------------------------------------------
if nargout > 1
    C = {};
    while ~isempty(E)
        c = E(1,:);
        E = E(2:end,:);
        while true
            i = any(ismember(E,c),2);
            if any(i)
                c = unique([c reshape(E(i,:),1,[])]);
                E(i,:) = [];
            else
                break;
            end
        end
        C{end+1} = c;
    end
end

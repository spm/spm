function [indField] = spm_getvec(P,fieldName,varargin)
% Get the indices of a specific field in a parameter structure
% FORMAT function [indField] = spm_getvec(P,fieldName,cellInd)
%
% P             - parameter structure
% fieldName     - field name in the parameter structure that is enquired
% cellInd       - index of a specific cell array of parameters (optional)
% indField      - indices of the enquired field 
%
% When using spm_vec, a parameter structure is vectorized field by field,
% cell array by cell array, if any. This function aims at finding the
% indices, in the vectorized parameter, of a particular field.
% When provided with cellInd, the function returns the indices of a
% specific cell array in the different cells that all belong to the given
% field, ie it returns the indices of P.fieldName{cellInd}.
%
%__________________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% Jean Daunizeau
% $Id: spm_getvec.m 2315 2008-10-08 14:43:18Z jean $

if ~isstruct(P)
    error('parameter P must be a struct!')
    return
end

if nargin >2
    cellInd                 = varargin{1};
else
    cellInd                 = [];
end

indField                    = [];
ii                          = 0;
fnames                      = fieldnames(P);
for i=1:length(fnames)
    pf                      = getfield(P,fnames{i});
    vpf                     = spm_vec(pf);
    if ~isequal(fnames{i},fieldName)
        ii                  = ii + size(vpf,1);
    else
        if ~isempty(cellInd)
            if ~iscell(pf)
                str         = ['P.',fnames{i}];
                error(['parameter ',str,' must be a cell array!'])
                return
            elseif numel(pf)< cellInd
                str         = ['P.',fnames{i}];
                error(['parameter ',str,...
                    ' must be a cell array of size at least equal to ',...
                    num2str(cellInd),'!'])
                return
            else
                for j =1:cellInd-1
                    vpf     = spm_vec(pf{j});
                    ii      = ii + size(vpf,1);
                end
                vpf         = spm_vec(pf{cellInd});
                indField    = [ii+1:ii+size(vpf,1)]';
            end
        else
            indField        = [ii+1:ii+size(vpf,1)]';
        end
    end
end


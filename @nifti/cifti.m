function [C,OutputFiles] = cifti(obj)
% Extract CIFTI-2 extension from a NIfTI-2 file and export data
%__________________________________________________________________________
% Copyright (C) 2019 Wellcome Trust Centre for Neuroimaging


C = [];
OutputFiles = {};

%-Check input is a CIFTI file
%--------------------------------------------------------------------------
if ~isfield(obj.hdr,'ext') || isempty(obj.hdr.ext) || obj.hdr.ext.ecode ~= 32
    return;
end

%-Parse CIFTI extension
%--------------------------------------------------------------------------
C.xml = char(obj.hdr.ext.edata(:)');
C.hdr = convert(xmltree(C.xml),[],'attributes',true);
s     = size(obj.dat); s = s(5:end);
C.dat = reshape(obj.dat,s);

%-Check input is a CIFTI-2 file
%--------------------------------------------------------------------------
switch C.hdr.attributes.Version
    case {'2'}
        % noop
    case {'1','1.0'}
        warning('CIFTI-1 file not supported.');
        return;
    otherwise
        warning('Unknown CIFTI version: %s.',C.hdr.attributes.Version);
        return;
end

%-Standard CIFTI Mapping Combinations
%==========================================================================

%-Dense Scalar (.dscalar.nii) or Dense Data Series (.dtseries.nii)
%--------------------------------------------------------------------------
if ismember(obj.hdr.intent_code,[3002, 3006]) % CONNECTIVITY_DENSE_SCALARS
                                              % CONNECTIVITY_DENSE_SERIES
    MapNames = {};
    Surfaces = {};
    Volumes  = {};
    vol      = struct;
    
    %-Read mapping information for each dimension
    %----------------------------------------------------------------------
    for i=1:numel(C.hdr.Matrix.MatrixIndicesMap)
        map = C.hdr.Matrix.MatrixIndicesMap{i};
        switch map.attributes.AppliesToMatrixDimension
            case '0'
                switch map.attributes.IndicesMapToDataType
                    case 'CIFTI_INDEX_TYPE_SCALARS'
                        MapNames = map.NamedMap;
                        if isstruct(MapNames), MapNames = {MapNames}; end
                        MapNames = cellfun(@(x) x.MapName,MapNames,'UniformOutput',false);
                    case 'CIFTI_INDEX_TYPE_SERIES'
                        MapNames = str2double(map.attributes.NumberOfSeriesPoints);
                        MapNames = arrayfun(@(x) sprintf('%04d',x),1:MapNames,'UniformOutput',false);
                    otherwise
                        error('Mismatch between NIfTI intent code and CIFTI index type.');
                end
            case '1'
                if ~strcmp(map.attributes.IndicesMapToDataType,'CIFTI_INDEX_TYPE_BRAIN_MODELS')
                    error('Mismatch between NIfTI intent code and CIFTI index type.');
                end
                if isstruct(map.BrainModel), map.BrainModel = {map.BrainModel}; end
                for j=1:numel(map.BrainModel)
                    switch map.BrainModel{j}.attributes.ModelType
                        case 'CIFTI_MODEL_TYPE_SURFACE'
                            nV  = str2double(map.BrainModel{j}.attributes.SurfaceNumberOfVertices);
                            off = str2double(map.BrainModel{j}.attributes.IndexOffset);
                            iV  = str2num(map.BrainModel{j}.VertexIndices) + 1;
                            if numel(iV) ~= str2double(map.BrainModel{j}.attributes.IndexCount)
                                error('Problem with number of vertices.');
                            end
                            BrainStructure = map.BrainModel{j}.attributes.BrainStructure(17:end);
                            Surfaces{end+1} = struct('nV',nV,'off',off,'iV',iV,'brain',BrainStructure);
                        case 'CIFTI_MODEL_TYPE_VOXELS'
                            off = str2double(map.BrainModel{j}.attributes.IndexOffset);
                            BrainStructure = map.BrainModel{j}.attributes.BrainStructure(17:end);
                            iV  = str2num(map.BrainModel{j}.VoxelIndicesIJK) + 1;
                            iV  = reshape(iV,3,[]);
                            if size(iV,2) ~= str2double(map.BrainModel{j}.attributes.IndexCount)
                                error('Problem with number of voxels.');
                            end
                            Volumes{end+1} = struct('off',off,'iV',iV,'brain',BrainStructure);
                    end
                end
                if ~isempty(Volumes)
                    vol.dim = str2num(map.Volume.attributes.VolumeDimensions);
                    vol.mat = str2num(map.Volume.TransformationMatrixVoxelIndicesIJKtoXYZ);
                    vol.mat = reshape(vol.mat,[4 4])' * [eye(4,3) [-1 -1 -1 1]'];
                    % assumes MeterExponent is -3
                end
            otherwise
                error('Data have to be dense scalars or series.');
        end
    end
    
    %-Export as GIfTI and NIfTI
    %----------------------------------------------------------------------
    for i=1:numel(MapNames)
        for j=1:numel(Surfaces)
            D = NaN(Surfaces{j}.nV,1);
            D(Surfaces{j}.iV) = C.dat(i,Surfaces{j}.off+(1:numel(Surfaces{j}.iV)));
            g = gifti(D);
            OutputFiles{end+1} = sprintf('%s_%s.gii',MapNames{i},Surfaces{j}.brain);
            save(g,OutputFiles{end},'ExternalFileBinary');
        end
        if ~isempty(Volumes)
            OutputFiles{end+1} = sprintf('%s.nii',MapNames{i});
            N = nifti; % using subsasgn as being used within the @nifti class
            N = subsasgn(N,substruct('.','dat'),file_array(OutputFiles{end},vol.dim,'float32-le',352));
            N = subsasgn(N,substruct('.','mat'),vol.mat);
            N = subsasgn(N,substruct('.','mat0'),vol.mat);
            N = subsasgn(N,substruct('.','mat_intent'),'Aligned');
            N = subsasgn(N,substruct('.','mat0_intent'),'Aligned');
            create(N);
            dat = NaN(vol.dim);
            for j=1:numel(Volumes)
                ind = sub2ind(vol.dim,Volumes{j}.iV(1,:),Volumes{j}.iV(2,:),Volumes{j}.iV(3,:));
                dat(ind) = C.dat(i,Volumes{j}.off+(1:numel(ind)));
            end
            N.dat = subsasgn(N.dat,substruct('()',repmat({':'},1,numel(vol.dim))),dat);
        end
    end
end

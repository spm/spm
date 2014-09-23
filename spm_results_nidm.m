function [outdir, prov] = spm_results_nidm(SPM,xSPM,TabDat)
% Export SPM stats results using the NIDASH Data Model (NIDM)
% FORMAT [outdir, prov] = spm_results_nidm(SPM,xSPM,TabDat)
% SPM     - structure containing analysis details (see spm_spm.m)
% xSPM    - structure containing inference details (see spm_getSPM.m)
% TabDat  - structure containing results details (see spm_list.m)
%
% outdir  - output directory
% prov    - provenance object (see spm_provenance.m)
%__________________________________________________________________________
% References:
% 
% Neuroimaging Data Sharing (NIDASH) Data Model (NIDM):
%   http://nidm.nidash.org/
%
% PROV-DM: The PROV Data Model:
%   http://www.w3.org/TR/prov-dm/
%__________________________________________________________________________
% Copyright (C) 2013-2014 Wellcome Trust Centre for Neuroimaging

% Guillaume Flandin
% $Id: spm_results_nidm.m 6193 2014-09-23 20:05:59Z guillaume $


%-Get input parameters, interactively if needed
%==========================================================================
if nargin < 1
    [SPM,xSPM] = spm_getSPM;
elseif nargin < 2
    if isstruct(SPM)
        xSPM = struct('swd',SPM.swd);
    else
        xSPM = struct('swd',spm_file('cpath',SPM));
    end
    [SPM,xSPM] = spm_getSPM(xSPM);
end
if nargin < 3
    % Consider Inf local maxima more than 0mm apart (i.e. all)
    TabDat = spm_list('Table',xSPM,Inf,0);
end

%-Options
%--------------------------------------------------------------------------
gz           = '.gz'; %-Compressed NIfTI {'.gz', ''}
coordsys     = 'nidm:MNICoordinateSystem'; %-Assuming MNI space
NIDMversion  =  '0.2.0';


%==========================================================================
%-Populate output directory
%==========================================================================
if ~exist(SPM.swd,'dir'), SPM.swd = pwd; end
outdir = fullfile(SPM.swd,'nidm');
outdir = spm_file(outdir,'uniquedir');
sts    = mkdir(outdir);
if ~sts, error('Cannot create directory "%s".',outdir); end

%-Design Matrix image (as png and csv)
%--------------------------------------------------------------------------
files.desimg = fullfile(outdir,'DesignMatrix.png');
DesMtx       = (SPM.xX.nKX + 1)*32;
ml           = floor(size(DesMtx,1)/size(DesMtx,2));
DesMtx       = reshape(repmat(DesMtx,ml,1),size(DesMtx,1),[]);
imwrite(DesMtx,gray(64),files.desimg,'png');
files.descsv = fullfile(outdir,'DesignMatrix.csv');
csvwrite(files.descsv,SPM.xX.xKXs.X);

%-Maximum Intensity Projection image (as png)
%--------------------------------------------------------------------------
files.mip    = fullfile(outdir,'MaximumIntensityProjection.png');
MIP          = spm_mip(xSPM.Z,xSPM.XYZmm,xSPM.M,xSPM.units);
imwrite(MIP,gray(64),files.mip,'png');

%-Beta images (as NIfTI)
%--------------------------------------------------------------------------
% for i=1:numel(SPM.Vbeta)
%     files.beta{i} = fullfile(outdir,[sprintf('Beta_%04d',i) '.nii' gz]);
%     img2nii(fullfile(xSPM.swd,SPM.Vbeta(i).fname), files.beta{i});
% end
files.beta = {};

%-SPM{.}, contrast, standard error and ESS images (as NIfTI)
%--------------------------------------------------------------------------
for i=1:numel(xSPM.Ic)
    if numel(xSPM.Ic) == 1, postfix = '';
    else                    postfix = sprintf('_%04d',i); end
    files.spm{i}  = fullfile(outdir,[xSPM.STAT 'Statistic' postfix '.nii' gz]);
    img2nii(fullfile(xSPM.swd,SPM.xCon(xSPM.Ic(i)).Vspm.fname), files.spm{i}, xSPM);
    if xSPM.STAT == 'T'
        files.con{i} = fullfile(outdir,['Contrast' postfix '.nii' gz]);
        img2nii(fullfile(xSPM.swd,SPM.xCon(xSPM.Ic(i)).Vcon.fname), files.con{i},...
            struct('STAT','con'));
        files.conse{i} = fullfile(outdir,['ContrastStandardError' postfix '.nii' gz]);
        Vc = SPM.xCon(xSPM.Ic(i)).c' * SPM.xX.Bcov * SPM.xCon(xSPM.Ic(i)).c;
    img2nii(fullfile(xSPM.swd,SPM.VResMS.fname), files.conse{i}, struct('fcn',@(x) sqrt(x*Vc)));
    elseif xSPM.STAT == 'F'
        files.ess{i} = fullfile(outdir,['ExtraSumOfSquares' postfix '.nii' gz]);
        img2nii(fullfile(xSPM.swd,SPM.xCon(xSPM.Ic(i)).Vcon.fname), files.ess{i});
    end
end

%-Thresholded SPM{.} image (as NIfTI)
%--------------------------------------------------------------------------
files.tspm = fullfile(outdir,['ExcursionSet.nii' gz]);
if ~isempty(gz), files.tspm = spm_file(files.tspm,'ext',''); end
evalc('spm_write_filtered(xSPM.Z,xSPM.XYZ,xSPM.DIM,xSPM.M,'''',files.tspm);');
if ~isempty(gz), gzip(files.tspm); spm_unlink(files.tspm); files.tspm = [files.tspm gz]; end

%-Residual Mean Squares image (as NIfTI)
%--------------------------------------------------------------------------
files.resms = fullfile(outdir,['ResidualMeanSquares.nii' gz]);
img2nii(fullfile(xSPM.swd,SPM.VResMS.fname), files.resms);

%-Resels per Voxel image (as NIfTI)
%--------------------------------------------------------------------------
files.rpv = fullfile(outdir,['ReselsPerVoxel.nii' gz]);
img2nii(fullfile(xSPM.swd,SPM.xVol.VRpv.fname), files.rpv);

%-Analysis mask image (as NIfTI)
%--------------------------------------------------------------------------
files.mask = fullfile(outdir,['Mask.nii' gz]);
img2nii(fullfile(xSPM.swd,SPM.VM.fname), files.mask);

%-Grand mean image (as NIfTI)
%--------------------------------------------------------------------------
files.grandmean = fullfile(outdir,'GrandMean.nii');
sf  = mean(SPM.xX.xKXs.X,1);
Vb  = SPM.Vbeta;
for i=1:numel(Vb), Vb(i).pinfo(1:2,:) = Vb(i).pinfo(1:2,:) * sf(i); end
Vgm = struct(...
    'fname',   files.grandmean,...
    'dim',     Vb(1).dim,...
    'dt',      [spm_type('float32') spm_platform('bigend')],...
    'mat',     Vb(1).mat,...
    'pinfo',   [1 0 0]',...
    'descrip', 'Grand Mean');
Vgm = spm_create_vol(Vgm);
Vgm.pinfo(1,1) = spm_add(Vb,Vgm);
Vgm = spm_create_vol(Vgm);
grandMeanMedian = spm_summarise(files.grandmean,SPM.VM.fname,@median);
if ~isempty(gz), gzip(files.grandmean); spm_unlink(files.grandmean); files.grandmean = [files.grandmean gz]; end

%-Explicit mask image (as NIfTI)
%--------------------------------------------------------------------------
if ~isempty(SPM.xM.VM)
    files.emask = fullfile(outdir,['CustomMask.nii' gz]);
    img2nii(fullfile(xSPM.swd,SPM.xM.VM.fname), files.emask);
else
    files.emask = '';
end

%-Clusters n-ary image (as NIfTI)
%--------------------------------------------------------------------------
files.clust = fullfile(outdir,['ClusterLabels.nii' gz]);
if ~isempty(gz), files.clust = spm_file(files.clust,'ext',''); end
Z   = spm_clusters(xSPM.XYZ);
idx = find(~cellfun(@isempty,{TabDat.dat{:,5}}));
n   = zeros(1,numel(idx));
for i=1:numel(idx)
    [unused,j] = spm_XYZreg('NearestXYZ',TabDat.dat{idx(i),12}',xSPM.XYZmm);
    n(i) = Z(j);
end
n(n) = 1:numel(n);
if max(Z) ~= numel(idx)
    warning('Small Volume Correction not handled yet.');
    n(numel(idx)+1:max(Z)) = 0;
end
Z    = n(Z);
evalc('spm_write_filtered(Z,xSPM.XYZ,xSPM.DIM,xSPM.M,'''',files.clust);');
if ~isempty(gz), gzip(files.clust); spm_unlink(files.clust); files.clust = [files.clust gz]; end

%-Display mask images (as NIfTI)
%--------------------------------------------------------------------------
for i=1:numel(xSPM.Im)
    files.dmask{i} = fullfile(outdir,[sprintf('DisplayMask_%04d.nii',i) gz]);
    if isnumeric(xSPM.Im)
        um = spm_u(xSPM.pm,[SPM.xCon(xSPM.Im(i)).eidf,SPM.xX.erdf],...
            SPM.xCon(xSPM.Im(i)).STAT);
        if ~xSPM.Ex, fcn = @(x) x > um;
        else         fcn = @(x) x <= um; end
        img2nii(SPM.xCon(xSPM.Im(i)).Vspm.fname, files.dmask{i}, struct('fcn',fcn));
    else
        if ~xSPM.Ex, fcn = @(x) x~=0 & ~isnan(x);
        else         fcn = @(x) ~(x~=0 & ~isnan(x)); end
        img2nii(xSPM.Im{i}, files.dmask{i}, struct('fcn',fcn));
    end
end
if numel(xSPM.Im) == 0, files.dmask = {}; end

%-SVC Mask (as NIfTI)
%--------------------------------------------------------------------------
if strcmp(TabDat.tit,'p-values adjusted for search volume')
    files.svcmask = '';
elseif strncmp(TabDat.tit,'search volume: ',15)
    warning('Small Volume Correction not handled yet.'); % see spm_VOI.m
    % '%0.1fmm sphere at [%.0f,%.0f,%.0f]'
    % '%0.1f x %0.1f x %0.1f mm box at [%.0f,%.0f,%.0f]'
    % 'image mask: %s'
    files.svcmask = '';
else
    warning('Unable to retrieve search volume details: assuming whole brain search.');
    files.svcmask = '';
end

%-Search Space mask image (as NIfTI)
%--------------------------------------------------------------------------
files.searchspace = fullfile(outdir,['SearchSpace.nii' gz]);
img2nii(fullfile(xSPM.swd,SPM.VM.fname), files.searchspace);


%==========================================================================
%-                          D A T A   M O D E L
%==========================================================================

clear coordspace originalfile isHumanReadable

niifmt = {'image/nifti','xsd:string'};
isHumanReadable(true);

pp = spm_provenance;

%-Namespaces
%--------------------------------------------------------------------------
pp.add_namespace('neurolex','http://neurolex.org/wiki/');
pp.add_namespace('spm','http://www.incf.org/ns/nidash/spm#');
pp.add_namespace('nidm','http://www.incf.org/ns/nidash/nidm#');
pp.add_namespace('niiri','http://iri.nidash.org/');
pp.add_namespace('crypto','http://id.loc.gov/vocabulary/preservation/cryptographicHashFunctions#');
pp.add_namespace('dct','http://purl.org/dc/terms/');

%-Provenance
%--------------------------------------------------------------------------
idResults = getid('niiri:spm_results_id',isHumanReadable);
pp.entity(idResults,{...
  'prov:type','prov:Bundle',...
  'prov:label','SPM Results',...
  'nidm:objectModel','nidm:SPMResults',...
  'nidm:version',{NIDMversion,'xsd:string'},...
  });
pp.wasGeneratedBy(idResults,'-',now);

p = spm_provenance;

%-Agent: SPM
%--------------------------------------------------------------------------
[V,R] = spm('Ver');
idSoftware = getid('niiri:software_id',isHumanReadable);
p.agent(idSoftware,{...
    'prov:type','nidm:SPM',...
    'prov:type','prov:SoftwareAgent',...
    'prov:label',{'SPM','xsd:string'},...
    'nidm:softwareVersion',{V,'xsd:string'},...
    'spm:softwareRevision',{R,'xsd:string'},...
    });

%-Entity: Coordinate Space
%--------------------------------------------------------------------------
id_data_coordspace = coordspace(p,xSPM.M,xSPM.DIM,xSPM.units,coordsys,1);

%-Entity: Image Data
%--------------------------------------------------------------------------
if isfield(SPM,'Sess')
    extra_fields = {...
        'nidm:grandMeanScaling',{'true','xsd:boolean'},...
        'nidm:targetIntensity',{SPM.xGX.GM,'xsd:float'},...
        };
else
    extra_fields = {...
        'nidm:grandMeanScaling',{'false','xsd:boolean'},...
        };
end
idData = getid('niiri:data_id',isHumanReadable);
p.entity(idData,{...
    'prov:type','prov:Collection',...
    'prov:type','nidm:Data',...
    'prov:label',{'Data','xsd:string'},...
    extra_fields{:}});

%-Entity: Design Matrix
%--------------------------------------------------------------------------
idDesignMatrix = getid('niiri:design_matrix_id',isHumanReadable);
idDesignMatrixImage = getid('niiri:design_matrix_png_id',isHumanReadable);

p.entity(idDesignMatrix,{...
    'prov:type','nidm:DesignMatrix',...
    'prov:location',{uri(spm_file(files.descsv,'cpath')),'xsd:anyURI'},...
    'nidm:filename',{spm_file(files.descsv,'filename'),'xsd:string'},...
    'dct:format',{'text/csv','xsd:string'},...
    'nidm:visualisation',idDesignMatrixImage,...
    'prov:label',{'Design Matrix','xsd:string'},...
    });

p.entity(idDesignMatrixImage,{...
    'prov:type','nidm:Image',...
    'prov:location',{uri(spm_file(files.desimg,'cpath')),'xsd:anyURI'},...
    'nidm:filename',{spm_file(files.desimg,'filename'),'xsd:string'},...
    'dct:format',{'image/png','xsd:string'},...
    });

%-Entity: Explicit Mask
%--------------------------------------------------------------------------
if ~isempty(SPM.xM.VM)
    if ~spm_check_orientations(struct('dim',{xSPM.DIM',SPM.xM.VM.dim},...
            'mat',{xSPM.M,SPM.xM.VM.mat}),false)
        id_emask_coordspace = coordspace(p,SPM.xM.VM.mat,SPM.xM.VM.dim',...
            xSPM.units,coordsys);
    end
    idMask2 = getid('niiri:mask_id_2',isHumanReadable);
    p.entity(idMask2,{...
        'prov:type','nidm:CustomMaskMap',...
        'prov:location',{uri(spm_file(files.emask,'cpath')),'xsd:anyURI'},...
        'nidm:filename',{spm_file(files.emask,'filename'),'xsd:string'},...
        'dct:format',niifmt,...
        'prov:label',{'Custom Mask','xsd:string'},...
        'nidm:atCoordinateSpace',id_emask_coordspace,...
        'crypto:sha512',{sha512sum(spm_file(files.emask,'cpath')),'xsd:string'},...
        });
    id = originalfile(p,files.emask);
    p.wasDerivedFrom(idMask2,id);
end

%-Entity: Noise Model
%--------------------------------------------------------------------------
if isfield(SPM.xVi,'form')
    if strcmp(SPM.xVi.form,'i.i.d')
        extra_fields_NM = {...
            'nidm:hasNoiseDependence','nidm:IndependentNoise',...
            'nidm:noiseVarianceHomogeneous',{'true','xsd:boolean'},...
            };
        extra_fields_PE = {
            'nidm:withEstimationMethod','nidm:OrdinaryLeastSquares',...
            };
    else
        extra_fields_NM = {...
            'nidm:hasNoiseDependence','nidm:SeriallyCorrelatedNoise',...
            'nidm:dependenceSpatialModel','nidm:SpatiallyGlobalModel',...
            'nidm:noiseVarianceHomogeneous',{'true','xsd:boolean'},...
            'nidm:varianceSpatialModel','nidm:SpatiallyLocalModel',...
            };
        extra_fields_PE = {
            'nidm:withEstimationMethod','nidm:GeneralizedLeastSquares',...
            };
    end
else
    if numel(SPM.xVi.Vi) == 1 % assume it's identity
        extra_fields_NM = {...
            'nidm:hasNoiseDependence','nidm:IndependentNoise',...
            'nidm:noiseVarianceHomogeneous',{'true','xsd:boolean'},...
            };
        extra_fields_PE = {
            'nidm:withEstimationMethod','nidm:OrdinaryLeastSquares',...
            };
    else
        extra_fields_NM = {...
            'nidm:hasNoiseDependence','nidm:ArbitrarilyCorrelatedNoise',...
            'nidm:dependenceSpatialModel','nidm:SpatiallyGlobalModel',...
            'nidm:noiseVarianceHomogeneous',{'false','xsd:boolean'},...
            'nidm:varianceSpatialModel','nidm:SpatiallyLocalModel',...
            };
        extra_fields_PE = {
            'nidm:withEstimationMethod','nidm:GeneralizedLeastSquares',...
            };
    end
end
idNoiseModel = getid('niiri:noise_model_id',isHumanReadable);
p.entity(idNoiseModel,{...
    'prov:type','nidm_NoiseModel',...
    'nidm:hasNoiseDistribution','nidm:GaussianDistribution',...
    extra_fields_NM{:}});

%-Activity: Model Parameters Estimation
%==========================================================================
idModelPE = getid('niiri:model_pe_id',isHumanReadable);
p.activity(idModelPE,{...
    'prov:type','nidm:ModelParametersEstimation',...
    'prov:label','Model parameters estimation',...
    extra_fields_PE{:}});
p.wasAssociatedWith(idModelPE, idSoftware);
p.used(idModelPE, idDesignMatrix);
p.used(idModelPE, idData);
p.used(idModelPE, idNoiseModel);
if ~isempty(SPM.xM.VM)
    p.used(idModelPE, idMask2);
end

%-Entity: Mask Map
%--------------------------------------------------------------------------
idMask1 = getid('niiri:mask_id_1',isHumanReadable);
p.entity(idMask1,{...
    'prov:type','nidm:MaskMap',...
    'prov:location',{uri(spm_file(files.mask,'cpath')),'xsd:anyURI'},...
    'nidm:filename',{spm_file(files.mask,'filename'),'xsd:string'},...
    'dct:format',niifmt,...
    'prov:label',{'Mask','xsd:string'},...
    'nidm:atCoordinateSpace',id_data_coordspace,...
    'crypto:sha512',{sha512sum(spm_file(files.mask,'cpath')),'xsd:string'},...
    });
id = originalfile(p,SPM.VM.fname);
p.wasDerivedFrom(idMask1,id);
p.wasGeneratedBy(idMask1, idModelPE);

%-Entity: Grand Mean Map
%--------------------------------------------------------------------------
idGrandMean = getid('niiri:grand_mean_map_id',isHumanReadable);
p.entity(idGrandMean,{...
    'prov:type','nidm:GrandMeanMap',...
    'prov:location',{uri(spm_file(files.grandmean,'cpath')),'xsd:anyURI'},...
    'nidm:filename',{spm_file(files.grandmean,'filename'),'xsd:string'},...
    'dct:format',niifmt,...
    'prov:label',{'Grand Mean Map','xsd:string'},...
    'nidm:maskedMedian',{grandMeanMedian,'xsd:float'},...
    'nidm:atCoordinateSpace',id_data_coordspace,...
    'crypto:sha512',{sha512sum(spm_file(files.grandmean,'cpath')),'xsd:string'},...
    });
p.wasGeneratedBy(idGrandMean, idModelPE);

%-Entity: Beta Maps
%--------------------------------------------------------------------------
idBeta = cell(1,numel(SPM.Vbeta));
for i=1:numel(SPM.Vbeta)
    if ~isempty(files.beta)
        extra_fields = {...
            'prov:location',{uri(files.beta{i}),'xsd:anyURI'},...
            'nidm:filename',{spm_file(files.beta{i},'filename'),'xsd:string'},...
            'dct:format',niifmt,...
            'crypto:sha512',{sha512sum(files.beta{i}),'xsd:string'},...
        };
    else
        extra_fields = {};
    end
    idBeta{i} = getid(sprintf('niiri:beta_map_id_%d',i),isHumanReadable);
    p.entity(idBeta{i},{...
        'prov:type','nidm:ParameterEstimateMap',...
        'prov:label',{sprintf('Beta Map %d',i),'xsd:string'},...
        'nidm:atCoordinateSpace',id_data_coordspace,...
        extra_fields{:},...
    });
    id = originalfile(p,SPM.Vbeta(i).fname);
    p.wasDerivedFrom(idBeta{i},id);
    p.wasGeneratedBy(idBeta{i}, idModelPE);
end

%-Entity: ResMS Map
%--------------------------------------------------------------------------
idResMS = getid('niiri:residual_mean_squares_map_id',isHumanReadable);
p.entity(idResMS,{...
    'prov:type','nidm:ResidualMeanSquaresMap',...
    'prov:location',{uri(spm_file(files.resms,'cpath')),'xsd:anyURI'},...
    'nidm:filename',{spm_file(files.resms,'filename'),'xsd:string'},...
    'dct:format',niifmt,...
    'prov:label',{'Residual Mean Squares Map','xsd:string'},...
    'nidm:atCoordinateSpace',id_data_coordspace,...
    'crypto:sha512',{sha512sum(spm_file(files.resms,'cpath')),'xsd:string'},...
    });
id = originalfile(p,SPM.VResMS.fname);
p.wasDerivedFrom(idResMS,id);
p.wasGeneratedBy(idResMS, idModelPE);

%-Entity: RPV Map
%--------------------------------------------------------------------------
idRPV = getid('niiri:resels_per_voxel_map_id',isHumanReadable);
p.entity(idRPV,{...
    'prov:type','spm:ReselsPerVoxelMap',...
    'prov:location',{uri(spm_file(files.rpv,'cpath')),'xsd:anyURI'},...
    'nidm:filename',{spm_file(files.rpv,'filename'),'xsd:string'},...
    'dct:format',niifmt,...
    'prov:label',{'Resels per Voxel Map','xsd:string'},...
    'nidm:atCoordinateSpace',id_data_coordspace,...
    'crypto:sha512',{sha512sum(spm_file(files.rpv,'cpath')),'xsd:string'},...
    });
id = originalfile(p,SPM.xVol.VRpv.fname);
p.wasDerivedFrom(idRPV,id);
p.wasGeneratedBy(idRPV, idModelPE);

%-Activity: Contrast Estimation
%==========================================================================
for c=1:numel(xSPM.Ic)
    if numel(xSPM.Ic) == 1, postfix = '';
    else                    postfix = sprintf('_%d',c); end
    
    idConVec = getid(['niiri:contrast_id' postfix],isHumanReadable);
    p.entity(idConVec,{...
        'prov:type','nidm:ContrastWeights',...
        'nidm:statisticType',['nidm:' xSPM.STAT 'Statistic'],...
        'nidm:contrastName',{SPM.xCon(xSPM.Ic(c)).name,'xsd:string'},... %esc
        'prov:label',{['Contrast: ' SPM.xCon(xSPM.Ic(c)).name],'xsd:string'},... %esc
        'prov:value',{SPM.xCon(xSPM.Ic(c)).c','xsd:string'},...
        });

    idConEst = getid(['niiri:contrast_estimation_id' postfix],isHumanReadable);
    p.activity(idConEst,{...
        'prov:type','nidm:ContrastEstimation',...
        'prov:label',['Contrast estimation' strrep(postfix,'_',' ')],...
        });
    p.wasAssociatedWith(idConEst, idSoftware);
    p.used(idConEst, idMask1);
    p.used(idConEst, idResMS);
    p.used(idConEst, idDesignMatrix);
    p.used(idConEst,idConVec);
    for i=1:numel(SPM.Vbeta)
        p.used(idConEst, idBeta{i});
    end
    
    idSPM{c} = getid(['niiri:statistic_map_id' postfix],isHumanReadable);
    p.entity(idSPM{c},{...
        'prov:type','nidm:StatisticMap',...
        'prov:location',{uri(spm_file(files.spm{c},'cpath')),'xsd:anyURI'},...
        'nidm:filename',{spm_file(files.spm{c},'filename'),'xsd:string'},...
        'dct:format',niifmt,...
        'prov:label',{['Statistic Map: ' SPM.xCon(xSPM.Ic(c)).name],'xsd:string'},... %esc
        'nidm:statisticType',['nidm:' xSPM.STAT 'Statistic'],...
        'nidm:contrastName',{SPM.xCon(xSPM.Ic(c)).name,'xsd:string'},... %esc
        'nidm:errorDegreesOfFreedom',{xSPM.df(2),'xsd:float'},...
        'nidm:effectDegreesOfFreedom',{xSPM.df(1),'xsd:float'},...
        'nidm:atCoordinateSpace',id_data_coordspace,...
        'crypto:sha512',{sha512sum(files.spm{c}),'xsd:string'},...
        });
    id = originalfile(p,SPM.xCon(xSPM.Ic(c)).Vspm.fname);
    p.wasDerivedFrom(idSPM{c},id);
    p.wasGeneratedBy(idSPM{c},idConEst);
    
    if xSPM.STAT == 'T'
        idContrast = getid(['niiri:contrast_map_id' postfix],isHumanReadable);
        p.entity(idContrast,{...
            'prov:type','nidm:ContrastMap',...
            'prov:location',{uri(spm_file(files.con{c},'cpath')),'xsd:anyURI'},...
            'nidm:filename',{spm_file(files.con{c},'filename'),'xsd:string'},...
            'dct:format',niifmt,...
            'prov:label',{['Contrast Map: ' SPM.xCon(xSPM.Ic(c)).name],'xsd:string'},... %esc
            'nidm:contrastName',{SPM.xCon(xSPM.Ic(c)).name,'xsd:string'},... %esc
            'nidm:atCoordinateSpace',id_data_coordspace,...
            'crypto:sha512',{sha512sum(spm_file(files.con{c},'cpath')),'xsd:string'},...
            });
        id = originalfile(p,SPM.xCon(xSPM.Ic(c)).Vcon.fname);
        p.wasDerivedFrom(idContrast,id);
        p.wasGeneratedBy(idContrast,idConEst);
        
        idSE = getid(['niiri:contrast_standard_error_map_id' postfix],isHumanReadable);
        p.entity(idSE,{...
            'prov:type','nidm:ContrastStandardErrorMap',...
            'prov:location',{uri(spm_file(files.conse{c},'cpath')),'xsd:anyURI'},...
            'nidm:filename',{spm_file(files.conse{c},'filename'),'xsd:string'},...
            'dct:format',niifmt,...
            'prov:label',{'Contrast Standard Error Map','xsd:string'},...
            'nidm:atCoordinateSpace',id_data_coordspace,...
            'crypto:sha512',{sha512sum(spm_file(files.conse{c},'cpath')),'xsd:string'},...
            });
        p.wasGeneratedBy(idSE,idConEst);
    end
    if xSPM.STAT == 'F'
        idESS = getid(['niiri:contrast_extra_sum_of_squares_id' postfix],isHumanReadable);
        p.entity(idESS,{...
            'prov:type','nidm:ContrastExtraSumOfSquaresMap',...
            'prov:location',{uri(spm_file(files.ess{c},'cpath')),'xsd:anyURI'},...
            'nidm:filename',{spm_file(files.ess{c},'filename'),'xsd:string'},...
            'dct:format',niifmt,...
            'prov:label',{'Contrast Extra Sum of Squares Map','xsd:string'},...
            'nidm:atCoordinateSpace',id_data_coordspace,...
            'crypto:sha512',{sha512sum(spm_file(files.ess{c},'cpath')),'xsd:string'},...
            });
        id = originalfile(p,SPM.xCon(xSPM.Ic(c)).Vcon.fname);
        p.wasDerivedFrom(idESS,id);
        p.wasGeneratedBy(idESS,idConEst);
    end
end

%-Entity: Height Threshold
%--------------------------------------------------------------------------
td = regexp(xSPM.thresDesc,'p\D?(?<u>[\.\d]+) \((?<thresDesc>\S+)\)','names');
if isempty(td)
    td = regexp(xSPM.thresDesc,'\w=(?<u>[\.\d]+)','names');
    if ~isempty(td)
        ustt = [xSPM.STAT ' value'];
    else
        ustt = 'unknown';
    end
else
    ustt = ['p-value ' td.thresDesc];
end
idHeightThresh = getid('niiri:height_threshold_id',isHumanReadable);
p.entity(idHeightThresh,{...
    'prov:type','nidm:HeightThreshold',...
    'prov:label',{['Height Threshold: ' xSPM.thresDesc],'xsd:string'},... %esc
    'nidm:userSpecifiedThresholdType',{ustt,'xsd:string'},...
    'prov:value',{TabDat.ftr{1,2}(1),'xsd:float'},... % xSPM.u
    'nidm:pValueUncorrected',{TabDat.ftr{1,2}(2),'xsd:float'},...
    'nidm:pValueFWER',{TabDat.ftr{1,2}(3),'xsd:float'},...
    });

%-Entity: Extent Threshold
%--------------------------------------------------------------------------
if spm_get_defaults('stats.rft.nonstat')
    warning('Non-stationary RFT results not handled yet.');
end
V2R = 1 / prod(xSPM.FWHM);

if TabDat.ftr{2,2}(1) > 0
    extra_fields = {'nidm:pValueUncorrected',{TabDat.ftr{2,2}(2),'xsd:float'},...
    'nidm:pValueFWER',{TabDat.ftr{2,2}(3),'xsd:float'}};
else
    extra_fields = {'nidm:pValueUncorrected',{1,'xsd:float'},...
    'nidm:pValueFWER',{1,'xsd:float'}};
end
idExtentThresh = getid('niiri:extent_threshold_id',isHumanReadable);
p.entity(idExtentThresh,{...
    'prov:type','nidm:ExtentThreshold',...
    'prov:label',{['Extent Threshold: k>=' num2str(TabDat.ftr{2,2}(1))],'xsd:string'},...
    'nidm:clusterSizeInVoxels',{TabDat.ftr{2,2}(1),'xsd:int'},... % xSPM.k
    'spm:clusterSizeInResels',{TabDat.ftr{2,2}(1)*V2R,'xsd:float'},...
    extra_fields{:},...
    });

%-Entity: maxNumberOfPeaksPerCluster & minDistanceBetweenPeaks
%--------------------------------------------------------------------------
% TabDat.str = 'table shows %d local maxima more than %.1fmm apart'
% maxNumberOfPeaksPerCluster = spm_get_defaults('stats.results.volume.nbmax');
% minDistanceBetweenPeaks = spm_get_defaults('stats.results.volume.distmin');
% clusterConnectivityCriterion = 18;

%-Activity: Inference
%==========================================================================
if numel(xSPM.Ic) == 1
    st = {'prov:type','nidm:Inference', ...
          'nidm:hasAlternativeHypothesis','nidm:OneTailedTest',...
          'prov:label','Inference'};
else
    if xSPM.n == 1
        st = {'prov:type','nidm:ConjunctionInference', ...
              'prov:label','Conjunction Inference'};
    else
        st = {'prov:type','nidm:kConjunctionInference', ...
              'prov:label','k-Conjunction Inference', ...
              'nidm:globalNullDegree',xSPM.n};
    end
end
idInference = getid('niiri:inference_id',isHumanReadable);
p.activity(idInference,st);
p.wasAssociatedWith(idInference, idSoftware);
p.used(idInference, idHeightThresh);
p.used(idInference, idExtentThresh);
for c=1:numel(xSPM.Ic)
    if numel(xSPM.Ic) == 1, postfix = '';
    else                    postfix = sprintf('_%d',c); end
    p.used(idInference, idSPM{c});
end
p.used(idInference, idRPV);
p.used(idInference, idMask1);

%-Entity: Display Mask Maps
%--------------------------------------------------------------------------
for i=1:numel(files.dmask)
    V = spm_vol(files.dmask{i});
    if ~spm_check_orientations(struct('dim',{xSPM.DIM',V.dim},...
            'mat',{xSPM.M,V.mat}),false)
        currCoordSpace = coordspace(p,V.mat,V.dim',xSPM.units,coordsys);
    else
        currCoordSpace = id_data_coordspace;
    end
    
    if numel(files.dmask) == 1, postfix = '';
    else                    postfix = sprintf('_%d',i); end
    idDMask = getid(['niiri:filtering_mask_map_id' postfix],isHumanReadable);
    p.entity(idDMask,{...
        'prov:type','nidm:DisplayMaskMap',...
        'prov:location',{uri(spm_file(files.dmask{i},'cpath')),'xsd:anyURI'},...
        'nidm:filename',{spm_file(files.dmask{i},'filename'),'xsd:string'},...
        'dct:format',niifmt,...
        'prov:label',{'Display Mask Map','xsd:string'},...
        'nidm:atCoordinateSpace',currCoordSpace,...
        'crypto:sha512',{sha512sum(spm_file(files.dmask{i},'cpath')),'xsd:string'},...
        });
    id = originalfile(p,files.dmask{i});
    p.wasDerivedFrom(idDMask,id);
    p.used(idInference, idDMask);
end

%-Entity: SVC Mask Map
%--------------------------------------------------------------------------
if ~isempty(files.svcmask)
    V = spm_vol(files.svcmask);
    if ~spm_check_orientations(struct('dim',{xSPM.DIM',V.dim},...
            'mat',{xSPM.M,V.mat}),false)
        currCoordSpace = coordspace(p,V.mat,V.dim',xSPM.units,coordsys);
    else
        currCoordSpace = id_data_coordspace;
    end
    idSVC = getid('niiri:sub_volume_id',isHumanReadable);
    p.entity(idSVC,{...
        'prov:type','nidm:SubVolumeMap',...
        'prov:location',{uri(spm_file(files.svcmask,'cpath')),'xsd:anyURI'},...
        'nidm:filename',{spm_file(files.svcmask,'filename'),'xsd:string'},...
        'prov:label',{'Sub-volume','xsd:string'},...
        'nidm:atCoordinateSpace',currCoordSpace,...
        'crypto:sha512',{sha512sum(spm_file(files.svcmask,'cpath')),'xsd:string'},...
        });
    id = originalfile(p,files.svcmask);
    p.wasDerivedFrom(idSVC,id);
    p.used(idInference, idSVC);
end

%-Entity: Search Space
%--------------------------------------------------------------------------
if spm_get_defaults('stats.rft.nonstat'), rftstat = {'false','xsd:boolean'};
else                                      rftstat = {'true','xsd:boolean'}; end
idSearchSpace = getid('niiri:search_space_id',isHumanReadable);
p.entity(idSearchSpace,{...
    'prov:type','nidm:SearchSpaceMap',...
    'prov:location',{uri(spm_file(files.searchspace,'cpath')),'xsd:anyURI'},...
    'nidm:filename',{spm_file(files.searchspace,'filename'),'xsd:string'},...
    'dct:format',niifmt,...
    'prov:label',{'Search Space Map','xsd:string'}...
    'nidm:atCoordinateSpace',id_data_coordspace,...
    'spm:searchVolumeInVoxels',{xSPM.S,'xsd:int'},...
    'spm:searchVolumeInUnits',{TabDat.ftr{8,2}(1),'xsd:float'},...
    'spm:reselSize',{TabDat.ftr{9,2}(end),'xsd:float'},...
    'spm:searchVolumeInResels',{TabDat.ftr{8,2}(3),'xsd:float'},...
    'spm:searchVolumeReselsGeometry',{xSPM.R,'xsd:string'},...
    'spm:noiseFWHMInVoxels',{xSPM.FWHM,'xsd:string'},...
    'spm:noiseFWHMInUnits',{TabDat.ftr{7,2}(1:3),'xsd:string'},...
    'nidm:randomFieldStationarity',rftstat,...
    'spm:expectedNumberOfVoxelsPerCluster',{TabDat.ftr{3,2},'xsd:float'},...
    'spm:expectedNumberOfClusters',{TabDat.ftr{4,2},'xsd:float'},...
    'spm:heightCriticalThresholdFWE05',{xSPM.uc(1),'xsd:float'},...
    'spm:heightCriticalThresholdFDR05',{xSPM.uc(2),'xsd:float'},...
    'spm:smallestSignifClusterSizeInVoxelsFWE05',{xSPM.uc(3),'xsd:int'},...
    'spm:smallestSignifClusterSizeInVoxelsFDR05',{xSPM.uc(4),'xsd:int'},...
    'crypto:sha512',{sha512sum(spm_file(files.searchspace,'cpath')),'xsd:string'},...
    });
p.wasGeneratedBy(idSearchSpace, idInference);

%-Entity: Excursion Set
%--------------------------------------------------------------------------
if size(TabDat.dat,1) > 0
    c  = TabDat.dat{1,2};
    pc = TabDat.dat{1,1};
else
    c  = 0;
    pc = NaN;
end
idExcursionSet = getid('niiri:excursion_set_id',isHumanReadable);
idClusterLabelsMap = getid('niiri:cluster_label_map_id',isHumanReadable);
idMaximumIntensityProjection = getid('niiri:maximum_intensity_projection_id',isHumanReadable);
p.entity(idExcursionSet,{...
    'prov:type','nidm:ExcursionSet',...
    'prov:location',{uri(spm_file(files.tspm,'cpath')),'xsd:anyURI'},...
    'nidm:filename',{spm_file(files.tspm,'filename'),'xsd:string'},...
    'dct:format',niifmt,...
    'prov:label',{'Excursion Set','xsd:string'},...
    'nidm:numberOfClusters',{c,'xsd:int'},...
    'nidm:pValue',{pc,'xsd:float'},...
    'nidm:hasClusterLabelsMap',idClusterLabelsMap,...
    'spm:hasMaximumIntensityProjection',idMaximumIntensityProjection,...
    'nidm:atCoordinateSpace',id_data_coordspace,...
    'crypto:sha512',{sha512sum(spm_file(files.tspm,'cpath')),'xsd:string'},...
    });
p.entity(idClusterLabelsMap,{...
    'prov:type','nidm:ClusterLabelsMap',...
    'prov:location',{uri(spm_file(files.clust,'cpath')),'xsd:anyURI'},...
    'nidm:filename',{spm_file(files.clust,'filename'),'xsd:string'},...
    'dct:format',niifmt,...
    });
p.entity(idMaximumIntensityProjection,{...
    'prov:type','nidm:Image',...
    'prov:location',{uri(spm_file(files.mip,'cpath')),'xsd:anyURI'},...
    'nidm:filename',{spm_file(files.mip,'filename'),'xsd:string'},...
    'dct:format',{'image/png','xsd:string'}...
    });
p.wasGeneratedBy(idExcursionSet, idInference); 

%-Entity: Clusters
%--------------------------------------------------------------------------
idx = find(~cellfun(@isempty,{TabDat.dat{:,5}}));
idCluster = cell(1,numel(idx));
for i=1:numel(idx)
    iClus = sprintf('%04d',i);
    idCluster{i} = getid(['niiri:cluster_' iClus],isHumanReadable);
    p.entity(idCluster{i},{...
        'prov:type','nidm:Cluster',...
        'prov:label',{['Cluster: ' iClus],'xsd:string'},...
        'nidm:clusterSizeInVoxels',{TabDat.dat{idx(i),5},'xsd:int'},...
        'spm:clusterSizeInResels',{TabDat.dat{idx(i),5}*V2R,'xsd:float'},...
        'nidm:pValueUncorrected',{TabDat.dat{idx(i),6},'xsd:float'},...
        'nidm:pValueFWER',{TabDat.dat{idx(i),3},'xsd:float'},...
        'nidm:qValueFDR',{TabDat.dat{idx(i),4},'xsd:float'},...
        'nidm:clusterLabelId',{num2str(i),'xsd:int'},...
        });
    p.wasDerivedFrom(idCluster{i}, idExcursionSet);
end

%-Entity: Peaks
%--------------------------------------------------------------------------
idx = cumsum(~cellfun(@isempty,{TabDat.dat{:,5}}));
for i=1:size(TabDat.dat,1)
    iPeak  = sprintf('%04d',i);
    idPeak = getid(['niiri:peak_' iPeak],isHumanReadable);
    idCoordinate = getid(['niiri:coordinate_' iPeak],isHumanReadable);
    
    p.entity(idPeak,{...
        'prov:type','nidm:Peak',...
        'prov:label',{['Peak: ' iPeak],'xsd:string'},...
        'prov:location',idCoordinate,...
        'prov:value',{TabDat.dat{i,9},'xsd:float'},...
        'nidm:equivalentZStatistic',{xsdfloat(TabDat.dat{i,10}),'xsd:float'},...
        'nidm:pValueUncorrected',{TabDat.dat{i,11},'xsd:float'},...
        'nidm:pValueFWER',{TabDat.dat{i,7},'xsd:float'},...
        'nidm:qValueFDR',{TabDat.dat{i,8},'xsd:float'},...
        });

    p.entity(idCoordinate,{...
        'prov:type','prov:Location',...
        'prov:type','nidm:Coordinate',...
        'prov:label',{['Coordinate: ' iPeak],'xsd:string'},...
        'nidm:coordinate1',{TabDat.dat{i,12}(1),'xsd:float'},...
        'nidm:coordinate2',{TabDat.dat{i,12}(2),'xsd:float'},...
        'nidm:coordinate3',{TabDat.dat{i,12}(3),'xsd:float'},...
        });
    
    p.wasDerivedFrom(idPeak, idCluster{idx(i)});
end

pp.bundle(idResults,p);

%==========================================================================
%-                  P R O V   S E R I A L I Z A T I O N
%==========================================================================
serialize(pp,fullfile(outdir,'nidm.provn'));
serialize(pp,fullfile(outdir,'nidm.ttl'));
%serialize(pp,fullfile(outdir,'nidm.json'));
%serialize(pp,fullfile(outdir,'nidm.pdf'));
%zip(fullfile(SPM.swd,[spm_file(outdir,'basename'),'.nidm.zip']),outdir)

prov = pp;


%==========================================================================
% function v = xsdfloat(v)
%==========================================================================
function v = xsdfloat(v)
% See http://books.xmlschemata.org/relaxng/ch19-77095.html
if numel(v) == 1 && isinf(v) && v > 0, v = 'INF';  end
if numel(v) == 1 && isinf(v) && v < 0, v = '-INF'; end
if numel(v) == 1 && isnan(v),          v = 'NaN';  end


%==========================================================================
% function str = esc(str)
%==========================================================================
function str = esc(str)
%-Escape
% See http://www.w3.org/TR/html4/charset.html#h-5.3.2
str = strrep(str,'&','&amp;');
str = strrep(str,'<','&lt;');
str = strrep(str,'>','&gt;');
str = strrep(str,'"','&quot;');


%==========================================================================
% function u = uri(u)
%==========================================================================
function u = uri(u)
%-File URI scheme
if ispc, s='/'; else s=''; end
u = ['file://' s strrep(spm_file(u,'cpath'),'\','/')];
e = ' ';
for i=1:length(e)
    u = strrep(u,e(i),['%' dec2hex(e(i))]);
end
u = ['file://./' spm_file(u,'filename')];


%==========================================================================
% function checksum = sha512sum(file)
%==========================================================================
function checksum = sha512sum(file)
%checksum = 'e43b6e01b0463fe7d40782137867a...'; return % TEMP (short)
md   = java.security.MessageDigest.getInstance('SHA-512');
file = spm_file(file,'cpath');
fid  = fopen(file,'rb');
if fid == -1, error('Cannot open "%s".',file); end
md.update(fread(fid,Inf,'*uint8'));
fclose(fid);
checksum = typecast(md.digest,'uint8');
checksum = lower(reshape(dec2hex(checksum)',1,[]));


%==========================================================================
% function checksum = md5sum(data)
%==========================================================================
function checksum = md5sum(data)
%checksum = 'd41d8cd98f00b204e9800998ecf8427e'; return % TEMP (short)
if ~nargin
    data = char(java.util.UUID.randomUUID);
end
md   = java.security.MessageDigest.getInstance('MD5');
if ischar(data)
    md.update(uint8(data));
else
    md.update(typecast(data,'uint8'));
end
checksum = typecast(md.digest,'uint8');
checksum = lower(reshape(dec2hex(checksum)',1,[]));


%==========================================================================
% function img2nii(img,nii,xSPM)
%==========================================================================
function img2nii(img,nii,xSPM)
if nargin == 2, xSPM = struct; end
if ~isfield(xSPM,'STAT'), xSPM.STAT = ''; end
if ~isfield(xSPM,'fcn'), xSPM.fcn = @(x) x; end
if nargin == 1, nii = spm_file(img,'ext','.nii'); end
gz = strcmp(spm_file(nii,'ext'),'gz');
if gz, nii = spm_file(nii,'ext',''); end
ni     = nifti(img);
no     = nifti;
no.dat = file_array(nii,...
                    ni.dat.dim,...
                    ni.dat.dtype,...
                    0,...
                    ni.dat.scl_slope,...
                    ni.dat.scl_inter);
no.mat  = ni.mat;
no.mat_intent = ni.mat_intent;
no.mat0 = ni.mat0;
no.mat0_intent = ni.mat0_intent;
no.descrip = ni.descrip;
switch xSPM.STAT
    case 'T'
        no.intent.name  = ['spm' xSPM.STATstr];
        no.intent.code  = 3;
        no.intent.param = xSPM.df(2);
    case 'F'
        no.intent.name  = ['spm' xSPM.STATstr];
        no.intent.code  = 4;
        no.intent.param = xSPM.df;
    case 'con'
        no.intent.name  = 'SPM contrast';
        no.intent.code  = 1001;
end

create(no);
no.dat(:,:,:) = xSPM.fcn(ni.dat(:,:,:));
if gz
    gzip(nii);
    spm_unlink(nii);
end


%==========================================================================
% function make_ROI(fname,DIM,M,xY)
%==========================================================================
function make_ROI(fname,DIM,M,xY)
gz = strcmp(spm_file(fname,'ext'),'gz');
if gz, fname = spm_file(fname,'ext',''); end
R = struct(...
    'fname',  fname,...
    'dim',    DIM,...
    'dt',     [spm_type('uint8'), spm_platform('bigend')],...
    'mat',    M,...
    'pinfo',  [1,0,0]',...
    'descrip','ROI');
Q    = zeros(DIM);
[xY, XYZmm, j] = spm_ROI(xY, struct('dim',DIM,'mat',M));
Q(j) = 1;
R    = spm_write_vol(R,Q);
if gz
    gzip(R.fname);
    spm_unlink(R.fname);
end


%==========================================================================
% function id = coordspace(p,M,DIM,units,coordsys,idx)
%==========================================================================
function id = coordspace(p,M,DIM,units,coordsys,idx)
persistent index
if nargin == 6
    index = idx;
else
    if isempty(index)
        index = 1;
    else
        index = index + 1;
    end
end
% Convert from first voxel at [1,1,1] to first voxel at [0,0,0]
v2wm = M * [eye(4,3) [1 1 1 1]'];
M    = M(1:3,1:3);
id = getid(['niiri:coordinate_space_id_' num2str(index)],isHumanReadable);
p.entity(id,{...
    'prov:type','nidm:CoordinateSpace',...
    'prov:label',{['Coordinate space ' num2str(index)],'xsd:string'},...
    'nidm:voxelToWorldMapping',{v2wm,'xsd:string'},...
    'nidm:voxelUnits',{units,'xsd:string'},...
    'nidm:voxelSize',{sqrt(diag(M'*M))','xsd:string'},...
    'nidm:inWorldCoordinateSystem',coordsys,...
    'nidm:numberOfDimensions',{numel(DIM),'xsd:int'},...
    'nidm:dimensionsInVoxels',{DIM,'xsd:string'}
    });

%==========================================================================
% function id = originalfile(p,file,idx)
%==========================================================================
function id = originalfile(p,file,idx)
persistent index
if nargin == 3
    index = idx;
else
    if isempty(index)
        index = 1;
    else
        index = index + 1;
    end
end
id = getid(['niiri:map_id_' num2str(index)],isHumanReadable);
p.entity(id,{...
    'prov:type','nidm:Map',...
    'nidm:filename',{spm_file(file,'filename'),'xsd:string'},...
    'dct:format',{'image/nifti','xsd:string'},...
    'crypto:sha512',{sha512sum(spm_file(file,'cpath')),'xsd:string'},...
    });

%==========================================================================
% function id = getid(id,humanReadable,checksum)
%==========================================================================
function id = getid(id,humanReadable,checksum)
if ~humanReadable
    if nargin == 2
        id = md5sum;
    else
        id = md5sum(checksum);
    end
end

%==========================================================================
% function i = isHumanReadable(i)
%==========================================================================
function i = isHumanReadable(i)
persistent isHR
if nargin, isHR = i; end
if isempty(isHR), error('Default not set.'); end
i = isHR;

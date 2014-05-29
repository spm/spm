function [outdir, prov] = spm_results_nidm(SPM,xSPM,TabDat)
% Export SPM stats results using the NIDASH Data Model (NI-DM)
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
% NIDASH Data Model Working Group (NI-DMWG):
%   http://nidm.nidash.org/
%
% PROV-DM: The PROV Data Model:
%   http://www.w3.org/TR/prov-dm/
%__________________________________________________________________________
% Copyright (C) 2013-2014 Wellcome Trust Centre for Neuroimaging

% Guillaume Flandin
% $Id: spm_results_nidm.m 6025 2014-05-29 13:35:51Z guillaume $


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
coordsys     = 'nidm:mniCoordinateSystem'; %-Assuming MNI space
NIDMversion  =  '0.1.0';


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
files.desimg = fullfile(outdir,'DesMtx.png');
DesMtx       = (SPM.xX.nKX + 1)*32;
ml           = floor(size(DesMtx,1)/size(DesMtx,2));
DesMtx       = reshape(repmat(DesMtx,ml,1),size(DesMtx,1),[]);
imwrite(DesMtx,gray(64),files.desimg,'png');
files.descsv = fullfile(outdir,'DesMtx.csv');
csvwrite(files.descsv,SPM.xX.xKXs.X);

%-Maximum Intensity Projection image (as png)
%--------------------------------------------------------------------------
files.mip    = fullfile(outdir,'MIP.png');
MIP          = spm_mip(xSPM.Z,xSPM.XYZmm,xSPM.M,xSPM.units);
imwrite(MIP,gray(64),files.mip,'png');

%-Beta images (as NIfTI)
%--------------------------------------------------------------------------
% for i=1:numel(SPM.Vbeta)
%     files.beta{i} = fullfile(outdir,[spm_file(SPM.Vbeta(i).fname,'basename') '.nii' gz]);
%     img2nii(fullfile(xSPM.swd,SPM.Vbeta(i).fname), files.beta{i});
% end
files.beta = {};

%-SPM{.}, contrast, standard error and ESS images (as NIfTI)
%--------------------------------------------------------------------------
for i=1:numel(xSPM.Ic)
    if numel(xSPM.Ic) == 1, postfix = '';
    else                    postfix = sprintf('_%d',i); end
    files.spm{i}  = fullfile(outdir,['spm' xSPM.STAT postfix '.nii' gz]);
    img2nii(fullfile(xSPM.swd,SPM.xCon(xSPM.Ic(i)).Vspm.fname), files.spm{i}, xSPM);
    if xSPM.STAT == 'T'
        files.con{i} = fullfile(outdir,['contrast' postfix '.nii' gz]);
        img2nii(fullfile(xSPM.swd,SPM.xCon(xSPM.Ic(i)).Vcon.fname), files.con{i},...
            struct('STAT','con'));
        files.conse{i} = fullfile(outdir,['contrastSE' postfix '.nii' gz]);
        Vc = SPM.xCon(xSPM.Ic(i)).c' * SPM.xX.Bcov * SPM.xCon(xSPM.Ic(i)).c;
    img2nii(fullfile(xSPM.swd,SPM.VResMS.fname), files.conse{i}, struct('fcn',@(x) sqrt(x*Vc)));
    elseif xSPM.STAT == 'F'
        files.ess{i} = fullfile(outdir,['ess' postfix '.nii' gz]);
        img2nii(fullfile(xSPM.swd,SPM.xCon(xSPM.Ic(i)).Vcon.fname), files.ess{i});
    end
end

%-Thresholded SPM{.} image (as NIfTI)
%--------------------------------------------------------------------------
files.tspm = fullfile(outdir,['thresh_spm' xSPM.STAT postfix '.nii' gz]);
if ~isempty(gz), files.tspm = spm_file(files.tspm,'ext',''); end
evalc('spm_write_filtered(xSPM.Z,xSPM.XYZ,xSPM.DIM,xSPM.M,'''',files.tspm);');
if ~isempty(gz), gzip(files.tspm); spm_unlink(files.tspm); files.tspm = [files.tspm gz]; end

%-Residual Mean Squares image (as NIfTI)
%--------------------------------------------------------------------------
files.resms = fullfile(outdir,['ResMS.nii' gz]);
img2nii(fullfile(xSPM.swd,SPM.VResMS.fname), files.resms);

%-Resels per Voxel image (as NIfTI)
%--------------------------------------------------------------------------
files.rpv = fullfile(outdir,['RPV.nii' gz]);
img2nii(fullfile(xSPM.swd,SPM.xVol.VRpv.fname), files.rpv);

%-Analysis mask image (as NIfTI)
%--------------------------------------------------------------------------
files.mask = fullfile(outdir,['mask.nii' gz]);
img2nii(fullfile(xSPM.swd,SPM.VM.fname), files.mask);

%-Explicit mask image (as NIfTI)
%--------------------------------------------------------------------------
if ~isempty(SPM.xM.VM)
    files.emask = fullfile(outdir,['custom_mask.nii' gz]);
    img2nii(fullfile(xSPM.swd,SPM.xM.VM.fname), files.emask);
else
    files.emask = '';
end

%-Clusters n-ary image (as NIfTI)
%--------------------------------------------------------------------------
files.clust = fullfile(outdir,['clusters.nii' gz]);
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

%-Filtering mask images (as NIfTI)
%--------------------------------------------------------------------------
for i=1:numel(xSPM.Im)
    files.fmask{i} = fullfile(outdir,[sprintf('fmask_%02d.nii',i) gz]);
    if isnumeric(xSPM.Im)
        um = spm_u(xSPM.pm,[SPM.xCon(xSPM.Im(i)).eidf,SPM.xX.erdf],...
            SPM.xCon(xSPM.Im(i)).STAT);
        if ~xSPM.Ex, fcn = @(x) x > um;
        else         fcn = @(x) x <= um; end
        img2nii(SPM.xCon(xSPM.Im(i)).Vspm.fname, files.fmask{i}, struct('fcn',fcn));
    else
        if ~xSPM.Ex, fcn = @(x) x~=0 & ~isnan(x);
        else         fcn = @(x) ~(x~=0 & ~isnan(x)); end
        img2nii(xSPM.Im{i}, files.fmask{i}, struct('fcn',fcn));
    end
end
if numel(xSPM.Im) == 0, files.fmask = {}; end

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


%==========================================================================
%-                          D A T A   M O D E L
%==========================================================================

pp = spm_provenance;

%-Namespaces
%--------------------------------------------------------------------------
pp.add_namespace('neurolex','http://neurolex.org/wiki/');
pp.add_namespace('spm','http://www.incf.org/ns/nidash/spm#');
pp.add_namespace('nidm','http://www.incf.org/ns/nidash/nidm#');
pp.add_namespace('niiri','http://iri.nidash.org/');
pp.add_namespace('crypto','http://id.loc.gov/vocabulary/preservation/cryptographicHashFunctions#');

%-Provenance
%--------------------------------------------------------------------------
pp.entity('niiri:spm_results_id',{...
  'prov:type','prov:Bundle',...
  'prov:label','SPM Results',...
  'nidm:objectModel','nidm:SPMResults',...
  'nidm:version',NIDMversion,...
  });
pp.wasGeneratedBy('niiri:spm_results_id','-',now);

p = spm_provenance;

%-Agent: SPM
%--------------------------------------------------------------------------
[V,R] = spm('Ver');
p.agent('niiri:software_id',{...
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
        %'nidm:medianIntensity',{SPM.xGX.GM,'xsd:float'},...
        };
else
    extra_fields = {...
        'nidm:grandMeanScaling',{'false','xsd:boolean'},...
        %'nidm:medianIntensity',{SPM.xGX.GM,'xsd:float'},...
        };
end
p.entity('niiri:data_id',{...
    'prov:type','prov:Collection',...
    'prov:type','nidm:Data',...
    'prov:label',{'Data','xsd:string'},...
    extra_fields{:}});

%-Entity: Design Matrix
%--------------------------------------------------------------------------
p.entity('niiri:design_matrix_id',{...
    'prov:type','nidm:DesignMatrix',...
    'prov:location',{uri(spm_file(files.descsv,'cpath')),'xsd:anyURI'},...
    'nidm:visualisation',{uri(spm_file(files.desimg,'cpath')),'xsd:anyURI'},...
    'prov:label',{'Design Matrix','xsd:string'},...
    'nidm:originalFileName',{spm_file(files.descsv,'filename'),'xsd:string'},...
    });

%-Entity: Explicit Mask
%--------------------------------------------------------------------------
if ~isempty(SPM.xM.VM)
    if ~spm_check_orientations(struct('dim',{xSPM.DIM',SPM.xM.VM.dim},...
            'mat',{xSPM.M,SPM.xM.VM.mat}),false)
        id_emask_coordspace = coordspace(p,SPM.xM.VM.mat,SPM.xM.VM.dim',...
            xSPM.units,coordsys);
    end
    
    p.entity('niiri:mask_id_2',{...
        'prov:type','nidm:CustomMaskMap',...
        'prov:location',{uri(spm_file(files.emask,'cpath')),'xsd:anyURI'},...
        'prov:label',{'Custom Mask','xsd:string'},...
        'nidm:originalFileName',{spm_file(files.emask,'filename'),'xsd:string'},...
        'nidm:atCoordinateSpace',id_emask_coordspace,...
        'crypto:sha512',{sha512sum(spm_file(files.emask,'cpath')),'xsd:string'},...
        });
end

%-Activity: Model Parameters Estimation
%==========================================================================
p.activity('niiri:model_pe_id',{...
    'prov:type','nidm:ModelParametersEstimation',...
    'prov:label','Model parameters estimation',...
    });
p.wasAssociatedWith('niiri:model_pe_id', 'niiri:software_id');
p.used('niiri:model_pe_id', 'niiri:design_matrix_id');
p.used('niiri:model_pe_id','niiri:data_id');
if ~isempty(SPM.xM.VM)
    p.used('niiri:model_pe_id','niiri:mask_id_2');
end

%-Entity: Mask Map
%--------------------------------------------------------------------------
p.entity('niiri:mask_id_1',{...
    'prov:type','nidm:MaskMap',...
    'prov:location',{uri(spm_file(files.mask,'cpath')),'xsd:anyURI'},...
    'prov:label',{'Mask','xsd:string'},...
    'nidm:originalFileName',{spm_file(files.mask,'filename'),'xsd:string'},...
    'nidm:atCoordinateSpace',id_data_coordspace,...
    'crypto:sha512',{sha512sum(spm_file(files.mask,'cpath')),'xsd:string'},...
    });
p.wasGeneratedBy('niiri:mask_id_1','niiri:model_pe_id');

%-Entity: Beta Maps
%--------------------------------------------------------------------------
for i=1:numel(SPM.Vbeta)
    if ~isempty(files.beta)
        extra_fields = {...
            'nidm:originalFileName',{spm_file(files.beta{i},'filename'),'xsd:string'},...
            'prov:location',{uri(files.beta{i}),'xsd:anyURI'},...
            'crypto:sha512',{sha512sum(files.beta{i}),'xsd:string'},...
        };
    else
        extra_fields = {};
    end
    p.entity(sprintf('niiri:beta_map_id_%d',i),{...
        'prov:type','nidm:ParameterEstimateMap',...
        'prov:label',{sprintf('Beta Map %d',i),'xsd:string'},...
        'nidm:atCoordinateSpace',id_data_coordspace,...
        extra_fields{:},...
    });
    p.wasGeneratedBy(sprintf('niiri:beta_map_id_%d',i), 'niiri:model_pe_id');
end

%-Entity: ResMS Map
%--------------------------------------------------------------------------
p.entity('niiri:residual_mean_squares_map_id',{...
    'prov:type','nidm:ResidualMeanSquaresMap',...
    'prov:location',{uri(spm_file(files.resms,'cpath')),'xsd:anyURI'},...
    'prov:label',{'Residual Mean Squares Map','xsd:string'},...
    'nidm:originalFileName',{spm_file(files.resms,'filename'),'xsd:string'},...
    'nidm:atCoordinateSpace',id_data_coordspace,...
    'crypto:sha512',{sha512sum(spm_file(files.resms,'cpath')),'xsd:string'},...
    });
p.wasGeneratedBy('niiri:residual_mean_squares_map_id','niiri:model_pe_id');

%-Entity: RPV Map
%--------------------------------------------------------------------------
p.entity('niiri:resels_per_voxel_map_id',{...
    'prov:type','spm:ReselsPerVoxelMap',...
    'prov:location',{uri(spm_file(files.rpv,'cpath')),'xsd:anyURI'},...
    'prov:label',{'Resels per Voxel Map','xsd:string'},...
    'nidm:originalFileName',{spm_file(files.rpv,'filename'),'xsd:string'},...
    'nidm:atCoordinateSpace',id_data_coordspace,...
    'crypto:sha512',{sha512sum(spm_file(files.rpv,'cpath')),'xsd:string'},...
    });
p.wasGeneratedBy('niiri:resels_per_voxel_map_id','niiri:model_pe_id');

%-Activity: Contrast Estimation
%==========================================================================
for c=1:numel(xSPM.Ic)
    if numel(xSPM.Ic) == 1, postfix = '';
    else                    postfix = sprintf('_%d',c); end
    
    idConVec = ['niiri:contrast_id' postfix];
    p.entity(idConVec,{...
        'prov:type',['nidm:' xSPM.STAT 'ContrastWeights'],...
        'nidm:contrastName',{SPM.xCon(xSPM.Ic(c)).name,'xsd:string'},... %esc
        'prov:label',{['Contrast: ' SPM.xCon(xSPM.Ic(c)).name],'xsd:string'},... %esc
        'prov:value',{SPM.xCon(xSPM.Ic(c)).c','xsd:string'},...
        });

    idConEst = ['niiri:contrast_estimation_id' postfix];
    p.activity(idConEst,{...
        'prov:type','nidm:ContrastEstimation',...
        'prov:label',['Contrast estimation' strrep(postfix,'_',' ')],...
        });
    p.wasAssociatedWith(idConEst, 'niiri:software_id');
    p.used(idConEst,'niiri:mask_id_1');
    p.used(idConEst,'niiri:residual_mean_squares_map_id');
    p.used(idConEst,'niiri:design_matrix_id');
    p.used(idConEst,idConVec);
    for i=1:numel(SPM.Vbeta)
        p.used(idConEst,sprintf('niiri:beta_map_id_%d',i));
    end
    
    idSPM = ['niiri:statistical_map_id' postfix];
    p.entity(idSPM,{...
        'prov:type',['nidm:' xSPM.STAT 'StatisticalMap'],...
        'prov:location',{uri(spm_file(files.spm{c},'cpath')),'xsd:anyURI'},...
        'prov:label',{['Statistical Map: ' SPM.xCon(xSPM.Ic(c)).name],'xsd:string'},... %esc
        'nidm:contrastName',{SPM.xCon(xSPM.Ic(c)).name,'xsd:string'},... %esc
        'nidm:originalFileName',{spm_file(files.spm{c},'filename'),'xsd:string'},...
        'nidm:errorDegreesOfFreedom',{xSPM.df(2),'xsd:float'},...
        'nidm:effectDegreesOfFreedom',{xSPM.df(1),'xsd:float'},...
        'nidm:atCoordinateSpace',id_data_coordspace,...
        'crypto:sha512',{sha512sum(files.spm{c}),'xsd:string'},...
        });
    p.wasGeneratedBy(idSPM,idConEst);
    
    if xSPM.STAT == 'T'
        idContrast = ['niiri:contrast_map_id' postfix];
        p.entity(idContrast,{...
            'prov:type','nidm:ContrastMap',...
            'prov:location',{uri(spm_file(files.con{c},'cpath')),'xsd:anyURI'},...
            'prov:label',{['Contrast Map: ' SPM.xCon(xSPM.Ic(c)).name],'xsd:string'},... %esc
            'nidm:contrastName',{SPM.xCon(xSPM.Ic(c)).name,'xsd:string'},... %esc
            'nidm:originalFileName',{spm_file(files.con{c},'filename'),'xsd:string'},...
            'nidm:atCoordinateSpace',id_data_coordspace,...
            'crypto:sha512',{sha512sum(spm_file(files.con{c},'cpath')),'xsd:string'},...
            });
        p.wasGeneratedBy(idContrast,idConEst);
        
        idSE = ['niiri:contrast_standard_error_map_id' postfix];
        p.entity(idSE,{...
            'prov:type','nidm:ContrastStandardErrorMap',...
            'prov:location',{uri(spm_file(files.conse{c},'cpath')),'xsd:anyURI'},...
            'prov:label',{'Contrast Standard Error Map','xsd:string'},...
            'nidm:originalFileName',{spm_file(files.conse{c},'filename'),'xsd:string'},...
            'nidm:atCoordinateSpace',id_data_coordspace,...
            'crypto:sha512',{sha512sum(spm_file(files.conse{c},'cpath')),'xsd:string'},...
            });
        p.wasGeneratedBy(idSE,idConEst);
    end
    if xSPM.STAT == 'F'
        idESS = ['niiri:contrast_extra_sum_of_squares_id' postfix];
        p.entity(idESS,{...
            'prov:type','nidm:ContrastExtraSumOfSquaresMap',...
            'prov:location',{uri(spm_file(files.ess{c},'cpath')),'xsd:anyURI'},...
            'prov:label',{'Contrast Extra Sum of Squares Map','xsd:string'},...
            'nidm:originalFileName',{spm_file(files.ess{c},'filename'),'xsd:string'},...
            'nidm:atCoordinateSpace',id_data_coordspace,...
            'crypto:sha512',{sha512sum(spm_file(files.ess{c},'cpath')),'xsd:string'},...
            });
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
p.entity('niiri:height_threshold_id',{...
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
p.entity('niiri:extent_threshold_id',{...
    'prov:type','nidm:ExtentThreshold',...
    'prov:label',{['Extent Threshold: k>=' num2str(TabDat.ftr{2,2}(1))],'xsd:string'},...
    'nidm:clusterSizeInVoxels',{TabDat.ftr{2,2}(1),'xsd:float'},... % xSPM.k
    'spm:clusterSizeInResels',{TabDat.ftr{2,2}(1)*V2R,'xsd:float'},...
    extra_fields{:},...
    });

%-Activity: Inference
%==========================================================================
if numel(xSPM.Ic) == 1
    st = {'prov:type','nidm:InferenceOneTailed', ...
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
p.activity('niiri:inference_id',st);
p.wasAssociatedWith('niiri:inference_id', 'niiri:software_id');
p.used('niiri:inference_id', 'niiri:height_threshold_id');
p.used('niiri:inference_id', 'niiri:extent_threshold_id');
for c=1:numel(xSPM.Ic)
    if numel(xSPM.Ic) == 1, postfix = '';
    else                    postfix = sprintf('_%d',c); end
    idSPM = ['niiri:statistical_map_id' postfix];
    p.used('niiri:inference_id', idSPM);
end
p.used('niiri:inference_id','niiri:resels_per_voxel_map_id');
p.used('niiri:inference_id','niiri:mask_id_1');

%-Entity: Filtering Mask Maps
%--------------------------------------------------------------------------
for i=1:numel(files.fmask)
    V = spm_vol(files.fmask{i});
    if ~spm_check_orientations(struct('dim',{xSPM.DIM',V.dim},...
            'mat',{xSPM.M,V.mat}),false)
        currCoordSpace = coordspace(p,V.mat,V.dim',xSPM.units,coordsys);
    else
        currCoordSpace = id_data_coordspace;
    end
    
    if numel(files.fmask) == 1, postfix = '';
    else                    postfix = sprintf('_%d',i); end
    idFMask = ['niiri:filtering_mask_map_id' postfix];
    p.entity(idFMask,{...
        'prov:type','nidm:MaskMap',...
        'prov:location',{uri(spm_file(files.fmask{i},'cpath')),'xsd:anyURI'},...
        'prov:label',{'Filtering Mask Map','xsd:string'},...
        'nidm:originalFileName',{spm_file(files.fmask{i},'filename'),'xsd:string'},...
        'nidm:atCoordinateSpace',currCoordSpace,...
        'crypto:sha512',{sha512sum(spm_file(files.fmask{i},'cpath')),'xsd:string'},...
        });
    p.used('niiri:inference_id',idFMask);
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
    p.entity('niiri:sub_volume_id',{...
        'prov:type','nidm:SubVolumeMap',...
        'prov:location',{uri(spm_file(files.svcmask,'cpath')),'xsd:anyURI'},...
        'prov:label',{'Sub-volume','xsd:string'},...
        'nidm:originalFileName',{spm_file(files.svcmask,'filename'),'xsd:string'},...
        'nidm:atCoordinateSpace',currCoordSpace,...
        'crypto:sha512',{sha512sum(spm_file(files.svcmask,'cpath')),'xsd:string'},...
        });
    p.used('niiri:inference_id','niiri:sub_volume_id');
end

%-Entity: Search Space
%--------------------------------------------------------------------------
if spm_get_defaults('stats.rft.nonstat'), rftstat = 'spm:nonStationaryRandomField';
else                                      rftstat = 'spm:stationaryRandomField'; end
p.entity('niiri:search_space_id',{...
    'prov:type','nidm:SearchSpaceMap',...
    'prov:location',{uri(spm_file(files.mask,'cpath')),'xsd:anyURI'},...
    'prov:label',{'Search Space','xsd:string'}...
    'nidm:originalFileName',{spm_file(files.mask,'filename'),'xsd:string'},...
    'nidm:atCoordinateSpace',id_data_coordspace,...
    'spm:searchVolumeInVoxels',{xSPM.S,'xsd:int'},...
    'spm:searchVolumeInProductOfUnits',{TabDat.ftr{8,2}(1),'xsd:float'},...
    'spm:reselSize',{TabDat.ftr{9,2}(end),'xsd:float'},...
    'spm:searchVolumeInResels',{TabDat.ftr{8,2}(3),'xsd:float'},...
    'spm:searchVolumeReselsGeometry',{xSPM.R,'xsd:string'},...
    'spm:noiseFWHMInVoxels',{xSPM.FWHM,'xsd:string'},...
    'spm:noiseFWHMInUnits',{TabDat.ftr{7,2}(1:3),'xsd:string'},...
    'nidm:randomFieldStationarity',rftstat,...
    'crypto:sha512',{sha512sum(spm_file(files.mask,'cpath')),'xsd:string'},...
    });
p.wasGeneratedBy('niiri:search_space_id','niiri:inference_id');

%-Entity: Statistical Image Properties
%--------------------------------------------------------------------------
p.entity('niiri:stat_image_properties_id',{...
    'prov:type','spm:StatisticImageProperties',...
    'prov:label','Statistical image properties',...
    'spm:expectedNumberOfVoxelsPerCluster',{TabDat.ftr{3,2},'xsd:float'},...
    'spm:expectedNumberOfClusters',{TabDat.ftr{4,2},'xsd:float'},...
    'spm:heightCriticalThresholdFWE05',{xSPM.uc(1),'xsd:float'},...
    'spm:heightCriticalThresholdFDR05',{xSPM.uc(2),'xsd:float'},...
    'spm:smallestSignifClusterSizeInVoxelsFWE05',{xSPM.uc(3),'xsd:float'},...
    'spm:smallestSignifClusterSizeInVoxelsFDR05',{xSPM.uc(4),'xsd:float'},...
    });
p.wasGeneratedBy('niiri:stat_image_properties_id', 'niiri:inference_id');

%-Entity: Excursion Set
%--------------------------------------------------------------------------
p.entity('niiri:excursion_set_id',{...
    'prov:type','nidm:ExcursionSet',...
    'prov:location',{uri(spm_file(files.tspm,'cpath')),'xsd:anyURI'},...
    'prov:label',{'Excursion Set','xsd:string'},...
    'nidm:clusterLabelsMap',{uri(spm_file(files.clust,'cpath')),'xsd:anyURI'},...
    'spm:maximumIntensityProjection',{uri(spm_file(files.mip,'cpath')),'xsd:anyURI'},...
    'nidm:originalFileName',{spm_file(files.tspm,'filename'),'xsd:string'},...
    'nidm:atCoordinateSpace',id_data_coordspace,...
    'crypto:sha512',{sha512sum(spm_file(files.tspm,'cpath')),'xsd:string'},...
    });
p.wasGeneratedBy('niiri:excursion_set_id', 'niiri:inference_id'); 

%-Entity: Set
%--------------------------------------------------------------------------
if size(TabDat.dat,1) > 0
    p.entity('niiri:set_statistic_id',{...
        'prov:type','spm:SetStatistic',...
        'prov:label',{'Set Statistic','xsd:string'},...
        'prov:value',{TabDat.dat{1,2},'xsd:int'},...
        'nidm:pValue',{TabDat.dat{1,1},'xsd:float'},...
        });
    p.wasDerivedFrom('niiri:set_statistic_id','niiri:excursion_set_id');
end

%-Entity: Clusters
%--------------------------------------------------------------------------
idx = find(~cellfun(@isempty,{TabDat.dat{:,5}}));
for i=1:numel(idx)
    iClus = sprintf('%04d',i);
    p.entity(['niiri:cluster_' iClus],{...
        'prov:type','nidm:ClusterStatistic',...
        'prov:label',{['Cluster Statistic: ' iClus],'xsd:string'},...
        'nidm:clusterSizeInVoxels',{TabDat.dat{idx(i),5},'xsd:int'},...
        'spm:clusterSizeInResels',{TabDat.dat{idx(i),5}*V2R,'xsd:float'},...
        'nidm:pValueUncorrected',{TabDat.dat{idx(i),6},'xsd:float'},...
        'nidm:pValueFWER',{TabDat.dat{idx(i),3},'xsd:float'},...
        'nidm:qValueFDR',{TabDat.dat{idx(i),4},'xsd:float'},...
        'nidm:clusterLabelId',num2str(i),...
        });
    p.wasDerivedFrom(['niiri:cluster_' iClus],'niiri:set_statistic_id');
end

%-Entity: Peaks
%--------------------------------------------------------------------------
idx = cumsum(~cellfun(@isempty,{TabDat.dat{:,5}}));
for i=1:size(TabDat.dat,1)
    iPeak = sprintf('%04d',i);
    iClus = sprintf('%04d',idx(i));
    p.entity(['niiri:peak_' iPeak],{...
        'prov:type','nidm:PeakStatistic',...
        'prov:label',{['Peak Statistic: ' iPeak],'xsd:string'},...
        'prov:location',['niiri:coordinate_' iPeak],...
        'prov:value',{TabDat.dat{i,9},'xsd:float'},...
        'nidm:equivalentZStatistic',{xsdfloat(TabDat.dat{i,10}),'xsd:float'},...
        'nidm:pValueUncorrected',{TabDat.dat{i,11},'xsd:float'},...
        'nidm:pValueFWER',{TabDat.dat{i,7},'xsd:float'},...
        'nidm:qValueFDR',{TabDat.dat{i,8},'xsd:float'},...
        });

    p.entity(['niiri:coordinate_' iPeak],{...
        'prov:type','prov:Location',...
        'prov:type','nidm:Coordinate',...
        'prov:label',{['Coordinate: ' iPeak],'xsd:string'},...
        'nidm:coordinate1',{TabDat.dat{i,12}(1),'xsd:float'},...
        'nidm:coordinate2',{TabDat.dat{i,12}(2),'xsd:float'},...
        'nidm:coordinate3',{TabDat.dat{i,12}(3),'xsd:float'},...
        });
    
    p.wasDerivedFrom(['niiri:peak_' iPeak], ['niiri:cluster_' iClus]);
end

pp.bundle('niiri:spm_results_id',p);

%==========================================================================
%-                  P R O V   S E R I A L I Z A T I O N
%==========================================================================
serialize(pp,fullfile(outdir,'spm_nidm.provn'));
serialize(pp,fullfile(outdir,'spm_nidm.json'));
serialize(pp,fullfile(outdir,'spm_nidm.ttl'));
%serialize(pp,fullfile(outdir,'spm_nidm.pdf'));
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
u = spm_file(u,'filename'); % TEMP (short)


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
id = ['niiri:coordinate_space_id_' num2str(index)];
p.entity(id,{...
    'prov:type','nidm:CoordinateSpace',...
    'prov:label',{['Coordinate space ' num2str(index)],'xsd:string'},...
    'nidm:voxelToWorldMapping',{v2wm,'xsd:string'},...
    'nidm:voxelUnits',{units,'xsd:string'},...
    'nidm:voxelSize',{sqrt(diag(M'*M))','xsd:string'},...
    'nidm:coordinateSystem',coordsys,...
    'nidm:numberOfDimensions',{numel(DIM),'xsd:int'},...
    'nidm:dimensions',{DIM,'xsd:string'}
    });

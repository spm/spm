function [D] = spm_opm_opreg(S)
% Read magnetometer data and optionally set up forward model
% FORMAT D = spm_opm_create(S)
%   S               - input structure
% Optional fields of S:
%   S.headfile     - path to headshape file      - Default:required
%   S.helmetref      - 3 x 3 matrix of fiducials  - Default:required
%   S.headhelmetref  - 3 x 3 matrix of fiducials  - Default:required
%   S.headfid        - 3 x 3 matrix of fiducials  - Default:required
%   S.headhelmetfid  - 3 x 3 matrix of fiducials  - Default:required
% Output:
%  tHelm       - transformed helmet object
%__________________________________________________________________________
% Copyright (C) 2018-2022 Wellcome Centre for Human Neuroimaging

% Tim Tierney

if ~isfield(S, 'affine'),                S.affine = 1; end

%-read and convert meshes to mm
%--------------------------------------------------------------------------
Native= gifti(S.headfile);

inmetres = sqrt(sum((max(Native.vertices)-min(Native.vertices)).^2))<10;
if inmetres
Native.vertices=Native.vertices*1000;
end


%%-  helmet --> head with helmet
%--------------------------------------------------------------------------
helm2headhelm = spm_eeg_inv_rigidreg(S.headhelmetref',S.helmetref');

% helm_in_headhelm = spm_mesh_transform(helmet,helm2headhelm);
% mesh_plot(helm_in_headhelm,headhelmet);


%-  head with helmet  --> head
%--------------------------------------------------------------------------

headhelm2head = spm_eeg_inv_rigidreg(S.headfid',S.headhelmetfid');
% headhelm_in_head = spm_mesh_transform(headhelmet,headhelm2head);
% mesh_plot(headhelm_in_head,Native);

%- templates 
%--------------------------------------------------------------------------
scalp = gifti(fullfile(spm('dir'), 'canonical','scalp_2562.surf.gii'));
brain = gifti(fullfile(spm('dir'), 'canonical','cortex_5124.surf.gii'));
iskull = gifti(fullfile(spm('dir'), 'canonical','iskull_2562.surf.gii'));
oskull = gifti(fullfile(spm('dir'), 'canonical','oskull_2562.surf.gii'));

%- head space --> template (6 paramater)
%--------------------------------------------------------------------------

if ~isfield(S, 'templatefid')
  fid_template = spm_eeg_fixpnt(ft_read_headshape(fullfile(spm('dir'), 'EEGtemplates','fiducials.sfp')));
  fid_template = fid_template.fid.pnt(1:3,:)';
else
  fid_template = S.templatefid;
end

head2templatescalp = spm_eeg_inv_rigidreg(fid_template,S.headfid');
scalpTemplate6 = spm_mesh_transform(Native,head2templatescalp);
%mesh_plot(scalpTemplate6,scalp);

%- sensors --> template (6 paramater)
%--------------------------------------------------------------------------
sens = S.D.sensors('MEG');
s = sens.coilpos;
senstemplate6 = spm_mesh_transform(s,head2templatescalp*headhelm2head*helm2headhelm);

%- cut everything below the ears off of scanned mesh
%--------------------------------------------------------------------------
p = double(scalpTemplate6.vertices);
cut = min(fid_template(1:3,3));
p(p(:,3)<cut,:) =[];

%- cut everything  left and right of scalp mesh
%--------------------------------------------------------------------------
ml = min(scalp.vertices);
mr = max(scalp.vertices);
cut = (p(:,1)<ml(1)) | (p(:,1)>mr(1)) ;
p(cut,:) =[];

%- cut everything  behind the mesh
%--------------------------------------------------------------------------
ml = min(scalp.vertices);
mr = max(scalp.vertices);
cut = (p(:,2)<ml(2)) | (p(:,2)>mr(2)) ;
p(cut,:) =[];


%- affine registration (12 paramater) 

%--------------------------------------------------------------------------
hmm= spm_eeg_inv_icp(double(scalp.vertices'),p',[],[],[],[],S.affine);
scalpTemplate12 = spm_mesh_transform(scalpTemplate6,hmm);


%-  sensors to template 
%--------------------------------------------------------------------------
sens2temp= hmm*head2templatescalp*headhelm2head*helm2headhelm;
sens=[];
sens.vertices =  s;
senstemplate = spm_mesh_transform(sens,sens2temp);

S.D.inv{1}.mesh = spm_eeg_inv_mesh([], 1);
S.D.inv{1}.mesh.Affine=sens2temp;

%-  registration
%--------------------------------------------------------------------------
temp2sens= inv(sens2temp);
S.D.inv{1}.datareg(1).sensors = S.D.sensors('MEG');
S.D.inv{1}.datareg(1).toMNI =sens2temp;
S.D.inv{1}.datareg(1).fromMNI = temp2sens;
S.D.inv{1}.datareg(1).modality = 'MEG';

%-  save new meshes
%--------------------------------------------------------------------------
tbrain = spm_mesh_transform(brain,temp2sens);
save(tbrain, fullfile(path(S.D),['y_','cortex_5124.surf.gii']));
S.D.inv{1}.mesh.tess_ctx = fullfile(path(S.D),['y_','cortex_5124.surf.gii']);
 
tscalp = spm_mesh_transform(scalp,temp2sens);
save(tscalp, fullfile(path(S.D),['y_','scalp_2562.surf.gii']));
S.D.inv{1}.mesh.tess_scalp  = fullfile(path(S.D),['y_','scalp_2562.surf.gii']);


tiskull = spm_mesh_transform(iskull,temp2sens);
save(tiskull, fullfile(path(S.D),['y_','iskull_2562.surf.gii']));
S.D.inv{1}.mesh.tess_iskull = fullfile(path(S.D),['y_','iskull_2562.surf.gii']);

toskull = spm_mesh_transform(oskull,temp2sens);
save(toskull, fullfile(path(S.D),['y_','oskull_2562.surf.gii']));
S.D.inv{1}.mesh.tess_oskull = fullfile(path(S.D),['y_','oskull_2562.surf.gii']);
save(S.D)

%-  forward model
%--------------------------------------------------------------------------
S.D.inv{1}.forward.voltype = 'Single Shell';
S.D = spm_eeg_inv_forward(S.D);
spm_eeg_inv_checkforward(S.D,1,1);


D=S.D;
save(D);


end
function [] = spm_dcm_voi (DCM_filename,voi_filenames)
% Insert new regions into a DCM model
% FORMAT [] = spm_dcm_voi (DCM_filename,voi_filenames)
%
% DCM_filename      Name of DCM file
% voi_filenames     Cell array of new VOI filenames eg. {'VOI_V1','VOI_V5','VOI_PPC'}
%
% The RT is assumed to be the same as before
%
% This function can be used, for example, to replace subject X's data by subject Y's.
% The model can then be re-estimated without having to go through
% model specification again.
%
% Copyright (C) 2005 Wellcome Department of Imaging Neuroscience

% Will Penny
% $Id$


load(DCM_filename);

% Check we have matching number of regions
n=length(voi_filenames);
if ~(n==DCM.n)
    disp('Error in spm_dcm_voi: mismatching number of regions');
    return
end

for i = 1:n
    load(voi_filenames{i});
    
    DCM.v=size(xY.u,1);
    if (i==1)
        DCM.Y.y=zeros(DCM.v,n);
    end
    DCM.Y.y(:,i)  = xY.u;
    
    DCM.Y.name{i} = xY.name;
    DCM.Y.X0 = xY.X0;
    DCM.Y.Ce = spm_Ce(ones(1,DCM.n)*DCM.v);
    % Store new response variable structure in DCM 
    DCM.xY(i) = xY;
end


%instr=['save ',DCM_filename,' DCM'];
%eval(instr);
if str2num(version('-release'))>=14,
    save(DCM_filename, 'DCM', '-V6');
else
    save(DCM_filename, 'DCM');
end;


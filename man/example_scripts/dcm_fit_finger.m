function [DCM] = dcm_fit_finger (yy,M,U,m)
% Fit DCM model to finger data
% FORMAT [DCM] = dcm_fit_finger (yy,M,U,m)
%
% yy        yy{n} for nth trial data
% M         model structure
% U         input structure
% m         PIF order
%
% DCM       o/p data structure

Nr=length(M.x);
Nt=length(yy);
DCM=[];
for n=1:Nt,
    DCM.xY.y{n}=yy{n};
end
DCM.xY.dt=U.dt;
DCM.xY.pst=U.tims;

DCM.options.Fdcm=[4 8];
% PIF order
switch m,
    case 1,
        DCM.As(:,:,1)=[0 0; 1 0];
        DCM.Bs{1}=zeros(Nr,Nr,1);
    case 2,
        DCM.As(:,:,1)=[0 0; 1 0];
        DCM.As(:,:,2)=[0 0; 1 0];
        DCM.Bs{1}=zeros(Nr,Nr,1);
        DCM.Bs{1}=zeros(Nr,Nr,2);

end

U.X=zeros(Nt,1);
DCM.xU=U;
DCM.xU.oldX=U.X;

% Give DCM estimation the observed initial state
DCM.M=M;

DCM=spm_dcm_phase(DCM);

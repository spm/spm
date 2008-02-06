function [DEM] = spm_dem_initialise(DEM)
% Initialises parameter estimates for DEM structures
% FORMAT [DEM] = spm_dem_initialise(DEM)
%
% DEM.M  - hierarchical model
% DEM.Y  - inputs or data
% DEM.U  - prior expectation of causes
% DEM.X  - observation confounds
%__________________________________________________________________________


% check model class
% -------------------------------------------------------------------------
switch lower(DEM.class)

    % set initial parameter estimates to eigenvectors of data
    % ---------------------------------------------------------------------
    case{'principal component analysis','pca','svd', ...
         'factor analysis','fa','independent component analysis','ica', ...
         'sparse coding','sc'}

        Y     = DEM.Y;
        for i = 1:(length(DEM.M) - 1)
            
            % get eigenspace
            % -------------------------------------------------------------
            u           = spm_svd(Y);
            DEM.M(i).pE = u(:,1:DEM.M(i).m);
            
        end
        
        
    % set initial estimates to ...
    % ---------------------------------------------------------------------
    case{'nonlinear model','nlfa'}

        for i = 1:(length(DEM.M) - 1)
            
            % get dimensions
            % -------------------------------------------------------------
            DEM.M(i).pE = [randn(DEM.M(i).l,DEM.M(i).m)/128 ...
                           randn(DEM.M(i).l,DEM.M(i).m)/128];
            
        end


    otherwise
        
        % Unknown model
        % -----------------------------------------------------------------
        rmfield(DEM.M,'P');
        msgbox('No initialization for this class of model')
end

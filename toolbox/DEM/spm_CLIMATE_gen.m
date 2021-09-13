function [Y,X,dates] = spm_CLIMATE_gen(P,M,U,NPI)
% Generate predictions and hidden states of a CLIMATE model
% FORMAT [Y,X,T] = spm_CLIMATE_gen(P,M,U,NPI)
% P    - model parameters
% M    - model structure (M.T - length of timeseries or data structure)
% U    - number of output variables [default: 2] or indices e.g., [4 5]
% NPI  - nonpharmaceutical intervention
%     NPI(i).period = {'dd-mm-yyyy','dd-mm-yyyy'}; % dates of epidemic
%     NPI(i).param  = {'xyz',...};                 % parameter name
%     NPI(i).Q      = (value1,...);                % parameter name
%     NPI(i).dates  = {'dd-mm-yyyy','dd-mm-yyyy'}; % dates of interevention
%
% Y:  {'Average_temperature'                         }
%     {'Area_total_under_SPEI_1_6'                   }
%     {'normalisedDrought_area_underSPEI_1_6'        }
%     {'Area_total_under_SPEI_2'                     }
%     {'normalisedDrought_area_underSPEI_2'          }
%     {'rainfall_mm'                                 }
%     {'yield_z_scores'                              }
%     {'production_z_scores'                         }
%     {'crop_production_ktons'                       }
%     {'Irrigated_total_area'                        }
%     {'Fertiliser_Consumption_tons'                 }
%     {'Milk__000T_'                                 }
%     {'Egg_millions_'                               }
%     {'Meat__000T_'                                 }
%     {'price_z_scores'                              }
%     {'PRIMARYSECTORCONSTANTPRICES_MillionsInRs_'   }
%     {'TOTALCONSTANTPRICES_MillionsInRs_'           }
%     {'PRIMARYSECTORCONSTANTSHARES_Percent_'        }
%     {'STATEMALEAVERAGEFIELDLABOUR_RsPerDay_'       }
%     {'STATEFEMALEAVERAGEFIELDLABOUR_RsPerDay_'     }
%     {'growth_failure_5_14yo_dalys_per_100k'        }
%     {'child_materNaNl_malnutrition_dalys_per_100k' }
%     {'under5_stunt_dalys_per_100k'                 }
%     {'food_subsidy_million_R'                      }
%     {'ChildWasting5_14YearOldsDALYSPer100000'      }
%     {'ChildUnderweight5_14YearsOld_DALYsPer100_000'}
%
% X       - (M.T x 4) marginal densities over four factors
% location   : {'home','out','ccu','removed','isolated','hospital'};
% infection  : {'susceptible','infected','infectious','Ab +ve','Ab -ve','Vac +ve','Vac Inf'};
% clinical   : {'asymptomatic','symptoms','ARDS','death'};
% diagnostic : {'untested','waiting','PCR +ve','PCR -ve','LFD +ve','LFD -ve'}
%
% Z{t} - joint density over hidden states at the time t
% W    - structure containing time varying parameters
%
% This function returns data Y and their latent states or causes X, given
% the parameters of a generative model. This model is a mean field
% approximation based upon population or density dynamics with certain
% conditional dependencies among the marginal densities over four factors.
% See SPM_covid_priors details. In brief, this routine transforms model
% parameters to (exponentiated) scale parameters and then generates a
% sequence of jointed densities over four factors, after assembling a state
% dependent probability transition matrix. The number in the timeseries is
% specified by M.T.
%
% Equipped with a time-dependent ensemble density, outcome measures are
% then generated as expected values. These include the rate of (new) deaths
% and cases per day. This routine can be extended to generate other
% outcomes, or indeed consider other factorisations of the probability
% transition matrices. The subroutine (spm_COVID_T) creating the
% probability transition matrices given the current states and model
% parameters defines the generative model. This model structure rests upon
% a mean field approximation to the transition probabilities that,
% crucially, depends upon (usually the marginal) densities in question.
% Working through the code below will show how this model is constructed.
%
% A more detailed description of the generative model can be found in the
% body of the script.
%__________________________________________________________________________
% Copyright (C) 2020 Wellcome Centre for Human Neuroimaging

% Karl Friston
% $Id: spm_CLIMATE_gen.m 8151 2021-09-13 09:12:37Z karl $


% The generative model:
%==========================================================================
% In brief,


% setup and defaults (assume new deaths and cases as outcome variables)
%--------------------------------------------------------------------------
if (nargin < 3) || isempty(U), U = 1:3; end         % two outcomes
if (nargin < 4), NPI = [];              end         % interventions

% deal with data structures (asynchronous timeseries)
%--------------------------------------------------------------------------
if isstruct(M.T)
    
    % extract data structure and specify temporal domain
    %----------------------------------------------------------------------
    D   = M.T;
    dT  = max(spm_vec(D.date));
    U   = [D.U];
    
else
    try
        dT = datenum(M.T,'dd-mm-yyyy');
    catch
        dT = M.T;
    end
end

try
    d0 = datenum(M.date,'dd-mm-yyyy');
catch
    d0 = M.date;
end


% dates
%--------------------------------------------------------------------------
dates(1) = d0;
for t = 2:1024
    dates(t) = datenum(datestr(dates(t - 1) + 34,'mm-yyyy'),'mm-yyyy'); %#ok<AGROW>
    if dates(t) > dT,  break,  end
end

% unpack and exponentiate parameters
%--------------------------------------------------------------------------
Q    = spm_vecfun(P,@exp);
R    = Q;                 % baseline parameters

% initial states
%--------------------------------------------------------------------------
x     = P.x(:);
f     = @(x,Q)spm_CLIMATE_f(x,Q);
dt    = gradient(dates)/32;
for t = 1:numel(dates)
    
    
    % time-dependent parameters
    %==================================================================
    
    % interventions (NPI)
    %------------------------------------------------------------------
    for j = 1:numel(NPI)
        
        % start and end dates
        %--------------------------------------------------------------
        dstart = datenum(NPI(j).dates{1},'dd-mm-yyyy') - datenum(NPI(j).period{1},'dd-mm-yyyy');
        dfinal = datenum(NPI(j).dates{2},'dd-mm-yyyy') - datenum(NPI(j).period{1},'dd-mm-yyyy');
        if (dates(t) > dstart) && (dates(t) <= dfinal)
            Q.(NPI(j).param) = NPI(j).Q;
        else
            Q.(NPI(j).param) = R.(NPI(j).param);
        end
    end
    
    % update ensemble density (x)
    %==================================================================
    dfdx = spm_cat(spm_diff(f,x,Q,1)); 
    x    = x + spm_dx(dfdx,f(x,Q),dt(t));
    
    X(t,:) = x(:)';

end


% outcomes
%==================================================================
%     {'Average_temperature'                        }
%     {'normalisedDrought_area_underSPEI_1_6'       }
%     {'normalisedDrought_area_underSPEI_2'         }
%     {'rainfall_mm'                                }
%     {'yield_z_scores'                             }
%     {'crop_production_ktons'                      }
%     {'Irrigated_total_area'                       }
%     {'Fertiliser_Consumption_tons'                }
%     {'Milk__000T_'                                }
%     {'Egg_millions_'                              }
%     {'Meat__000T_'                                }
%     {'price_z_scores'                             }
%     {'PRIMARYSECTORCONSTANTPRICES_MillionsInRs_'  }
%     {'TOTALCONSTANTPRICES_MillionsInRs_'          }
%     {'PRIMARYSECTORCONSTANTSHARES_Percent_'       }
%     {'STATEMALEAVERAGEFIELDLABOUR_RsPerDay_'      }
%     {'growth_failure_5_14yo_dalys_per_100k'       }
%     {'child_materNaNl_malnutrition_dalys_per_100k'}
%     {'food_subsidy_million_R'                     }

% {'Meteorological (1)'      }
% {'Meteorological (2)'      }
% {'Meteorological (3)'      }
% {'Anthropological activity'}
% {'Primary sector activity' }
% {'Precipitation'           }
% {'crop production'         }
% {'Irrigation'              }
% {'Fertilisation'           }
% {'Food production'         }
% {'Food price'              }
% {'Malnutrition'            }

% Average temperature
%------------------------------------------------------------------
Y(:,1) = P.y(1) + X(:,1:3)*P.y(2:4);

% Extreme drought
%------------------------------------------------------------------
Y(:,2) = P.y(7)*erf((P.y(5) - X(:,6)) * P.y(6));

% Exceptional drought
%------------------------------------------------------------------
Y(:,3) = P.y(7)*erf((P.y(8) - X(:,6)) * P.y(6));

% rainfall
%------------------------------------------------------------------
Y(:,4) = P.y(1) + X(:,1:3)*P.y(2:4);

    
% retain specified output variables
%==========================================================================
Y = Y(:,U);

% vectorise if data are asynchronous
%==========================================================================
if exist('D','var')
    for t = 1:numel(D)
        j      = ismember(dates,D(t).date);
        D(t).Y = Y(j,t);
    end
    Y  = spm_vec(D.Y);
end

return

function dxdt = spm_CLIMATE_f(x,P)

% P.x(1)  = 8;               % 'Meteorological (1)'                         ;
% P.x(2)  = 8;               % 'Meteorological (2)'                         ;
% P.x(3)  = 24;              % 'Meteorological (3)'                         ;
% P.x(4)  = 1;               % 'Anthropological activity'                   ;
% P.x(5)  = 1;               % 'Primary sector activity'                    ;
% P.x(6)  = 0;               % 'Precipitation'                              ;
% P.x(7)  = 0;               % 'Crop production'                            ;
% P.x(8)  = 0;               % 'Irrigation'                                 ;
% P.x(9)  = 0;               % 'Fertiliser'                                 ;
% P.x(10) = 0;               % 'Food production'                            ;
% P.x(11) = 0;               % 'Food price'                                 ;
% P.x(12) = 0;               % 'Malnutrition'                               ;

% dynamics and parameters of a Lorentz system (with Jacobian)
%==========================================================================
% dxdt = f(x):  see notes at the end of this script
%--------------------------------------------------------------------------
dxdt     = spm_zeros(x);

% Meteorological
%--------------------------------------------------------------------------
rho      = P.L(3) + P.LR*x(4);
dxdt(1)  = ( P.L(1)*x(2) - P.L(1)*x(1))/P.L(4);
dxdt(2)  = ( rho*x(1) - x(2) - x(1)*x(3))/P.L(4);
dxdt(3)  = (-P.L(2)*x(3) + x(1)*x(2))/P.L(4);

% Anthropological activity
%--------------------------------------------------------------------------
dxdt(4)  = x(4)*P.R(2) - x(4)*(P.R(1) + P.R(3)*x(12));
dxdt(5)  = x(5)*P.R(2) - x(5)*(P.R(1) + P.R(4)*x(12));

% Precipitation
%--------------------------------------------------------------------------
dxdt(6)  = erf((x(3) - P.P(1))*P.P(2)) - x(6)*P.P(3);

% Crop production
%--------------------------------------------------------------------------
dxdt(7)  = x(6)*P.C(1) + x(8)*P.C(2) + x(9)*P.C(3) - x(7)*P.C(4);

% Irrigation
%--------------------------------------------------------------------------
dxdt(8)  = erf((P.I(1) - x(6))*P.I(2)) - x(8)*P.I(3);

% Fertilisation
%--------------------------------------------------------------------------
demand   = x(4) - (x(7) + 1);
dxdt(9)  = demand*P.E(1) + x(9)*P.E(2);

% Food production
%--------------------------------------------------------------------------
demand   = x(4) - (x(10) + 1);
climate  = erf((x(3) - P.FL(1))*P.FL(2));
dxdt(10) = demand*P.F(1) + climate*P.F(2) - x(10)*P.F(3);

% Food price
%--------------------------------------------------------------------------
dxdt(11) = x(7)*P.Z(1) + x(9)*P.Z(2) + x(10)*P.Z(3) - x(11)*P.Z(4);

% Malnutrition
%--------------------------------------------------------------------------
dxdt(12) = x(11)*P.M(1) - x(11)*P.M(2);



return

% Jacobian
%------------------------------------------------------------------
dfdx = [];








function [on,rt,ac] = spm_ADEM_cue_rt(DEM)
% returns reaction times and accuracy for ADEM_cued_response demo
% FORMAT [on,rt,ac] = spm_ADEM_cue_rt(DEM)
%
% DEM - DEM structure from ADEM_cued_response.m
%
% on  - cue onset
% ac  - accuracy
% rt  - reaction time
%__________________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging
 
% Karl Friston
% $Id: spm_ADEM_cue_rt.m 4231 2011-03-07 21:00:02Z karl $
 
% distance from target and cue contrast
%--------------------------------------------------------------------------
n   = length(DEM.M(1).x.a);                   % number of targets
F   = DEM.pU.v{1}((1:2) + 2,:);               % location of finger
L   = DEM.pP.P{1};                            % location of targets
for i = 1:n
    D(i,:) = (F(1,:) - L(1,i)).^2;
    D(i,:) = (F(2,:) - L(2,i)).^2 + D(i,:);   % distance from targets
end        
C   = DEM.pU.v{1}((1:n) + 4,:);               % contrast of targets
 
r   = 1/32;                                   % radius of proximity
c   = diff(C > 1,1,2) > 0;                    % target onset
on  = {};                                     % cue onset
ac  = {};                                     % accuracy
rt  = {};                                     % reaction time
 
% get performance
%--------------------------------------------------------------------------
for i = 1:n
    
    on{i} = find(c(i,:));
    for j = 1:length(on{i})
        try
            
            % minimum distance
            %--------------------------------------------------------------
            d        = D(i,(1:8) + on{i}(j))';
            ac{i}(j) = sqrt(min(d));
            
            % estimated reaction time
            %--------------------------------------------------------------
            X        = (1:length(d))' - 1;
            B        = pinv([X.^0 X.^1])*log(d);
            rt{i}(j) = (log(r) - B(1))/B(2);
            
        catch
            ac{i}(j) = NaN;
            rt{i}(j) = NaN;
        end
    end
end
 
% sort trials (over all targets
%--------------------------------------------------------------------------
on = spm_vec(on); [i j] = sort(on,1,'ascend'); on = on(j);
ac = spm_vec(ac); ac = ac(j);
rt = spm_vec(rt); rt = rt(j);
 
% remove first trial
%--------------------------------------------------------------------------
on(1) = [];
rt(1) = [];
ac(1) = [];
 
% convert spatial error to accuracy
%--------------------------------------------------------------------------
ac    = 1./ac;
 
% convert reaction time to ms
%--------------------------------------------------------------------------
dt    = 64/1000;
rt    = rt*1000*dt;
on    = on*dt;

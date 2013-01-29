function spm_dcm_display(varargin)
% Region and anatomical graph display
% FORMAT spm_dcm_display(xY,a,c,h)
% xY    - cell of region structures (see spm_regions)
% a     - connections of directed graph a(i,j,1) = p value; 
%                                       a(i,j,2) = MAP estimate value
% c     - node-specific inputs
% h     - figure handle [default: Graphics window]
%__________________________________________________________________________
% Copyright (C) 2002-2011 Wellcome Trust Centre for Neuroimaging
 
% Karl Friston
% $Id: spm_dcm_display.m 5219 2013-01-29 17:07:07Z spm $
 
 
% input arguments
%--------------------------------------------------------------------------
n = length(varargin);
 
% get dimensions
%--------------------------------------------------------------------------
if n < 1, xY = [];  else xY = varargin{1}; end
if n < 2, a  = [];  else a  = varargin{2}; end
if n < 3, c  = [];  else c  = varargin{3}; end
if n < 4
    Fgraph  = spm_figure('GetWin','Graphics');
    ha = gca;
else
    ha = varargin{4};
    Fgraph = get(ha,'parent');
end

% graphics parameters
%--------------------------------------------------------------------------
col   = {'r','g','b','c','y','m'};
rad   = 16;                   % radius of self-connections
w     = 4;                    % line width
M     = 32;                   % MarkerSize for regions
U     = 0.9;                  % Display threshold on p-value

% get dimensions, locations and names
%--------------------------------------------------------------------------
m     = size(xY,2);
L     = [];
for i = 1:m
    L       = [L xY(i).xyz];
    name{i} = xY(i).name(1:min(end,3));
end
L     = [L; ones(1,m)];
o     = mean(L,2);
M1    = spm_matrix(-o(1:3)');
 
 
% get orientation with the greatest dispersion
%--------------------------------------------------------------------------
for i = 1:3
    u{i}      = eye(4,3);
    u{i}(:,i) = [];
    s(i)      = det(u{i}'*M1*L*L'*M'*u{i});
end
[i,j]   = max(s);
u       = u{j};
 
% compute projection matrix for 'principal' plane
%--------------------------------------------------------------------------
M2      = u';
M2(4,4) = 1;
M1(j,4) = 0;
L       = M2*M1*L;
 
% coordinates
%--------------------------------------------------------------------------
i       = (min(L(1,:)) - M):(max(L(1,:)) + M);
j       = (min(L(2,:)) - M):(max(L(2,:)) + M);
x       = kron(ones(size(j)),i);
y       = kron(j,ones(size(i)));
xyz     = [x; y; zeros(1,length(x)); ones(1,length(x))];
xyz     = pinv(M1)*pinv(M2)*xyz;
M3      = [1 0 0 -min(i); 0 -1 0 max(j); 0 0 0 0; 0 0 0 1];
L       = M3*L;
 
 
% get T1 background
%--------------------------------------------------------------------------
V       = spm_vol(fullfile(spm('Dir'),'canonical','single_subj_T1.nii'));
ijk     = V.mat \ xyz;
t1      = spm_sample_vol(V,ijk(1,:),ijk(2,:),ijk(3,:),2);
t1      = (64 - 16) + 16*t1/max(t1(:));
 
 
% Watermark and regions
%--------------------------------------------------------------------------
str     = get(get(ha,'Title'),'String');
image(rot90(reshape(t1,length(i),length(j))),'parent',ha)
axis(ha,'image','off')
title(ha,str)
 
 
% Connections
%--------------------------------------------------------------------------
Q     = [-pi:pi/32:pi];
Q     = rad*[sin(Q); cos(Q)];
q     = 1/3;
for i = 1:length(a)
    for j = 1:length(a)
        if ~isnan(a(i,j,1))
 
            % show connections
            %--------------------------------------------------------------
            if i ~= j
 
                % line
                %----------------------------------------------------------
                k = rem(j - 1,length(col)) + 1;
                h = line(L(1,[i j]),L(2,[i j]),'Color',col{k},...
                        'LineStyle',':',...
                        'LineWidth',w);
 
                % if significant
                %----------------------------------------------------------
                if a(i,j,1) > U
                    set(h,'LineStyle','-','LineWidth',w)
 
                    % text
                    %------------------------------------------------------
                    u     = q*(L(1,j) - L(1,i)) + L(1,i);
                    v     = q*(L(2,j) - L(2,i)) + L(2,i);
                    str   = {};
                    for k = 1:size(a,3)
                        str{k} = sprintf('%0.2f ',a(i,j,k));
                    end
                    h     = text(u,v,1,str(:),'FontSize',12,...
                                'HorizontalAlignment','Center');
                end
 
            % self-connection
            %--------------------------------------------------------------
            else
 
                % line
                %----------------------------------------------------------
                k     = rem(i - 1,length(col)) + 1;
                u     = L(1,i);
                v     = L(2,i);
                u     = Q(1,:) + u;
                v     = Q(2,:) + v;
                h     = line(u,v,'Color',col{k},'LineStyle',':','LineWidth',w);
 
                % if significant
                %----------------------------------------------------------
                if a(i,j,1) > U
                    set(h,'LineStyle','-','LineWidth',w)
 
                    % text
                    %------------------------------------------------------
                    u     = u(48);
                    v     = v(48);
                    str   = {};
                    for k = 1:size(a,3)
                        str{k} = sprintf('%0.2f ',a(i,j,k));
                    end
                    h     = text(u,v,1,str(:),'FontSize',12,...
                                'HorizontalAlignment','Center');
                end
            end
        end
    end
end
 
% Exogenous inputs
%--------------------------------------------------------------------------
for i = 1:size(c,1)
    if ~isnan(c(i,1))
 
        % line
        %------------------------------------------------------------------
        k     = rem(i - 1,length(col)) + 1;
        u     = L(1,i);
        v     = L(2,i);
        u     = [u (rad + u)];
        v     = [v v];
        h     = line(u,v,'Color',col{k},'LineStyle',':','LineWidth',w);
 
        % if significant
        %------------------------------------------------------------------
        if c(i,1) > U
            set(h,'LineStyle','-','LineWidth',w)
 
            % patch
            %--------------------------------------------------------------
            u     = u(2);
            v     = v(2);
            str   = {};
            for k = 1:size(c,2)
                str{k} = sprintf('%0.2f ',c(i,k));
            end
            h     = text(u,v,str(:),'FontSize',12,...
                                    'HorizontalAlignment','Center');
        end
    end
end
 
 
% projected coordinates of voxels within region[s]
%--------------------------------------------------------------------------
hold on
for i = 1:m
    k = rem(i - 1,length(col)) + 1;
    line(L(1,i),L(2,i),...
        'Color',col{k},...
        'Marker','.',...
        'LineStyle','none',...
        'MarkerSize',98);
 
    text(L(1,i),L(2,i),name{i},'FontSize',16,...
        'FontWeight','Bold',...
        'Color','w',...
        'HorizontalAlignment','center')
end
hold off

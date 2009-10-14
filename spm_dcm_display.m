function spm_dcm_display(varargin)
% Region and anatomical graph display
% FORMAT spm_dcm_display(xY,a,c,u,M,U)
% xY    - cell of region structures (see spm_regions)
% a     - connections of directed graph a(i,j,1) = p value
% c     - node-specific inputs
% u     - (3 x 2) projection matrix     [default = [] ]
% M     - margin (mm)               [default = 24 ]
% U     - theshold for plotting connections     [default = 0.9]
%__________________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% Karl Friston
% $Id: spm_dcm_display.m 3465 2009-10-14 15:14:29Z guillaume $


% null input arguments
%--------------------------------------------------------------------------
Fgraph  = spm_figure('GetWin','Graphics');
n       = length(varargin);

if n < 1; xY = [];  else,   xY = varargin{1};    end
if n < 2; a  = [];  else,   a  = varargin{2};    end
if n < 3; c  = [];  else,   c  = varargin{3};    end
if n < 4; u  = [];  else,   u  = varargin{4};    end
if n < 5; M  = 16;  else,   M  = varargin{5};    end
if n < 6; U  = .9;  else,   U  = varargin{6};    end

col     = [0 0 0];
rad     = 6;
w       = 2;

% enforce orientation for 1 or 2 regions
%--------------------------------------------------------------------------
if length(xY) < 3
    u = [[1 0 0];[0 1 0]]';
end

% get dimensions
%--------------------------------------------------------------------------
m       = size(xY,2);
L       = [];
for i = 1:m
    L       = [L xY(i).xyz];
    name{i} = xY(i).name;
end
L       = [L; ones(1,m)];
o       = mean(L,2);
M1      = spm_matrix(-o(1:3)');

% compute projection matrix for 'principal' plane
%--------------------------------------------------------------------------
if isempty(u)
    [u s v] = svd(M1*L);
    u       = u(:,[2 1]);
    [i j]   = max(abs(u));
    u       = u*diag(sign([u(j(1),1) u(j(2),2)]));
end
M2      = u';
M2(4,4) = 1;
L       = M2*M1*L;

% coordinates
%--------------------------------------------------------------------------
i       = (min(L(1,:)) - M):(max(L(1,:)) + M);
j       = (min(L(2,:)) - M):(max(L(2,:)) + M);
x       = kron(ones(size(j)),i);
y       = kron(j,ones(size(i)));
xyz     = [x; y; zeros(1,length(x)); ones(1,length(x))];
xyz     = pinv(M1)*pinv(M2)*xyz;
M3      = [1 0 0 -min(i); 0 -1 0 -min(j); 0 0 0 0; 0 0 0 1];
L       = M3*L;


% get T1 background
%--------------------------------------------------------------------------
V       = spm_vol(fullfile(spm('Dir'),'canonical','single_subj_T1.nii'));
ijk     = inv(V.mat)*xyz;
t1      = spm_sample_vol(V,ijk(1,:),ijk(2,:),ijk(3,:),2);
t1      = (64 - 16) + 16*t1/max(t1(:));


% Watermark and regions
%--------------------------------------------------------------------------
str     = get(get(gca,'Title'),'String');
image(rot90(reshape(t1,length(i),length(j))))
axis image off
title(str)


% Connections
%--------------------------------------------------------------------------
Q     = [-pi:pi/32:pi];
Q     = rad*[sin(Q); cos(Q)];
x     = mean(L(1,:));
y     = mean(L(2,:));
q     = 1/3;
for i = 1:length(a)
for j = 1:length(a)
    if ~isnan(a(i,j,1))

        % show connection
        %------------------------------------------------------------------
        if i ~= j

            % line
            %--------------------------------------------------------------
            h = line(L(1,[i j]),L(2,[i j]),...
                'Color',col,...
                'LineStyle',':',...
                'LineWidth',w);

            % if significant
            %--------------------------------------------------------------
            if a(i,j,1) > U
                set(h,'LineStyle','-','LineWidth',w)

                % text
                %----------------------------------------------------------
                u     = q*(L(1,j) - L(1,i)) + L(1,i);
                v     = q*(L(2,j) - L(2,i)) + L(2,i);
                str   = {};
                for k = 1:size(a,3)
                    str{k}   = sprintf('%0.2f ',a(i,j,k));
                end
                h     = text(u,v,1,str(:),'FontSize',10,...
                        'HorizontalAlignment','Center');
                if a(i,j,2) < 0
                    set(h,  'Color','r')
                end
            end

        % self connection
        %------------------------------------------------------------------
        else

            % line
            %--------------------------------------------------------------
            u     = (L(1,i) - x);
            v     = (L(2,i) - y);
            l     = sqrt(u^2 + v^2);
            l     = (l + rad)/l;
            u     = Q(1,:) + x + l*u;
            v     = Q(2,:) + y + l*v;
            h     = line(u,v,...
                'Color',col,...
                'LineStyle',':',...
                'LineWidth',w);

            % if significant
            %--------------------------------------------------------------
            if a(i,j,1) > U
                set(h,'LineStyle','-','LineWidth',w)

                % text
                %----------------------------------------------------------
                u     = u(48);
                v     = v(48);
                str   = {};
                for k = 1:size(a,3)
                    str{k}   = sprintf('%0.2f ',a(i,j,k));
                end     
                h     = text(u,v,1,str(:),'FontSize',10,...
                    'HorizontalAlignment','Center');
                if a(i,j,2) < 0
                    set(h,  'Color','r')
                end
            end
        end
    end
end
end

% Extrinsic influences
%--------------------------------------------------------------------------
for i = 1:size(c,1)
    if ~isnan(c(i,1))

        % line
        %------------------------------------------------------------------
        u     = (L(1,i) - x);
        v     = (L(2,i) - y);
        l     = sqrt(u^2 + v^2);
        l     = (l + rad)/l;
        u     = [u l*u] + x;
        v     = [v l*v] + y;
        h     = line(u,v,...
            'Color',col,...
            'LineStyle',':',...
            'LineWidth',w);
        
        % if significant
        %------------------------------------------------------------------
        if c(i,1) > U
            set(h,'LineStyle','-','LineWidth',w)

            % patch
            %--------------------------------------------------------------
            u     = u(2);
            v     = v(2);

            % text
            %--------------------------------------------------------------
            str   = {};
            for k = 1:size(c,2)
                str{k}   = sprintf('%0.2f ',c(i,k));
            end
            h     = text(u,v,str(:),'FontSize',10,...
                'HorizontalAlignment','Center');
            if c(i,2) < 0
                set(h,  'Color','r')
            end

        end
    end
end


% projected coordinates of voxels within region[s]
%--------------------------------------------------------------------------
hold on
for i = 1:m
    l      = xY(i).XYZmm;
    n      = size(l,2);
    l      = [l; ones(1,n)];
    l      = M3*M2*M1*l;
    plot(l(1,:),l(2,:),'.r','MarkerSize',4)
end

line(L(1,:),L(2,:),...
        'Color',[0 0 0],...
        'Marker','.',...
        'LineStyle','none',...
        'MarkerSize',64);

text(L(1,:),L(2,:),name,'FontSize',8,...
            'FontWeight','Bold',...
            'Color','w',...
            'HorizontalAlignment','center',...
            'FontAngle','italic')
hold off

%% Can we recreate the breathing and bifurcation modes?

% addpath(genpath('/Users/samst/Dropbox/Hakai/'));
% addpath(genpath('/Users/samst/Dropbox/UBC/GSW/'));
% addpath(genpath('/Users/samst/Dropbox/UBC/m_map/'));
% addpath(genpath('/Users/samst/Dropbox/UBC/Misc/'));

%% Load data
clear

% Initialize structure
B = struct();
C= struct();
G= struct();

info = ncinfo('G1.nc');
vars = {info.Variables.Name};

s={'B','C','G'};

% Loop through all variables
for f=1:length(s)
    for i=1:length(vars)
        varname = vars{i};

        % Read data from both files
        data1 = ncread([s{f},'1.nc'],varname);
        data2= ncread([s{f},'2.nc'],varname);

        % Get dimension names
        dims = info.Variables(i).Dimensions;
        dim_names = string({dims.Name});

        % Concatenate along 4th dimension if 'time' is present as 4th dim
        if any(dim_names == "time")
            if length(dims)==1
                eval([s{f},'.(varname) = double(cat(1, data1, data2));'])
            elseif length(dims)==3
                eval([s{f},'.(varname) = double(squeeze(cat(3, data1(1,1,:), data2(1,1,:))));'])
            elseif length(dims)==4
                eval([s{f},'.(varname) = double(squeeze(cat(4, data1(1,1,:,:), data2(1,1,:,:))));'])
            end
        else
                eval([s{f},'.(varname) = double(data1);'])
        end
    end
    eval([s{f}, '.time = datenum(datetime(1970,1,1) + seconds(', s{f}, '.time));']);
    eval([s{f}, '.pr = gsw_p_from_z(',s{f}, '.depth.*-1, ',s{f}, '.latitude(1));']);
    eval(['SA = gsw_SA_from_SP(', s{f}, '.so, ', s{f}, '.pr, ', s{f}, '.longitude, ', s{f}, '.latitude(1));']);
    eval(['CT = gsw_CT_from_pt(SA, ', s{f}, '.thetao);']);
    eval([s{f}, '.s_t = gsw_sigma0(SA, CT);']);
    eval([s{f}, '.dyn_height = gsw_geo_strf_dyn_height(SA, CT, ', s{f}, '.pr, 1000);']);
end

D_NPC=(C.dyn_height(1,:)-G.dyn_height(1,:))./9.7963;
D_GAk=(B.dyn_height(1,:)-G.dyn_height(1,:))./9.7963;
D_CC=(C.dyn_height(1,:)-B.dyn_height(1,:))./9.7963;
%% Plot deltas
col=[255 214 140]/255; % YELLOW!
lat_lim=[50 53]; lon_lim=[-132 -127];

figure('units','centimeters','outerposition',[0 0 13 17],'color','w');
axGlob=axes('position',[0.1 0.6 0.75 0.35], 'clipping' , 'on');
m_proj('equidistant', 'lon', [190 260], 'lat', [15 60]); % North Pacific bounds
hold on;

gO=load('globalOxI.mat');
[gO.lat_grid,gO.lon_grid]=meshgrid(gO.lat,gO.lon);
lon_grid_corrected= mod(gO.lon_grid + 360, 360);
cmap = cmocean('balance');
cmap = [[0.5 0.5 0.5]; cmap];  % Add grey for NaNs at the start
colormap(gca, cmap);

% Set NaN values to a specific value outside the contour range
data = gO.dyn_heightT(:,:,1)./9.81; 
msk=gO.lon_grid<-100 | gO.lon_grid>120;
data(~msk)=NaN;
data(isnan(data)) = -0.1;  % Replace NaNs with a value outside the plotted range

% Plot the data with adjusted NaN handling
[CS, CH] = m_contourf(lon_grid_corrected, gO.lat_grid(:,:,1), data, ...
    [0:0.05:2.5 -0.1], 'linestyle', 'none');

% Set the color axis to include the NaN color
clim([0 2.5]);  % Ensure the colormap includes the NaN value (-1)
m_coast('patch', col, 'edgecolor', 'none'); % Light gray continents

lbls={'BC','CA','OR','WA'};

count=0;
Mp = m_shaperead('ne_50m_admin_1_states_provinces'); 
for k = [38 54 87 97]
    lonC = mod(Mp.ncst{k}(:,1) + 360, 360);  % Convert longitudes to 0-360
    latC = Mp.ncst{k}(:,2);

    % Find NaN indices to break the segments
    nanIdx = isnan(lonC) | isnan(latC);
    segmentStart = [1; find(nanIdx) + 1];  % Start of each segment
    segmentEnd = [find(nanIdx) - 1; length(lonC)];  % End of each segment

    % Remove invalid indices (where start > end)
    validSeg = segmentStart <= segmentEnd;
    segmentStart = segmentStart(validSeg);
    segmentEnd = segmentEnd(validSeg);

    % Loop through segments and plot them separately
    for j = 1:length(segmentStart)
        m_plot(lonC(segmentStart(j):segmentEnd(j)), ...
             latC(segmentStart(j):segmentEnd(j)), 'color',[0.6 0 0], 'LineWidth', 0.1);
    end
    count=count+1;
    % m_text(mean(lonC,'omitnan'),mean(latC,'omitnan'), lbls{count});
end

Mp = m_shaperead('ne_50m_admin_0_countries'); 
for k = [17 56 76 90 203]
    lonC = mod(Mp.ncst{k}(:,1) + 360, 360);  % Convert longitudes to 0-360
    latC = Mp.ncst{k}(:,2);

    % Find NaN indices to break the segments
    nanIdx = isnan(lonC) | isnan(latC);
    segmentStart = [1; find(nanIdx) + 1];  % Start of each segment
    segmentEnd = [find(nanIdx) - 1; length(lonC)];  % End of each segment

    % Remove invalid indices (where start > end)
    validSeg = segmentStart <= segmentEnd;
    segmentStart = segmentStart(validSeg);
    segmentEnd = segmentEnd(validSeg);

    % Loop through segments and plot them separately
    for j = 1:length(segmentStart)
        m_plot(lonC(segmentStart(j):segmentEnd(j)), ...
             latC(segmentStart(j):segmentEnd(j)), 'color',[0.6 0.6 0.6], 'LineWidth', 0.2);
    end
end

% Highlight Queen Charlotte Sound region
m_patch(mod(lon_lim([1 2 2 1 1]) + 360, 360), lat_lim([1 1 2 2 1]), 'r', ...
    'facealpha', 0, 'edgecolor', 'r','linewidth',2); % Semi-transparent red box
% 
% m_text(247.97,53.24,'Canada','fontweight','bold')
% m_text(250.51,42.28,'U.S.A.','fontweight','bold')

% m_scatter(mod(lon_lim(1) + 360, 360),lat_lim(1),100,'ro','linewidth',2);

% Add annotation for Queen Charlotte Sound
% m_text(mean(lon_lim), mean(lat_lim), 'QCS', 'color', 'k', ...
%     'fontsize', 6, 'fontweight', 'bold', 'horizontalalignment', 'center');

% Add grid lines for the orthographic map
m_grid('xtick',4 , 'ytick', 3, 'linestyle', 'none'); % Hide grid lines
text(0.05, 0.95,'a)','fontsize',9,'units','normalized');
%
ax=m_contfbar(0.9,[.1 .9],CS,CH);
ylabel(ax,'Dynamic metres','Interpreter','latex');

m_scatter(mod(G.longitude(1) + 360, 360),G.latitude(1),100,'filled','kp');
m_text(mod(G.longitude(1) + 360, 360),G.latitude(1),'\bf G');

m_scatter(mod(C.longitude(1) + 360, 360),C.latitude(1),100,'filled','kp');
m_text(mod(C.longitude(1) + 360, 360),C.latitude(1),'\bf C');

m_scatter(mod(B.longitude(1) + 360, 360),B.latitude(1),100,'filled','kp');
m_text(mod(B.longitude(1) + 360, 360),B.latitude(1),'\bf B');

ax2=axes('position',[0.1 0.3 0.75 0.25]);
plot(G.time,D_NPC,'k','linewidth',0.7);
hold on
% scatter(G.time,D_NPC,1,'filled','markerfacecolor','k');
plot(G.time,D_GAk./D_NPC,'color',rgb_x('pinkish red'),'linewidth',0.7);
% scatter(G.time,D_GAk./D_NPC,1,'filled','markerfacecolor','r');
plot(G.time,D_GAk,'color',rgb_x('blue'),'linewidth',0.7);
% scatter(G.time,D_GAk,1,'filled','markerfacecolor','b');
plot(G.time,D_CC,'color',rgb_x('moss'),'linewidth',0.7);
% scatter(G.time,D_CC,1,'filled','markerfacecolor','g');
grid on; axis tight; axdate(20); xticklabels([]);
ylabel('Dynamic metres');
text(0.05, 0.95,'b)','fontsize',9,'units','normalized');
% legend('\DeltaD_{NPC}','\DeltaG_{Ak}/\DeltaD_{NPC}','\DeltaG_{Ak}','\DeltaG_{CCurr}');

text(1.05, 0.6,'\DeltaD_{NPC}','fontsize',8,'units', ...
    'normalized','color',rgb_x('black'));
text(1.05, 0.55,{'\DeltaG_{Ak}/';'\DeltaD_{NPC}'},'fontsize',8,'units', ...
    'normalized','color',rgb_x('pinkish red'));
text(1.05, 0.4,'\DeltaG_{Ak}','fontsize',8,'units', ...
    'normalized','color',rgb_x('blue'));
text(1.05, 0.35,'\DeltaG_{CCurr}','fontsize',8,'units', ...
    'normalized','color',rgb_x('moss'));

load fRatio.mat
load fD_NPC.mat
load fD_CC.mat
load fD_GAk.mat

% ax3=axes('position',[0.25 0.1 0.5 0.15]);
% plot(G.time,smoothdata(D_GAk./D_NPC,6,'movmean'),'k','linewidth',1);
% hold on
% plot(G.time,smoothdata(D_NPC,6,'movmean'),'r','linewidth',1);
% plot(G.time,smoothdata(D_GAk,6,'movmean'),'b','linewidth',1);
% plot(G.time,smoothdata(D_CC,6,'movmean'),'g','linewidth',1);
% 
% plot(fRatio(:,1),smoothdata(fRatio(:,2),6,'movmean'),'k--');
% plot(fD_NPC(:,1),smoothdata(fD_NPC(:,2),6,'movmean'),'r--');
% plot(fD_GAk(:,1),smoothdata(fD_GAk(:,2),6,'movmean'),'b--');
% plot(fD_CC(:,1),smoothdata(fD_CC(:,2),6,'movmean'),'g--');
% grid on; xlim([datenum(2002,01,01),datenum(2006,01,01)]);
% axdate(4);

% Freeland method: calculate COV
% Demean
D_NPC_anom = D_NPC - mean(D_NPC,'omitnan');
D_GAk_anom = D_GAk - mean(D_GAk,'omitnan');
D_CC_anom = D_CC - mean(D_CC,'omitnan');

%Stack
X = [D_NPC_anom(:), D_GAk_anom(:), D_CC_anom(:)];

% Covariance
CC = cov(X,'omitrows');

% Eigen-decomposition
[V, D] = eig(CC);
[eigenvalues, idx] = sort(diag(D),'descend');
V = V(:,idx);
D = diag(eigenvalues);

% Variance explained
variance_explained = 100 * eigenvalues / sum(eigenvalues);

% Project onto modes
mode_time_series = X * V;


% % Cummins and Freeland method
% Q = [1,  0.5,  1;
%      1,  0.5,  0;
%      1, -0.5, -1];
% [Q_orth, ~] = qr(Q, 0);
% mode_time_series_fixed = X * Q_orth;
% var_total = sum(var(X,'omitnan'));
% var_modes = var(mode_time_series_fixed,'omitnan');
% variance_explained_fixed = 100 * var_modes / var_total;
% 
% % Define fixed-mode basis (3D, unit vectors)
% v_breathing = [1; 1; 1] / norm([1; 1; 1]);
% v_bifurcation = [0; 1; -1] / norm([0; 1; -1]);
% 
% % Dot products: correlation between each PCA mode and the physical modes
% dot_breathing = mode_time_series * v_breathing;
% dot_bifurcation = mode_time_series * v_bifurcation;
% r_breathing = corr(dot_breathing(:,1), mode_time_series_fixed(:,1), 'rows', 'pairwise');
% r_bifurcation = corr(dot_bifurcation(:,1), mode_time_series_fixed(:,2), 'rows', 'pairwise');
% 
% % PCA eigenvectors (columns)
% v1 = V(:,1);   % leading mode
% v2 = V(:,2);   % secondary mode
% 
% % Fixed physical vectors
% b = [1; 1; 1] / norm([1; 1; 1]);
% bf = [0; 1; -1] / norm([0; 1; -1]);
% 
% % Projections
% dot1_breathing = dot(v1, b);
% dot1_bifurcation = dot(v1, bf);
% dot2_breathing = dot(v2, b);
% dot2_bifurcation = dot(v2, bf);

% In this dataset, standard PCA yields eigenvectors nearly identical 
% to the physically defined breathing and bifurcation modes described by 
% Cummins & Freeland (2007), confirming the internal consistency of the 
% ΔD series and the robustness of the circulation patterns.

%In this dataset, the PCA identifies the same two patterns of variability
% described by Cummins & Freeland (2007): a breathing mode that reflects 
% changes in total transport, and a bifurcation mode that reflects how 
% transport is split between the Alaska and California Currents. This
% shows that the data are internally consistent, and that the PCA c
% aptures physically meaningful patterns in the flow.

ax3=axes('position',[0.1 0.135 0.75 0.15]);
plot(G.time, mode_time_series(:,1),'color',rgb_x('violet'));
hold on;
plot(G.time, mode_time_series(:,2),'color',rgb_x('teal'));
axdate(20)
% yyaxis right
% load NEPGpcsMONTH.mat
% plot(G.time,score(:,2));
grid on; axis tight; %xlim([datenum(2002,01,01),datenum(2006,01,01)]);
% legend('EOF1','EOF2');
text(1.05, 0.6,'Breathing','units', ...
    'normalized','color',rgb_x('violet'));
text(1.05, 0.5,'Bifurcation','units', ...
    'normalized','color',rgb_x('teal'));

set(findall(gcf,'-property','fontsize'),'fontsize',8);
set(findall(gcf,'-property','fontname'),'fontname','Latin Modern Roman');


%%
% export_fig /Users/samuelstevens/Dropbox/Hakai/gyre/supply/figs/Freeland.pdf -dpdf -nofontswap

%%


figure
plot(G.time, mode_time_series(:,1))
hold on;
plot(G.time, mode_time_series(:,2))
axdate(20)
yyaxis right
load NEPGpcsMONTH.mat
plot(G.time,score(:,2));
grid on; xlim([datenum(2002,01,01),datenum(2006,01,01)]);
legend('EOF1','EOF2','PCA')

% ylabel('Mode 2 amplitude')
% title('Bifurcation Mode (Mode 2)')

% save('/Users/samuelstevens/Dropbox/Hakai/gyre/supply/data/breathBifurcate.mat',...
%     'G','C','B',"D_CC",'D_GAk','D_NPC','variance_explained','mode_time_series');


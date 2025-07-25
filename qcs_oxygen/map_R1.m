%% Map 
addpath(genpath('/Users/samst/Dropbox/Hakai/'));
addpath(genpath('/Users/samst/Dropbox/UBC/GSW/'));
addpath(genpath('/Users/samst/Dropbox/UBC/m_map/'));
addpath(genpath('/Users/samst/Dropbox/UBC/Misc/'));

clear
load Dcorrected_allProfs_20241011.mat
load glider_allProfs.mat canyonX canyonZ canyonAX
load S2.mat 
load H1.mat
M.Scott2=S2;  M.Hak1=H1;
clearvars H1 S2
load CS09_ctd.mat
load CS01_ctd.mat

msk=CS01.lat>51;
CS01.lon(msk)=NaN;CS01.lat(msk)=NaN; 
CS01.ox(msk)=NaN;

%% find missions from Upwelling season '22 and 23

msk=allProfs.mtime>datenum(2022,04,01) &...
    allProfs.mtime<datenum(2022,11,01) & ...
    allProfs.missionIdx==1;
missions=unique(allProfs.missionNum(msk));


for i=1:length(missions)
    msk=allProfs.missionNum==missions(i);

    lon=allProfs.longitude(msk);
    lat=allProfs.latitude(msk);
    mt=allProfs.mtime(msk);

    mtMsk=mt<mt(lon==min(lon));

    if missions(i)==23
        G.S23.lon=lon(mtMsk);
        G.S23.lat=lat(mtMsk);
    end

    if missions(i)==33
        G.S33.lon=lon(mtMsk);
        G.S33.lat=lat(mtMsk);
    end

    if missions(i)==37
        G.S37.lon=lon(mtMsk);
        G.S37.lat=lat(mtMsk);
    end

    if missions(i)==36
        G.S36.lon=lon(~mtMsk);
        G.S36.lat=lat(~mtMsk);
    end

    if missions(i)==25
        G.S25.lon=lon(~mtMsk);
        G.S25.lat=lat(~mtMsk);
    end


end

% find and plot missions from Upwelling season '23
msk=allProfs.mtime>datenum(2023,04,01) &...
    allProfs.mtime<datenum(2023,11,01) & ...
    allProfs.missionIdx==1;
missions=unique(allProfs.missionNum(msk));

for i=1:length(missions)
    msk=allProfs.missionNum==missions(i);

    lon=allProfs.longitude(msk);
    lat=allProfs.latitude(msk);
    mt=allProfs.mtime(msk);

    mtMsk=mt<mt(lon==min(lon));

    if missions(i)==43
        G.S43.lon=lon(mtMsk);
        G.S43.lat=lat(mtMsk);
    end

    if missions(i)==49
        G.S49.lon=lon(mtMsk);
        G.S49.lat=lat(mtMsk);
    end

    if missions(i)==44
        G.S44.lon=lon(mtMsk);
        G.S44.lat=lat(mtMsk);
    end

    if missions(i)==45
        G.S45.lon=lon(~mtMsk);
        G.S45.lat=lat(~mtMsk);
    end

    if missions(i)==48
        G.S48.lon=lon(~mtMsk);
        G.S48.lat=lat(~mtMsk);
    end

end


%% maps
col=[255 214 140]/255; % YELLOW!
lat_lim=[50 53]; lon_lim=[-132 -127];
fname='/Users/samst/Dropbox/UBC/Misc/british_columbia_3_msl_2013.nc';
lat=ncread(fname,'lat');    
lon=ncread(fname,'lon');
ilon=lon>=lon_lim(1) & lon<=lon_lim(2);
ilat=lat>=lat_lim(1) & lat<=lat_lim(2);
Z=ncread(fname,...
    'Band1',[ find((ilon),1,'first') find((ilat),1,'first')],...
    [ sum(ilon) sum(ilat)],[1 1]);

figure('units','centimeters','outerposition',[0 0 13 17],'color','w');

% ax1=axes('position',[0.525 0.5-0.05 0.425 0.425]);
% t = tiledlayout(2,1,'TileSpacing', 'Compact', 'Padding', 'Compact');

%%%%%%%
axGlob=axes('position',[0.1 0.55 0.75 0.4], 'clipping' , 'on');
% axGlob=nexttile;
% Add small inset map highlighting Queen Charlotte Sound
m_proj('equidistant', 'lon', [120 260], 'lat', [0 60]); % North Pacific bounds
% m_proj('ortho', 'lon', -179, 'lat', 20); 
hold on;

gO=load('globalOxI.mat');
[gO.lat_grid,gO.lon_grid]=meshgrid(gO.lat,gO.lon);
lon_grid_corrected= mod(gO.lon_grid + 360, 360);
% cmap = cmocean('tempo'); cmap(1:30,:)=[];


% Set NaN values to a specific value outside the contour range
data = gO.oxI(:,:,1); 
msk=gO.lon_grid<-100 | gO.lon_grid>120;
data(~msk)=NaN;
data(isnan(data)) = -1;  % Replace NaNs with a value outside the plotted range

% Plot the data with adjusted NaN handling
[CS, CH] = m_contourf(lon_grid_corrected, gO.lat_grid(:,:,1), data, [0:30:360 -1], 'linestyle', 'none');
cmap=(cmocean('tarn','pivotpoint',61));
cmap = [[0.5 0.5 0.5]; cmap];  % Add grey for NaNs at the start
colormap(gca, cmap);


% Set the color axis to include the NaN color
caxis([-2 360]);  % Ensure the colormap includes the NaN value (-1)
% 
m_contour(lon_grid_corrected,gO.lat_grid(:,:,1),gO.geoI,-20:1.25:0,'linecolor',[0.3 0.3 0.3]);

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
    m_text(mean(lonC,'omitnan'),mean(latC,'omitnan'), lbls{count});
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

m_text(247.97,53.24,'Canada','fontweight','bold')
m_text(250.51,42.28,'U.S.A.','fontweight','bold')

% m_scatter(mod(lon_lim(1) + 360, 360),lat_lim(1),100,'ro','linewidth',2);

% Add annotation for Queen Charlotte Sound
% m_text(mean(lon_lim), mean(lat_lim), 'QCS', 'color', 'k', ...
%     'fontsize', 6, 'fontweight', 'bold', 'horizontalalignment', 'center');

% Add grid lines for the orthographic map
m_grid('xtick',4 , 'ytick', 3, 'linestyle', 'none'); % Hide grid lines
text(0.05, 0.95,'a)','fontsize',9,'units','normalized');
%
ax=m_contfbar(1.05,[.2 .8],CS,CH);
ylabel(ax,'O$_2$ ($\mu$mol kg$^{-1}$)','Interpreter','latex');

% ax1=nexttile;
ax1=axes('position',[0.1 0.125 0.8 0.4], 'clipping' , 'on');
m_proj('equidistant','lon',lon_lim,'lat',lat_lim);   % Projection
hold on
deep=-4000:1000:-1000; shelf=[-600:50:-100 0];
[CS,CH]=m_contourf(lon(ilon),lat(ilat),Z',[deep shelf],'linestyle','none');

g=cmocean('gray',600);g(1:390,:)=[];g(end-20:end,:)=[];
cm=[g;m_colmap('blues',50)];
colormap(ax1,cm);

ax = colorbar('location','east outside'); 
set(ax, 'ylim', [-375 0]);
ax.Ticks=[-300  -200  -100     0];
ax.TickLabels={'300';'200';'100';'0'};
xlabel(ax,'Depth (m)','Interpreter','latex')
CH.FaceAlpha=0.9;

[CS,CH]=m_etopo2('contour',[-200 -200],'edgecolor',rgb_x('grey'));
clabel(CS,CH,'labelspacing',300,'color',rgb_x('grey'),'fontsize',6,...
    'fontname','Latin Modern Roman');
m_gshhs_f('patch',col,'edgecolor',[0.3 0.3 0.3],'linewidth',0.1);
m_grid('tickdir','out','linestyle','none','xtick',5,...
    'YaxisLocation','left');

flds=fieldnames(G);
lnes=linspecer(length(flds));
for i=1:length(flds)
    m_plot(G.(flds{i}).lon,G.(flds{i}).lat,'color',lnes(i,:),'linewidth',0.5);
end

flds=fieldnames(M);
lnes=[rgb_x('steel blue');rgb_x('rose')];
for i=1:length(flds)
    m_scatter(M.(flds{i}).longitude,M.(flds{i}).latitude,150,'^',...
        'markeredgecolor','w','markerfacecolor',lnes(i,:),'linewidth',2);
end
% m_text(M.(flds{2}).longitude-.1,M.(flds{2}).latitude,'Hak1','color',lnes(2,:));
% m_text(M.(flds{1}).longitude-.1,M.(flds{1}).latitude,'Scott2','color',lnes(1,:));
% keyboard
m_annotation('textarrow',[-129.0631-1.5  -129.0631-0.25],[51.4114-0.3 51.4114-0.1],...
    'string',{'Sea Otter';'Trough'},'fontsize',6,'Fontweight','normal',...
    'horizontalalignment','right','headlength',5);
m_annotation('textarrow',[-129.2  -128.7],[50.65 51.1],...
    'string',{'Cook Trough'},'fontsize',6,'Fontweight','normal',...
    'horizontalalignment','right','headlength',5);
m_annotation('textarrow',[-129.0631-1.5  -129.9253],[51.4114 51.6550-0.1],...
    'string',{'Mitchell';'Trough'},'fontsize',6,'Fontweight','normal',...
    'horizontalalignment','right','headlength',5);

m_scatter(mean(CS09.lon),mean(CS09.lat),50,'filled','markerfacecolor',rgb_x('black'),...
    'markeredgecolor',rgb_x('white'),'markerfacealpha',0.8);
m_scatter(nanmean(CS01.lon),nanmean(CS01.lat),50,'filled','markerfacecolor',rgb_x('red'),...
    'markeredgecolor',rgb_x('white'),'markerfacealpha',0.8);
% m_text(nanmean(CS09.lon)-.1,nanmean(CS09.lat),'CS09','color','k');
% m_text(nanmean(CS01.lon)-.1,nanmean(CS01.lat),'CS01','color','r');

m_text(-127.35,51.99,{'BC'},'fontsize',6,'horizontalalignment','center');
text(0.05, 0.95, 'b)', 'Units', 'normalized');

set(findall(gcf,'-property','fontsize'),'fontsize',8);
set(findall(gcf,'-property','fontname'),'fontname','Latin Modern Roman');

%%
export_fig /Users/samuelstevens/Dropbox/Hakai/gliders/oxygen/figures/map_R1.pdf -dpdf -nofontswap
% exportgraphics(gcf, '/Users/samuelstevens/Dropbox/Hakai/gliders/oxygen/figures/historicalHR.pdf', ...
%     'ContentType', 'vector', 'Resolution', 300);

% export_fig /Users/samuelstevens/Dropbox/Hakai/gliders/oxygen/figures/poster/map.png -m6 -nofontswap -transparent

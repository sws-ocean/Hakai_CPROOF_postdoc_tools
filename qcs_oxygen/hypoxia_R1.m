%% Oxygen story
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
% load DO_1990-2024_QCS.mat
% load DO_1990-2024_QCS_v2.mat
load isoDO_DFOandHakai.mat
load CS09_ctd.mat
load QCSwind.mat

bottomDO.oxygen(bottomDO.oxygen>600)=NaN;

sig_thresh = [26 26.5 26.75];

allProfs.s_t=gsw_sigma0(allProfs.SA,allProfs.CT);
allProfs.spice=gsw_spiciness0(allProfs.SA,allProfs.CT);
%% find and plot missions from Upwelling season '22 and 23

msk=allProfs.mtime>datenum(2022,04,01) &...
    allProfs.mtime<datenum(2022,11,01) & ...
    allProfs.missionIdx==1;
missions=unique(allProfs.missionNum(msk));

smoothing={1 1};

for i=1:length(missions)
    msk=allProfs.missionNum==missions(i);

    T=allProfs.temperature(:,msk);
    DO=allProfs.oxygen_corrected(:,msk);
    s_t=allProfs.s_t(:,msk);
    X=allProfs.canyonX(msk);
    lon=allProfs.longitude(msk);
    lat=allProfs.latitude(msk);
    mt=allProfs.mtime(msk);
    idx=allProfs.canyonIDX(msk);
    Mis=allProfs.mission{msk};
    
    mtMsk=mt<mt(lon==min(lon));

    [xx,yy]=meshgrid(X(mtMsk),allProfs.z(:,1));
    [xi,yi]=meshgrid(0:0.5:max(canyonX),allProfs.z(:,1));
    zi=griddata(xx,yy,DO(:,mtMsk),xi,yi);
    B = smoothdata2(zi,"movmean",smoothing);
    zi=griddata(xx,yy,s_t(:,mtMsk),xi,yi);
    Bs = smoothdata2(zi,"movmean",smoothing);

    if missions(i)==23
       G.S23.xi=xi; G.S23.yi=yi; G.S23.B=B; G.S23.Bs=Bs;
       G.S23.tt=[min(mt(mtMsk)) max(mt(mtMsk))];
       G.S23.t=[datestr(min(mt(mtMsk)),'dd mmm yyyy'),'-',...
        datestr(max(mt(mtMsk)),'dd mmm yyyy')];
       G.S23.lon=lon(mtMsk);
       G.S23.lat=lat(mtMsk);

       count=0;
       idx=find(mtMsk);
       G.S23.BTMox=[];
       G.S23.BTMst=[];
       for ii=1:length(idx)
           if X(idx(ii))<118
               count=count+1;
               msk=find(~isnan(DO(:,idx(ii))),1,'last');
               G.S23.BTMox(ii)=DO(msk,idx(ii));
               G.S23.BTMst(ii)=s_t(msk,idx(ii));
               G.S23.BTMidx(ii)=msk;
           end
       end
    end

    if missions(i)==33
       G.S33.xi=xi; G.S33.yi=yi; G.S33.B=B; G.S33.Bs=Bs;
       G.S33.tt=[min(mt(mtMsk)) max(mt(mtMsk))];
       G.S33.t=[datestr(min(mt(mtMsk)),'dd mmm yyyy'),'-',...
        datestr(max(mt(mtMsk)),'dd mmm yyyy')];
       G.S33.lon=lon(mtMsk);
       G.S33.lat=lat(mtMsk);

       count=0;
       idx=find(mtMsk);
       G.S33.BTMox=[];
       G.S33.BTMst=[];
       for ii=1:length(idx)
           if X(idx(ii))<118 && X(idx(ii))>0
               count=count+1;
               msk=find(~isnan(DO(:,idx(ii))),1,'last');
               G.S33.BTMox(ii)=DO(msk,idx(ii));
               G.S33.BTMst(ii)=s_t(msk,idx(ii));
               G.S33.BTMidx(ii)=msk;
           end
       end
    end

     if missions(i)==37
       G.S37.xi=xi; G.S37.yi=yi; G.S37.B=B; G.S37.Bs=Bs;
       G.S37.tt=[min(mt(mtMsk)) max(mt(mtMsk))];
       G.S37.t=[datestr(min(mt(mtMsk)),'dd mmm yyyy'),'-',...
        datestr(max(mt(mtMsk)),'dd mmm yyyy')];
       G.S37.lon=lon(mtMsk);
       G.S37.lat=lat(mtMsk);

       count=0;
       idx=find(mtMsk);
       G.S37.BTMox=[];
       G.S37.BTMst=[];
       for ii=1:length(idx)
           if X(idx(ii))<118
               count=count+1;
               msk=find(~isnan(DO(:,idx(ii))),1,'last');
               G.S37.BTMox(ii)=DO(msk,idx(ii));
               G.S37.BTMst(ii)=s_t(msk,idx(ii));
               G.S37.BTMidx(ii)=msk;
           end
       end
     end

    [xx,yy]=meshgrid(X(~mtMsk),allProfs.z(:,1));
    [xi,yi]=meshgrid(0:0.5:max(canyonX),allProfs.z(:,1));
    zi=griddata(xx,yy,DO(:,~mtMsk),xi,yi);
    B = smoothdata2(zi,"movmean",smoothing);
    zi=griddata(xx,yy,s_t(:,~mtMsk),xi,yi);
    Bs = smoothdata2(zi,"movmean",smoothing);

    if missions(i)==36
       G.S36.xi=xi; G.S36.yi=yi; G.S36.B=B; G.S36.Bs=Bs;
       G.S36.tt=[min(mt(~mtMsk)) max(mt(~mtMsk))];
       G.S36.t=[datestr(min(mt(~mtMsk)),'dd mmm yyyy'),'-',...
        datestr(max(mt(~mtMsk)),'dd mmm yyyy')];
       G.S36.lon=lon(~mtMsk);
       G.S36.lat=lat(~mtMsk);

       count=0;
       idx=find(~mtMsk);
       G.S36.BTMox=[];
       G.S36.BTMst=[];
       for ii=1:length(idx)
           if X(idx(ii))<118 && sum(isfinite(DO(:,idx(ii))))>0
               count=count+1;
               msk=find(~isnan(DO(:,idx(ii))),1,'last');
               G.S36.BTMox(ii)=DO(msk,idx(ii));
               G.S36.BTMst(ii)=s_t(msk,idx(ii));
               G.S36.BTMidx(ii)=msk;
           end
       end
    end

    if missions(i)==25
       G.S25.xi=xi; G.S25.yi=yi; G.S25.B=B; G.S25.Bs=Bs;
       G.S25.tt=[min(mt(~mtMsk)) max(mt(~mtMsk))];
       G.S25.t=[datestr(min(mt(~mtMsk)),'dd mmm yyyy'),'-',...
        datestr(max(mt(~mtMsk)),'dd mmm yyyy')];
       G.S25.lon=lon(~mtMsk);
       G.S25.lat=lat(~mtMsk);

       count=0;
       idx=find(~mtMsk);
       G.S25.BTMox=[];
       G.S25.BTMst=[];
       for ii=1:length(idx)
           if X(idx(ii))<118 && sum(isfinite(DO(:,idx(ii))))>0
               count=count+1;
               msk=find(~isnan(DO(:,idx(ii))),1,'last');
               G.S25.BTMox(ii)=DO(msk,idx(ii));
               G.S25.BTMst(ii)=s_t(msk,idx(ii));
               G.S25.BTMidx(ii)=msk;
           end
       end
    end


end

% find and plot missions from Upwelling season '23
msk=allProfs.mtime>datenum(2023,04,01) &...
    allProfs.mtime<datenum(2023,11,01) & ...
    allProfs.missionIdx==1;
missions=unique(allProfs.missionNum(msk));

for i=1:length(missions)
    msk=allProfs.missionNum==missions(i);

    T=allProfs.temperature(:,msk);
    DO=allProfs.oxygen_corrected(:,msk);
    s_t=allProfs.s_t(:,msk);
    X=allProfs.canyonX(msk);
    lon=allProfs.longitude(msk);
    lat=allProfs.latitude(msk);
    mt=allProfs.mtime(msk);
    idx=allProfs.canyonIDX(msk);
    Mis=allProfs.mission{msk};

    mtMsk=mt<mt(lon==min(lon));

    [xx,yy]=meshgrid(X(mtMsk),allProfs.z(:,1));
    [xi,yi]=meshgrid(0:0.5:max(canyonX),allProfs.z(:,1));
    zi=griddata(xx,yy,DO(:,mtMsk),xi,yi);
    B = smoothdata2(zi,"movmean",smoothing);
    zi=griddata(xx,yy,s_t(:,mtMsk),xi,yi);
    Bs = smoothdata2(zi,"movmean",smoothing);

    if missions(i)==43
        G.S43.xi=xi; G.S43.yi=yi; G.S43.B=B; G.S43.Bs=Bs;
        G.S43.tt=[min(mt(mtMsk)) max(mt(mtMsk))];
        G.S43.t=[datestr(min(mt(mtMsk)),'dd mmm yyyy'),'-',...
            datestr(max(mt(mtMsk)),'dd mmm yyyy')];
        G.S43.lon=lon(mtMsk);
        G.S43.lat=lat(mtMsk);

        count=0;
        idx=find(mtMsk);
        G.S43.BTMox=[];
        G.S43.BTMst=[];
        for ii=1:length(idx)
            if X(idx(ii))<118
                count=count+1;
                msk=find(~isnan(DO(:,idx(ii))),1,'last');
                G.S43.BTMox(ii)=DO(msk,idx(ii));
                G.S43.BTMst(ii)=s_t(msk,idx(ii));
            end
        end
    end

    if missions(i)==49
        G.S49.xi=xi; G.S49.yi=yi; G.S49.B=B; G.S49.Bs=Bs;
        G.S49.tt=[min(mt(mtMsk)) max(mt(mtMsk))];
        G.S49.t=[datestr(min(mt(mtMsk)),'dd mmm yyyy'),'-',...
            datestr(max(mt(mtMsk)),'dd mmm yyyy')];
        G.S49.lon=lon(mtMsk);
        G.S49.lat=lat(mtMsk);

        count=0;
        idx=find(mtMsk);
        G.S49.BTMox=[];
        G.S49.BTMst=[];
        for ii=1:length(idx)
            if X(idx(ii))<118
                count=count+1;
                msk=find(~isnan(DO(:,idx(ii))),1,'last');
                G.S49.BTMox(ii)=DO(msk,idx(ii));
                G.S49.BTMst(ii)=s_t(msk,idx(ii));
            end
        end
    end

    if missions(i)==44
        G.S44.xi=xi; G.S44.yi=yi; G.S44.B=B; G.S44.Bs=Bs;
        G.S44.tt=[min(mt(mtMsk)) max(mt(mtMsk))];
        G.S44.t=[datestr(min(mt(mtMsk)),'dd mmm yyyy'),'-',...
            datestr(max(mt(mtMsk)),'dd mmm yyyy')];
        G.S44.lon=lon(mtMsk);
        G.S44.lat=lat(mtMsk);

        count=0;
        idx=find(mtMsk);
        G.S44.BTMox=[];
        G.S44.BTMst=[];
        for ii=1:length(idx)
            if X(idx(ii))<118
                count=count+1;
                msk=find(~isnan(DO(:,idx(ii))),1,'last');
                G.S44.BTMox(ii)=DO(msk,idx(ii));
                G.S44.BTMst(ii)=s_t(msk,idx(ii));
            end
        end
    end

     [xx,yy]=meshgrid(X(~mtMsk),allProfs.z(:,1));
    [xi,yi]=meshgrid(0:0.5:max(canyonX),allProfs.z(:,1));
    zi=griddata(xx,yy,DO(:,~mtMsk),xi,yi);
    B = smoothdata2(zi,"movmean",smoothing);
    zi=griddata(xx,yy,s_t(:,~mtMsk),xi,yi);
    Bs = smoothdata2(zi,"movmean",smoothing);

    if missions(i)==45
       G.S45.xi=xi; G.S45.yi=yi; G.S45.B=B; G.S45.Bs=Bs;
       G.S45.tt=[min(mt(~mtMsk)) max(mt(~mtMsk))];
       G.S45.t=[datestr(min(mt(~mtMsk)),'dd mmm yyyy'),'-',...
        datestr(max(mt(~mtMsk)),'dd mmm yyyy')];
       G.S45.lon=lon(~mtMsk);
       G.S45.lat=lat(~mtMsk);

       count=0;
       idx=find(~mtMsk);
       G.S45.BTMox=[];
       G.S45.BTMst=[];
       for ii=1:length(idx)
           if X(idx(ii))<118 && sum(isfinite(DO(:,idx(ii))))>0
               count=count+1;
               msk=find(~isnan(DO(:,idx(ii))),1,'last');
               G.S45.BTMox(ii)=DO(msk,idx(ii));
               G.S45.BTMst(ii)=s_t(msk,idx(ii));
           end
       end
    end

    if missions(i)==48
       G.S48.xi=xi; G.S48.yi=yi; G.S48.B=B; G.S48.Bs=Bs;
       G.S48.tt=[min(mt(~mtMsk)) max(mt(~mtMsk))];
       G.S48.t=[datestr(min(mt(~mtMsk)),'dd mmm yyyy'),'-',...
        datestr(max(mt(~mtMsk)),'dd mmm yyyy')];
       G.S48.lon=lon(~mtMsk);
       G.S48.lat=lat(~mtMsk);

       count=0;
       idx=find(~mtMsk);
       G.S48.BTMox=[];
       G.S48.BTMst=[];
       for ii=1:length(idx)
           if X(idx(ii))<118 && sum(isfinite(DO(:,idx(ii))))>0
               count=count+1;
               msk=find(~isnan(DO(:,idx(ii))),1,'last');
               G.S48.BTMox(ii)=DO(msk,idx(ii));
               G.S48.BTMst(ii)=s_t(msk,idx(ii));
           end
       end
    end

end


%% Make glider Ob time-series near moorings
Gwnd=7*1; % glider window average in days
allProfs.Hak1.time=[ceil(min(M.Hak1.time)):Gwnd:floor(floor(max(M.Hak1.time)))]';

allProfs.Hak1.temperature=NaN(size(allProfs.Hak1.time));
allProfs.Hak1.salinity=NaN(size(allProfs.Hak1.time));
allProfs.Hak1.oxygen=NaN(size(allProfs.Hak1.time));
allProfs.Hak1.sigma_theta=NaN(size(allProfs.Hak1.time));

allProfs.Hak1.dist=gsw_distance([ones(size(allProfs.longitude))'.*M.Hak1.longitude ...
    allProfs.longitude'],[ones(size(allProfs.longitude))'.*M.Hak1.latitude ...
    allProfs.latitude']);

h=waitbar(0,'Please wait...');
for i=1:length(allProfs.Hak1.time)-1
    waitbar(i / (length(allProfs.Hak1.time) - 1), h);

    msk=allProfs.Hak1.dist' < 2.5e3 & ...
        allProfs.mtime>allProfs.Hak1.time(i) &...
        allProfs.mtime<allProfs.Hak1.time(i)+Gwnd;

    if sum(msk)>0
        allProfs.Hak1.temperature(i)=mean(mean(allProfs.temperature(112:end,msk),1,'omitnan'),'omitnan');
        allProfs.Hak1.salinity(i)=mean(mean(allProfs.SA(112:end,msk),1,'omitnan'),'omitnan');
        allProfs.Hak1.oxygen(i)=mean(min(allProfs.oxygen_corrected(112:end,msk)),'omitnan');
        allProfs.Hak1.sigma_theta(i)=mean(mean(allProfs.s_t(112:end,msk),1,'omitnan'),'omitnan');
    end
end
delete(h);


allProfs.Scott2.time=[ceil(min(M.Scott2.time)):Gwnd:floor(floor(max(M.Scott2.time)))]';

allProfs.Scott2.temperature=NaN(size(allProfs.Scott2.time));
allProfs.Scott2.salinity=NaN(size(allProfs.Scott2.time));
allProfs.Scott2.oxygen=NaN(size(allProfs.Scott2.time));
allProfs.Scott2.sigma_theta=NaN(size(allProfs.Scott2.time));

allProfs.Scott2.dist=gsw_distance([ones(size(allProfs.longitude))'.*M.Scott2.longitude ...
    allProfs.longitude'],[ones(size(allProfs.longitude))'.*M.Scott2.latitude ...
    allProfs.latitude']);

h=waitbar(0,'Please wait...');
for i=1:length(allProfs.Scott2.time)-1
    waitbar(i / (length(allProfs.Scott2.time) - 1), h);
    msk=allProfs.Scott2.dist' < 15e3 & ...
    allProfs.mtime>allProfs.Scott2.time(i) &...
    allProfs.mtime<allProfs.Scott2.time(i)+Gwnd;

    if sum(msk)>0
        allProfs.Scott2.temperature(i)=mean(mean(allProfs.temperature(259:end,msk),1,'omitnan'),'omitnan');
        allProfs.Scott2.salinity(i)=mean(mean(allProfs.SA(259:end,msk),1,'omitnan'),'omitnan');
        allProfs.Scott2.oxygen(i)=mean(min(allProfs.oxygen_corrected(259:end,msk)),'omitnan');
        allProfs.Scott2.sigma_theta(i)=mean(mean(allProfs.s_t(259:end,msk),1,'omitnan'),'omitnan');
    end
end
delete(h);

% [tt,zz]=meshgrid(allProfs.mtime(msk),allProfs.z(:,1));

% Find moorings on transect
dd=gsw_distance([ones(size(canyonAX(:,1))).*M.Hak1.longitude canyonAX(:,1)],...
    [ones(size(canyonAX(:,1))).*M.Hak1.latitude canyonAX(:,2)]);
[~,M.Hak1.canyonIDX]=min(dd);

dd=gsw_distance([ones(size(canyonAX(:,1))).*M.Scott2.longitude canyonAX(:,1)],...
    [ones(size(canyonAX(:,1))).*M.Scott2.latitude canyonAX(:,2)]);
[~,M.Scott2.canyonIDX]=min(dd);
%% Interpolate in x

flds=fieldnames(G);

for i=1:length(flds)

    % Replace NaNs with a value outside the data range (e.g., -1)
    G.(flds{i}).Cdata = G.(flds{i}).B;
    % G.(flds{i}).Cdata=inpaint_nans(G.(flds{i}).Cdata);
    G.(flds{i}).Cdata(isnan(G.(flds{i}).Cdata)) = -1;

end

% 
% for i=1:length(flds)
%     for ii=1:285
%         msk1=find(~isnan(G.(flds{i}).B(ii,:)),1,'first');
%         msk2=find(isnan(G.(flds{i}).B(ii,:)),1,'last');
%         if ~isempty(msk2)
%             G.(flds{i}).B(ii,msk1:msk2)=inpaint_nans(G.(flds{i}).B(ii,msk1:msk2));
%         end
%         G.(flds{i}).Cdata=G.(flds{i}).B;
%         msk1=find(~isnan(G.(flds{i}).Bs(ii,:)),1,'first');
%         msk2=find(isnan(G.(flds{i}).Bs(ii,:)),1,'last');
%         if ~isempty(msk2)
%             G.(flds{i}).Bs(ii,msk1:msk2)=inpaint_nans(G.(flds{i}).Bs(ii,msk1:msk2));
%         end
%     end
% end


% for i=1:length(flds)
%     for ii=1:400
%         G.(flds{i}).Cdata(ii,:)=inpaint_nans(G.(flds{i}).B(ii,:));
%     end
% end

%% make pretty plot
figure('units','centimeters','outerposition',[0 0 18 17],'color','w');

% [XX,~]=meshgrid(128:-0.5:-30,1:1100);
% canyonXX=128.1:-0.1:-30;
XX=(xi-120).*-1;
canyonXX=(canyonX-120).*-1;

axMap=axes('position', [0.11 0.85 0.7 0.15-0.0125]);
load('BCcoast');
lat_lim=[50.75 52]; lon_lim=[-130.5 -127.5];

fname='/Users/samst/Dropbox/UBC/Misc/british_columbia_3_msl_2013.nc';
lat=ncread(fname,'lat');    
lon=ncread(fname,'lon');
ilon=lon>=lon_lim(1) & lon<=lon_lim(2);
ilat=lat>=lat_lim(1) & lat<=lat_lim(2);
Z=ncread(fname,...
    'Band1',[ find((ilon),1,'first') find((ilat),1,'first')],...
    [ sum(ilon) sum(ilat)],[1 1]);
col=[255 214 140]/255; % YELLOW!
m_proj('equidistant','lon',lon_lim,'lat',lat_lim);   % Projection
hold on
deep=-4000:1000:-1000; shelf=[-600:50:-100 0];
[CS,CH]=m_contourf(lon(ilon),lat(ilat),Z',[deep shelf],'linestyle','none');
g=cmocean('gray',600);g(1:390,:)=[];g(end-20:end,:)=[];
cm=[g;m_colmap('blues',50)];
colormap(cm);
ax = colorbar('Position', [0.9 0.8 0.01 0.1]); 
set(ax, 'ylim', [-375 0]);
ax.Ticks=[-300  -200  -100     0];
ax.TickLabels={'300';'200';'100';'0'};
title(ax,'Depth (m)','Interpreter','latex')
CH.FaceAlpha=0.9;
m_gshhs_f('patch',col,'edgecolor',[0.3 0.3 0.3],'linewidth',0.3);
m_grid('linestyle','none','linewidth',0.5,'tickdir','out','xtick',3);
m_plot(canyonAX(:,1),canyonAX(:,2),'k','linewidth',1);

% Build axis
lnes=linspecer(length(0:25:max(canyonX)));
count=0;

for i=0:30:max(canyonX)
    [~,idx]=min(abs(canyonX-i));
    count=count+1;
    m_scatter(canyonAX(idx,1),canyonAX(idx,2),30,'filled','markerfacecolor','k',...
        'markeredgecolor','w','marker','s');
    % m_text(canyonAX(idx,1)-0.4,canyonAX(idx,2),sprintf('%3.0f km',canyonXX(idx)));

end
m_scatter(w22.lon,w22.lat,50,'filled','markeredgecolor','w', ...
    'markerfacecolor',rgb_x('rose'),'linewidth',1,'markerfacealpha',0.5);
text(0.05, 0.9, 'a)', 'Units', 'normalized','fontweight','bold');

%%%%%%%%%%%%%%
lnes=[rgb_x('steel blue');rgb_x('rose')];
axG1 = axes('color', [0.75 0.75 0.75], 'position', [0.1 0.75-0.075 0.35 0.15-0.0125]);
Cdata = G.S23.Cdata;
Cdata(isnan(Cdata)) = -1;
[CC, hh] = contourf(XX, G.S23.yi, Cdata, [0:20:450, -1], 'LineStyle', 'none');
clim([0 400]);
cmap = cmocean('tarn', 'pivotpoint', 1.4*43.7);
cmap = [[0.9 0.9 0.9]; cmap];  % Add a color for NaNs at the beginning
colormap(gca, cmap);
axis ij;
ylim([0 400]);
hold on
contour(XX, G.S23.yi, G.S23.B, [1.4*43.7 1.4*43.7], 'linewidth', 0.75, ...
    'linecolor', [0.1 0.1 0.1], 'linestyle', '-');
l = area(canyonXX, canyonZ*-1, 1000, 'facecolor',[0.8 0.8 0.8], 'facealpha', 1, ...
    'edgecolor', 'none');
[CS, CH] = contour(XX, G.S23.yi, G.S23.Bs, [26 26.5 26.75], 'linewidth', 0.5, ...
    'linecolor', [0.5 0.5 0.5]);
[CS2, CH2] = contour(XX, G.S23.yi, G.S23.Bs, [26.65 26.65], 'linewidth', 0.5, ...
    'linecolor', [0.7 0.7 0.7]);
clabel(CS2, CH2, 'LabelSpacing', 200, 'color', [0.6 0.6 0.6], 'fontsize', 5, ...
    'fontname', 'Latin Modern Roman');
clabel(CS, CH, 'LabelSpacing', 200, 'color', [0.4 0.4 0.4], 'fontsize', 6, ...
    'fontname', 'Latin Modern Roman');
yticks(0:150:450); xticks(-30:30:120);
xticklabels([]);
set(gca, 'TickDir', 'out');
cl = clim;
text(10,350, G.S23.t);

CB=m_contfbar([0.05 0.6],1.1,CC,hh,'axfrac',0.07,'edgecolor','none');
yl=CB.YLim; axes(CB);
line([1.4*43.7 1.4*43.7],yl,'color','k','linewidth',1,'linestyle','-');
CB.XTick=[1 61 150 250 400]; CB.XTickLabelRotation=0;
title(CB,{'O$_2$ ($\mu$mol kg$^{-1}$)'},'Interpreter','latex');

msk=datenum(w22.t)>G.S23.tt(1) & datenum(w22.t)<G.S23.tt(2);
dir  = w22.dir(msk);            % meteorological ° (0 =N, 90 =E, …)
spd  = w22.speed(msk);          % m s-¹  (optional, not used below)
frac   = 0.40;          % inset side = 30 % of parent-axes width/height
pad    = 0.01;          % 1 % figure padding from the edges of axG1
posG1  = get(axG1,'Position');            % [left bottom width height]
rosePos = [ posG1(1) + posG1(3) - posG1(3)*frac - pad , ... % left
            posG1(2)                           + pad , ... % bottom
            posG1(3)*frac , ...                           % width
            posG1(4)*frac ];                              % height
rosePos(1)=rosePos(1)+0.047;

pax1 = polaraxes('Position',rosePos,'Color','none');
hold(pax1,'on')
nbins = 16;
polarhistogram(pax1,deg2rad(dir),nbins, ...
               'Normalization','probability', ...
               'FaceColor',rgb_x('rose'),'EdgeColor','none');
pax1.ThetaZeroLocation = 'top';     % 0° = North
pax1.ThetaDir          = 'clockwise';
pax1.RLim              = [0 0.25];  % shrink radial scale
pax1.RTick             = [];
pax1.GridAlpha         = 0.75;
pax1.FontSize          = 6;
pax1.ThetaTickLabel    = [];        % hide azimuth labels


%%%%%%%%%%%%
axG2=axes('color',[0.75 0.75 0.75],'position',[0.1 0.6-0.075 0.35 0.15-0.0125]);
[~,h]=contourf(XX,G.S25.yi,G.S25.Cdata,[0:20:450, -1],'LineStyle','none');
clim([0 400]);
% Set the colormap and add a transparent color for the NaNs
cmap = cmocean('tarn', 'pivotpoint', 1.4*43.7);
cmap = [[0.9 0.9 0.9]; cmap];  % Add a color for NaNs at the beginning
colormap(gca, cmap);
axis ij; ylim([0 400]);
hold on
contour(XX,G.S25.yi,G.S25.B,[1.4*43.7 1.4*43.7],'linewidth',0.75,...
    'linecolor',[0.1 0.1 0.1],'linestyle','-');
l=area(canyonXX,canyonZ*-1,1000,'facecolor',[0.8 0.8 0.8],'facealpha',1,...
    'edgecolor','none');
[CS,CH]=contour(XX,G.S25.yi,G.S25.Bs,[26 26.5 26.75],'linewidth',0.5,...
    'linecolor',[0.5 0.5 0.5]);
[CS2,CH2]=contour(XX,G.S25.yi,G.S25.Bs,[26.65 26.65],'linewidth',0.5,...
    'linecolor',[0.7 0.7 0.7]);
clabel(CS2,CH2,'LabelSpacing',200,'color',[0.6 0.6 0.6],'fontsize',5,...
    'fontname','Latin Modern Roman');
clabel(CS,CH,'LabelSpacing',200,'color',[0.4 0.4 0.4],'fontsize',6,...
    'fontname','Latin Modern Roman');
yticks(0:150:450);xticks(-30:30:120); xticklabels([]); set(gca, 'TickDir', 'out')
cl=clim;
text(10,350,G.S25.t);

msk=datenum(w22.t)>G.S25.tt(1) & datenum(w22.t)<G.S25.tt(2);
dir  = w22.dir(msk);            % meteorological ° (0 =N, 90 =E, …)
spd  = w22.speed(msk);          % m s-¹  (optional, not used below)
frac   = 0.40;          % inset side = 30 % of parent-axes width/height
pad    = 0.01;          % 1 % figure padding from the edges of axG1
posG1  = get(axG2,'Position');            % [left bottom width height]
rosePos = [ posG1(1) + posG1(3) - posG1(3)*frac - pad , ... % left
            posG1(2)                           + pad , ... % bottom
            posG1(3)*frac , ...                           % width
            posG1(4)*frac ];                              % height
rosePos(1)=rosePos(1)+0.047;

pax2 = polaraxes('Position',rosePos,'Color','none');
hold(pax2,'on')
nbins = 16;
polarhistogram(pax2,deg2rad(dir),nbins, ...
               'Normalization','probability', ...
               'FaceColor',rgb_x('rose'),'EdgeColor','none');
pax2.ThetaZeroLocation = 'top';     % 0° = North
pax2.ThetaDir          = 'clockwise';
pax2.RLim              = [0 0.25];  % shrink radial scale
pax2.RTick             = [];
pax2.GridAlpha         = 0.75;
pax2.FontSize          = 6;
pax2.ThetaTickLabel    = [];        % hide azimuth labels

%%%%%%%%%%%%
axG3=axes('color',[0.75 0.75 0.75],'position',[0.1 0.45-0.075 0.35 0.15-0.0125]);
axG3.Color = [0.6 0.6 0.6];
[CS,CH]=contourf(XX,G.S33.yi,G.S33.Cdata,[0:20:450, -1],'LineStyle','none');
clim([0 400]);
cmap = cmocean('tarn', 'pivotpoint', 1.4*43.7);
cmap = [[0.9 0.9 0.9]; cmap];
colormap(axG3, cmap);
axis ij; ylim([0 400]);
hold on
contour(XX,G.S33.yi,G.S33.B,[1.4*43.7 1.4*43.7],'linewidth',0.75,...
    'linecolor',[0.1 0.1 0.1],'linestyle','-');
l=area(canyonXX,canyonZ*-1,1000,'facecolor',[0.8 0.8 0.8],'facealpha',1,...
    'edgecolor','none');
contour(XX,G.S33.yi,G.S33.Bs,[26 26.5 26.75],'linewidth',0.5,...
    'linecolor',[0.5 0.5 0.5]);
[CS2,CH2]=contour(XX,G.S33.yi,G.S33.Bs,[26.65 26.65],'linewidth',0.5,...
    'linecolor',[0.7 0.7 0.7]);
clabel(CS2,CH2,'LabelSpacing',200,'color',[0.6 0.6 0.6],'fontsize',5,...
    'fontname','Latin Modern Roman');
clabel(CS,CH,'LabelSpacing',200,'color',[0.4 0.4 0.4],'fontsize',6,...
    'fontname','Latin Modern Roman');
yticks(0:150:450); xticks(-30:30:120); xticklabels([]); set(gca, 'TickDir', 'out')
cl=clim;
text(10,350,G.S33.t);
ylabel('Depth (m)');

msk=datenum(w22.t)>G.S33.tt(1) & datenum(w22.t)<G.S33.tt(2);
dir  = w22.dir(msk);            % meteorological ° (0 =N, 90 =E, …)
spd  = w22.speed(msk);          % m s-¹  (optional, not used below)
frac   = 0.40;          % inset side = 30 % of parent-axes width/height
pad    = 0.01;          % 1 % figure padding from the edges of axG1
posG1  = get(axG3,'Position');            % [left bottom width height]
rosePos = [ posG1(1) + posG1(3) - posG1(3)*frac - pad , ... % left
            posG1(2)                           + pad , ... % bottom
            posG1(3)*frac , ...                           % width
            posG1(4)*frac ];                              % height
rosePos(1)=rosePos(1)+0.047;

pax3 = polaraxes('Position',rosePos,'Color','none');
hold(pax3,'on')
nbins = 16;
polarhistogram(pax3,deg2rad(dir),nbins, ...
               'Normalization','probability', ...
               'FaceColor',rgb_x('rose'),'EdgeColor','none');
pax3.ThetaZeroLocation = 'top';     % 0° = North
pax3.ThetaDir          = 'clockwise';
pax3.RLim              = [0 0.25];  % shrink radial scale
pax3.RTick             = [];
pax3.GridAlpha         = 0.75;
pax3.FontSize          = 6;
pax3.ThetaTickLabel    = [];        % hide azimuth labels

%%%%%%%%%%%%
axG4=axes('color',[0.75 0.75 0.75],'position',[0.1 0.15-0.075 0.35 0.15-0.015]);
[CS,CH]=contourf(XX,G.S36.yi,G.S36.Cdata,[0:20:450, -1],'LineStyle','none');
clim(cl);
% Set the colormap and add a transparent color for the NaNs
cmap = cmocean('tarn', 'pivotpoint', 1.4*43.7);
cmap = [[0.9 0.9 0.9]; cmap];  % Add a color for NaNs at the beginning
colormap(gca, cmap);
axis ij; ylim([0 400]);
hold on
contour(XX,G.S36.yi,G.S36.B,[1.4*43.7 1.4*43.7],'linewidth',0.75,...
    'linecolor',[0.1 0.1 0.1],'linestyle','-');
l=area(canyonXX,canyonZ*-1,1000,'facecolor',[0.8 0.8 0.8],'facealpha',1,...
    'edgecolor','none');
[CSc,CHc]=contour(XX,G.S36.yi,G.S36.Bs,[26 26.5 26.75],'linewidth',0.5,...
    'linecolor',[0.5 0.5 0.5]);
[CS2,CH2]=contour(XX,G.S36.yi,G.S36.Bs,[26.65 26.65],'linewidth',0.5,...
    'linecolor',[0.7 0.7 0.7]);
clabel(CS2,CH2,'LabelSpacing',200,'color',[0.6 0.6 0.6],'fontsize',5,...
    'fontname','Latin Modern Roman');
clabel(CSc,CHc,'LabelSpacing',200,'color',[0.4 0.4 0.4],'fontsize',6,...
    'fontname','Latin Modern Roman');
yticks(0:150:450); xticks(-30:30:120); set(gca, 'TickDir', 'out')
text(10,350,G.S36.t);

msk=datenum(w22.t)>G.S36.tt(1) & datenum(w22.t)<G.S36.tt(2);
dir  = w22.dir(msk);            % meteorological ° (0 =N, 90 =E, …)
spd  = w22.speed(msk);          % m s-¹  (optional, not used below)
frac   = 0.40;          % inset side = 30 % of parent-axes width/height
pad    = 0.01;          % 1 % figure padding from the edges of axG1
posG1  = get(axG4,'Position');            % [left bottom width height]
rosePos = [ posG1(1) + posG1(3) - posG1(3)*frac - pad , ... % left
            posG1(2)                           + pad , ... % bottom
            posG1(3)*frac , ...                           % width
            posG1(4)*frac ];                              % height
rosePos(1)=rosePos(1)+0.047;

pax4 = polaraxes('Position',rosePos,'Color','none');
hold(pax4,'on')
nbins = 16;
polarhistogram(pax4,deg2rad(dir),nbins, ...
               'Normalization','probability', ...
               'FaceColor',rgb_x('rose'),'EdgeColor','none');
pax4.ThetaZeroLocation = 'top';     % 0° = North
pax4.ThetaDir          = 'clockwise';
pax4.RLim              = [0 0.25];  % shrink radial scale
pax4.RTick             = [];
pax4.GridAlpha         = 0.75;
pax4.FontSize          = 6;
pax4.ThetaTickLabel    = [];        % hide azimuth labels

%%%%%%% 2023 %%%%%%%%%%% 
%%%%%%%%%%

axG5=axes('color',[0.75 0.75 0.75],'position',[0.475 0.75-0.075 0.35 0.15-0.0125]);
[~,h]=contourf(XX,G.S43.yi,G.S43.Cdata,[0:20:450, -1],'LineStyle','none');

clim([0 400]);
% Set the colormap and add a transparent color for the NaNs
cmap = cmocean('tarn', 'pivotpoint', 1.4*43.7);
cmap = [[0.9 0.9 0.9]; cmap];  % Add a color for NaNs at the beginning
colormap(gca, cmap);

axis ij; ylim([0 400]);
hold on
contour(XX,G.S43.yi,G.S43.B,[1.4*43.7 1.4*43.7],'linewidth',0.75,...
    'linecolor',[0.1 0.1 0.1],'linestyle','-');
l=area(canyonXX,canyonZ*-1,1000,'facecolor',[0.8 0.8 0.8],'facealpha',1,...
    'edgecolor','none');
[CS,CH]=contour(XX,G.S43.yi,G.S43.Bs,[26 26.5 26.75],'linewidth',0.5,...
    'linecolor',[0.5 0.5 0.5]);
[CS2,CH2]=contour(XX,G.S43.yi,G.S43.Bs,[26.65 26.65],'linewidth',0.5,...
    'linecolor',[0.7 0.7 0.7]);
clabel(CS2,CH2,'LabelSpacing',200,'color',[0.6 0.6 0.6],'fontsize',5,...
    'fontname','Latin Modern Roman');
clabel(CS,CH,'LabelSpacing',200,'color',[0.4 0.4 0.4],'fontsize',6,...
    'fontname','Latin Modern Roman');
yticks(0:150:450); xticks(-30:30:120); xticklabels([]); set(gca, 'TickDir', 'out')
cl=clim;
text(10,350,G.S43.t);
yticklabels([]);

msk=datenum(w23.t)>G.S43.tt(1) & datenum(w23.t)<G.S43.tt(2);
dir  = w23.dir(msk);            % meteorological ° (0 =N, 90 =E, …)
spd  = w23.speed(msk);          % m s-¹  (optional, not used below)
frac   = 0.40;          % inset side = 30 % of parent-axes width/height
pad    = 0.01;          % 1 % figure padding from the edges of axG1
posG1  = get(axG5,'Position');            % [left bottom width height]
rosePos = [ posG1(1) + posG1(3) - posG1(3)*frac - pad , ... % left
            posG1(2)                           + pad , ... % bottom
            posG1(3)*frac , ...                           % width
            posG1(4)*frac ];                              % height
rosePos(1)=rosePos(1)+0.047;

pax5 = polaraxes('Position',rosePos,'Color','none');
hold(pax5,'on')
nbins = 16;
polarhistogram(pax5,deg2rad(dir),nbins, ...
               'Normalization','probability', ...
               'FaceColor',rgb_x('rose'),'EdgeColor','none');
pax5.ThetaZeroLocation = 'top';     % 0° = North
pax5.ThetaDir          = 'clockwise';
pax5.RLim              = [0 0.25];  % shrink radial scale
pax5.RTick             = [];
pax5.GridAlpha         = 0.75;
pax5.FontSize          = 6;
pax5.ThetaTickLabel    = [];        % hide azimuth labels

%%%%%%%%%%%%
axG6=axes('color',[0.75 0.75 0.75],'position',[0.475 0.6-0.075 0.35 0.15-0.0125]);
[~,h]=contourf(XX,G.S44.yi,G.S44.Cdata,[0:20:450, -1],'LineStyle','none');
clim([0 400]);
% Set the colormap and add a transparent color for the NaNs
cmap = cmocean('tarn', 'pivotpoint', 1.4*43.7);
cmap = [[0.9 0.9 0.9]; cmap];  % Add a color for NaNs at the beginning
colormap(gca, cmap);
axis ij; ylim([0 400]);
hold on
contour(XX,G.S44.yi,G.S44.B,[1.4*43.7 1.4*43.7],'linewidth',0.75,...
    'linecolor',[0.1 0.1 0.1],'linestyle','-');
l=area(canyonXX,canyonZ*-1,1000,'facecolor',[0.8 0.8 0.8],'facealpha',1,...
    'edgecolor','none');
[CS,CH]=contour(XX,G.S44.yi,G.S44.Bs,[26 26.5 26.75],'linewidth',0.5,...
    'linecolor',[0.5 0.5 0.5]);
[CS2,CH2]=contour(XX,G.S44.yi,G.S44.Bs,[26.65 26.65],'linewidth',0.5,...
    'linecolor',[0.7 0.7 0.7]);
clabel(CS2,CH2,'LabelSpacing',200,'color',[0.6 0.6 0.6],'fontsize',5,...
    'fontname','Latin Modern Roman');
clabel(CS,CH,'LabelSpacing',200,'color',[0.4 0.4 0.4],'fontsize',6,...
    'fontname','Latin Modern Roman');
yticks(0:150:450); xticks(-30:30:120); xticklabels([]); set(gca, 'TickDir', 'out')
cl=clim;
text(10,350,G.S44.t);
yticklabels([]);

msk=datenum(w23.t)>G.S44.tt(1) & datenum(w23.t)<G.S44.tt(2);
dir  = w23.dir(msk);            % meteorological ° (0 =N, 90 =E, …)
spd  = w23.speed(msk);          % m s-¹  (optional, not used below)
frac   = 0.40;          % inset side = 30 % of parent-axes width/height
pad    = 0.01;          % 1 % figure padding from the edges of axG1
posG1  = get(axG6,'Position');            % [left bottom width height]
rosePos = [ posG1(1) + posG1(3) - posG1(3)*frac - pad , ... % left
            posG1(2)                           + pad , ... % bottom
            posG1(3)*frac , ...                           % width
            posG1(4)*frac ];                              % height
rosePos(1)=rosePos(1)+0.047;

pax6 = polaraxes('Position',rosePos,'Color','none');
hold(pax6,'on')
nbins = 16;
polarhistogram(pax6,deg2rad(dir),nbins, ...
               'Normalization','probability', ...
               'FaceColor',rgb_x('rose'),'EdgeColor','none');
pax6.ThetaZeroLocation = 'top';     % 0° = North
pax6.ThetaDir          = 'clockwise';
pax6.RLim              = [0 0.25];  % shrink radial scale
pax6.RTick             = [];
pax6.GridAlpha         = 0.75;
pax6.FontSize          = 6;
pax6.ThetaTickLabel    = [];        % hide azimuth labels


%%%%%%%%%%%%
axG7=axes('color',[0.75 0.75 0.75],'position',[0.475 0.45-0.075 0.35 0.15-0.0125]);
axG7.Color = [0.6 0.6 0.6];
[~,h]=contourf(XX,G.S45.yi,G.S45.Cdata,[0:20:450, -1],'LineStyle','none');
clim([0 400]);
% Set the colormap and add a transparent color for the NaNs
cmap = cmocean('tarn', 'pivotpoint', 1.4*43.7);
cmap = [[0.9 0.9 0.9]; cmap];  % Add a color for NaNs at the beginning
colormap(gca, cmap);

axis ij; ylim([0 400]);
hold on
contour(XX,G.S45.yi,G.S45.B,[1.4*43.7 1.4*43.7],'linewidth',0.75,...
    'linecolor',[0.1 0.1 0.1],'linestyle','-');
l=area(canyonXX,canyonZ*-1,1000,'facecolor',[0.8 0.8 0.8],'facealpha',1,...
    'edgecolor','none');
[CS,CH]=contour(XX,G.S45.yi,G.S45.Bs,[26 26.5 26.75],'linewidth',0.5,...
    'linecolor',[0.5 0.5 0.5]);
[CS2,CH2]=contour(XX,G.S45.yi,G.S45.Bs,[26.65 26.65],'linewidth',0.5,...
    'linecolor',[0.7 0.7 0.7]);
clabel(CS2,CH2,'LabelSpacing',200,'color',[0.6 0.6 0.6],'fontsize',5,...
    'fontname','Latin Modern Roman');
clabel(CS,CH,'LabelSpacing',200,'color',[0.4 0.4 0.4],'fontsize',6,...
    'fontname','Latin Modern Roman');
yticks(0:150:450); xticks(-30:30:120); xticklabels([]); set(gca, 'TickDir', 'out')
cl=clim;
text(10,350,G.S45.t);
yticklabels([]);

msk=datenum(w23.t)>G.S45.tt(1) & datenum(w23.t)<G.S45.tt(2);
dir  = w23.dir(msk);            % meteorological ° (0 =N, 90 =E, …)
spd  = w23.speed(msk);          % m s-¹  (optional, not used below)
frac   = 0.40;          % inset side = 30 % of parent-axes width/height
pad    = 0.01;          % 1 % figure padding from the edges of axG1
posG1  = get(axG7,'Position');            % [left bottom width height]
rosePos = [ posG1(1) + posG1(3) - posG1(3)*frac - pad , ... % left
            posG1(2)                           + pad , ... % bottom
            posG1(3)*frac , ...                           % width
            posG1(4)*frac ];                              % height
rosePos(1)=rosePos(1)+0.047;

pax7 = polaraxes('Position',rosePos,'Color','none');
hold(pax7,'on')
nbins = 16;
polarhistogram(pax7,deg2rad(dir),nbins, ...
               'Normalization','probability', ...
               'FaceColor',rgb_x('rose'),'EdgeColor','none');
pax7.ThetaZeroLocation = 'top';     % 0° = North
pax7.ThetaDir          = 'clockwise';
pax7.RLim              = [0 0.25];  % shrink radial scale
pax7.RTick             = [];
pax7.GridAlpha         = 0.75;
pax7.FontSize          = 6;
pax7.ThetaTickLabel    = [];        % hide azimuth labels

%%%%%%%%%%%%
axG8=axes('color',[0.75 0.75 0.75],'position',[0.475 0.3-0.075 0.35 0.15-0.0125]);
[CS,CH]=contourf(XX,G.S48.yi,G.S48.Cdata,[0:20:450, -1],'LineStyle','none');
clim(cl);
% Set the colormap and add a transparent color for the NaNs
cmap = cmocean('tarn', 'pivotpoint', 1.4*43.7);
cmap = [[0.9 0.9 0.9]; cmap];  % Add a color for NaNs at the beginning
colormap(gca, cmap);
axis ij; ylim([0 400]);
hold on
contour(XX,G.S48.yi,G.S48.B,[1.4*43.7 1.4*43.7],'linewidth',0.75,...
    'linecolor',[0.1 0.1 0.1],'linestyle','-');
l=area(canyonXX,canyonZ*-1,1000,'facecolor',[0.8 0.8 0.8],'facealpha',1,...
    'edgecolor','none');
[CSc,CHc]=contour(XX,G.S48.yi,G.S48.Bs,[26 26.5 26.75],'linewidth',0.5,...
    'linecolor',[0.5 0.5 0.5]);
[CS2,CH2]=contour(XX,G.S48.yi,G.S48.Bs,[26.65 26.65],'linewidth',0.5,...
    'linecolor',[0.7 0.7 0.7]);
clabel(CS2,CH2,'LabelSpacing',200,'color',[0.6 0.6 0.6],'fontsize',5,...
    'fontname','Latin Modern Roman');
clabel(CSc,CHc,'LabelSpacing',200,'color',[0.4 0.4 0.4],'fontsize',6,...
    'fontname','Latin Modern Roman');
yticks(0:150:450); xticks(-30:30:120); xticklabels([]); set(gca, 'TickDir', 'out')
text(10,350,G.S48.t);
yticklabels([]);

msk=datenum(w23.t)>G.S48.tt(1) & datenum(w23.t)<G.S48.tt(2);
dir  = w23.dir(msk);            % meteorological ° (0 =N, 90 =E, …)
spd  = w23.speed(msk);          % m s-¹  (optional, not used below)
frac   = 0.40;          % inset side = 30 % of parent-axes width/height
pad    = 0.01;          % 1 % figure padding from the edges of axG1
posG1  = get(axG8,'Position');            % [left bottom width height]
rosePos = [ posG1(1) + posG1(3) - posG1(3)*frac - pad , ... % left
            posG1(2)                           + pad , ... % bottom
            posG1(3)*frac , ...                           % width
            posG1(4)*frac ];                              % height
rosePos(1)=rosePos(1)+0.047;

pax8 = polaraxes('Position',rosePos,'Color','none');
hold(pax8,'on')
nbins = 16;
polarhistogram(pax8,deg2rad(dir),nbins, ...
               'Normalization','probability', ...
               'FaceColor',rgb_x('rose'),'EdgeColor','none');
pax8.ThetaZeroLocation = 'top';     % 0° = North
pax8.ThetaDir          = 'clockwise';
pax8.RLim              = [0 0.25];  % shrink radial scale
pax8.RTick             = [];
pax8.GridAlpha         = 0.75;
pax8.FontSize          = 6;
pax8.ThetaTickLabel    = [];        % hide azimuth labels

%%%%%%%%%%%%
axG9=axes('color',[0.75 0.75 0.75],'position',[0.475 0.15-0.075 0.35 0.15-0.0125]);
contourf(XX,G.S49.yi,G.S49.Cdata,[0:20:450, -1],'LineStyle','none');
clim(cl);
cmap = cmocean('tarn', 'pivotpoint', 1.4*43.7);
cmap = [[0.9 0.9 0.9]; cmap];  % Add a color for NaNs at the beginning
colormap(gca, cmap);
axis ij; ylim([0 400]);
hold on
contour(XX,G.S49.yi,G.S49.B,[1.4*43.7 1.4*43.7],'linewidth',0.75,...
    'linecolor',[0.1 0.1 0.1],'linestyle','-');
l=area(canyonXX,canyonZ*-1,1000,'facecolor',[0.8 0.8 0.8],'facealpha',1,...
    'edgecolor','none');
[CS,CH]=contour(XX,G.S49.yi,G.S49.Bs,[26 26.5 26.75],'linewidth',0.5,...
    'linecolor',[0.5 0.5 0.5]);
[CS2,CH2]=contour(XX,G.S49.yi,G.S49.Bs,[26.65 26.65],'linewidth',0.5,...
    'linecolor',[0.7 0.7 0.7]);
clabel(CS2,CH2,'LabelSpacing',200,'color',[0.6 0.6 0.6],'fontsize',5,...
    'fontname','Latin Modern Roman');
clabel(CS,CH,'LabelSpacing',200,'color',[0.4 0.4 0.4],'fontsize',6,...
    'fontname','Latin Modern Roman');
yticks(0:150:450); xticks(-30:30:120); set(gca, 'TickDir', 'out')
xlabel('Distance from 300 m isobath (km)');
text(10,350,G.S49.t);
yticklabels([]);

msk=datenum(w23.t)>G.S49.tt(1) & datenum(w23.t)<G.S49.tt(2);
dir  = w23.dir(msk);            % meteorological ° (0 =N, 90 =E, …)
spd  = w23.speed(msk);          % m s-¹  (optional, not used below)
frac   = 0.40;          % inset side = 30 % of parent-axes width/height
pad    = 0.01;          % 1 % figure padding from the edges of axG1
posG1  = get(axG9,'Position');            % [left bottom width height]
rosePos = [ posG1(1) + posG1(3) - posG1(3)*frac - pad , ... % left
            posG1(2)                           + pad , ... % bottom
            posG1(3)*frac , ...                           % width
            posG1(4)*frac ];                              % height
rosePos(1)=rosePos(1)+0.047;

pax9 = polaraxes('Position',rosePos,'Color','none');
hold(pax9,'on')
nbins = 16;
polarhistogram(pax9,deg2rad(dir),nbins, ...
               'Normalization','probability', ...
               'FaceColor',rgb_x('rose'),'EdgeColor','none');
pax9.ThetaZeroLocation = 'top';     % 0° = North
pax9.ThetaDir          = 'clockwise';
pax9.RLim              = [0 0.25];  % shrink radial scale
pax9.RTick             = [];
pax9.GridAlpha         = 0.75;
pax9.FontSize          = 6;
pax9.ThetaTickLabel    = [];        % hide azimuth labels

%%%%%%%%%%%%
set(findall(gcf,'-property','fontsize'),'fontsize',8);

% Define labels
labels = {'a)', 'b)', 'c)', 'd)', 'e)', 'f)', 'g)', 'h)', 'i)', 'j)', 'k)', 'l)', 'm)', 'n)', 'o)'};

% Define an array to store the axes handles
axHandles = {axMap axG1 axG2 axG3 axG4 axG5 axG6 axG7 axG8 axG9};

% Position text labels on each subplot
for idx = 2:length(axHandles)
    axes(axHandles{idx}); % Switch to the current axes

    if idx==1

        % Add label to the bottom right corner
        text(0.985, 0.825, labels{idx}, 'Units', 'normalized', ...
            'HorizontalAlignment', 'right', 'VerticalAlignment', 'bottom', ...
            'FontSize', 7,'fontweight','bold');
    else
        % Add label to the bottom right corner
        text(0.075, 0.8, labels{idx}, 'Units', 'normalized', ...
            'HorizontalAlignment', 'right', 'VerticalAlignment', 'bottom', ...
            'FontSize', 7,'fontweight','bold');
    end
end
axes(CB);

uistack([pax1 pax2 pax3 pax4 pax5 pax6 pax7 pax8 pax9],'top');

set(findall(gcf, '-property', 'fontsize'), 'fontsize', 8);
set(findall(gcf,'-property','fontname'),'fontname','Latin Modern Roman');


%%
export_fig /Users/samuelstevens/Dropbox/Hakai/gliders/oxygen/figures/glider_only_R1.pdf -dpdf -nofontswap
% %
% 
% exportgraphics(gcf, '/Users/samuelstevens/Dropbox/Hakai/gliders/oxygen/figures/glider_only.pdf', ...
%                'ContentType', 'vector', 'Resolution', 300);

% export_fig /Users/samuelstevens/Dropbox/Hakai/gliders/oxygen/figures/poster/glider.png -m6 -nofontswap -transparent


%% Historical context maps
col=[255 214 140]/255; % YELLOW!
lat_lim=[50.5 52.5]; lon_lim=[-131 -127];
fname='/Users/samst/Dropbox/UBC/Misc/british_columbia_3_msl_2013.nc';
lat=ncread(fname,'lat');    
lon=ncread(fname,'lon');
ilon=lon>=lon_lim(1) & lon<=lon_lim(2);
ilat=lat>=lat_lim(1) & lat<=lat_lim(2);
Z=ncread(fname,...
    'Band1',[ find((ilon),1,'first') find((ilat),1,'first')],...
    [ sum(ilon) sum(ilat)],[1 1]);

figure('units','centimeters','outerposition',[0 0 18 17],'color','w');

ax1=axes('position',[0.075 0.55 0.375 0.375]);
m_proj('equidistant','lon',lon_lim,'lat',lat_lim);   % Projection
hold on
deep=-4000:1000:-1000; shelf=[-600:50:-100 0];
[CS,CH]=m_contourf(lon(ilon),lat(ilat),Z',[deep shelf],'linestyle','none');

g=cmocean('gray',600);g(1:390,:)=[];g(end-20:end,:)=[];
cm=[g;m_colmap('blues',50)];
colormap(ax1,cm);

ax = colorbar('Position', [0.5 0.65 0.01 0.1]); 
set(ax, 'ylim', [-375 0]);
ax.Ticks=[-300  -200  -100     0];
ax.TickLabels={'300';'200';'100';'0'};
title(ax,'Depth (m)','Interpreter','latex')
CH.FaceAlpha=0.9;

[CS,CH]=m_etopo2('contour',[-200 -200],'edgecolor',rgb_x('grey'));
clabel(CS,CH,'labelspacing',300,'color',rgb_x('grey'),'fontsize',6,...
    'fontname','Latin Modern Roman');
m_gshhs_f('patch',col,'edgecolor',[0.3 0.3 0.3],'linewidth',0.1);
m_grid('tickdir','out','linestyle','none','xtick',5,...
    'YaxisLocation','left');
m_annotation('textbox',[-130.4,50.7,0.01,0.01],'string','2003-2021','fontsize',6,...
    'FitBoxToText','on','horizontalalignment','center','color',rgb_x('sky blue'),...
    'backgroundcolor','w');

text(0.025, 0.95, 'a)', 'Units', 'normalized');

mnth=month(bottomDO.time);
% Add bottom DO data 
Cmsk1=find(bottomDO.time>datenum(2003,01,01) & bottomDO.time<datenum(2022,01,01) &...
    bottomDO.z<=350 & bottomDO.latitude<52.5 & mnth>=5 & mnth<=10);
Cmsk2=find(bottomDO.time>datenum(2022,01,01) & bottomDO.time<datenum(2024,01,01) & bottomDO.z<=350 &...
    bottomDO.latitude<52.5  & mnth>=5 & mnth<=10);

count1=0; count2=0;
for i=1:length(Cmsk1)
    if bottomDO.oxygen(Cmsk1(i))<1.4*43.7
        count1=count1+1;
        s1(count1)=m_scatter(bottomDO.longitude(Cmsk1(i)),bottomDO.latitude(Cmsk1(i))....
            ,50,'rx','linewidth',1.5);

    else
        count2=count2+1;
        s2(count2)=m_scatter(bottomDO.longitude(Cmsk1(i)),bottomDO.latitude(Cmsk1(i)),10,'x',...
            'markeredgecolor',[0.4 0.4 0.4]);
   
    end
end
uistack(s1,'top');


ax4=axes('position',[0.55 0.55 0.375 0.375]);
m_proj('equidistant','lon',lon_lim,'lat',lat_lim);   % Projection
hold on
deep=-4000:1000:-1000; shelf=[-600:50:-100 0];
[CS,CH]=m_contourf(lon(ilon),lat(ilat),Z',[deep shelf],'linestyle','none');
colormap(cm);
CH.FaceAlpha=0.9;
[CS,CH]=m_etopo2('contour',[-200 -200],'edgecolor',rgb_x('grey'));
clabel(CS,CH,'labelspacing',300,'color',rgb_x('grey'),'fontsize',6,...
    'fontname','Latin Modern Roman');

m_gshhs_f('patch',col,'edgecolor',[0.3 0.3 0.3],'linewidth',0.1);
m_grid('tickdir','out','linestyle','none','xtick',5,...
    'Yticklabels',[]);
m_annotation('textbox',[-130.4,50.7,0.01,0.01],'string','2022-2023','fontsize',6,...
    'FitBoxToText','on','horizontalalignment','center','color',rgb_x('deep pink'),...
    'backgroundcolor','w');
text(0.025, 0.95,'b)','fontsize',8,'units','normalized');

count3=0; count4=0;
for i=1:length(Cmsk2)
    if bottomDO.oxygen(Cmsk2(i))<1.4*43.7
        count3=count3+1;
        s3(count3)=m_scatter(bottomDO.longitude(Cmsk2(i)),bottomDO.latitude(Cmsk2(i)),...
            50,'rx','linewidth',1.5);

    else
        count4=count4+1;
        s4(count4)=m_scatter(bottomDO.longitude(Cmsk2(i)),bottomDO.latitude(Cmsk2(i)),10,'x',...
            'markeredgecolor',[0.4 0.4 0.4]);
    end
end
uistack(s3,'top');


ax3=axes('position',[0.075,0.1,0.36,0.375]); hold on;
yr=year(bottomDO.time);
yrs=2003:2023;
mnth=month(bottomDO.time);
mnths=NaN(length(yrs),12);
cmb=m_colmap('mBOD',size(mnths,2)); 

for i=[1:4 11:12]
    cmb(i,:)=[0.7 0.7 0.7];
end

for i=1:length(yrs)
    msk=yr==yrs(i);
    mnths(i,:)=histcounts(mnth(msk),0.5:1:12.5);
end
bh=bar(yrs,mnths,'stacked','facecolor','flat','edgecolor','none');
for i=1:size(mnths,2)
    bh(i).CData = cmb(i,:);
end
ylabel('\it{n}');
xticks(2005:5:2024)
xlabel('Year');
grid on
set(gca,'tickdir','out');
text(0.025, 0.95,'c)','fontsize',8,'units','normalized');


ax2=axes('position',[0.55,0.1,0.36,0.375]);
h1=histogram(bottomDO.oxygen(Cmsk1),'binedges',0:10:350,'facecolor',rgb_x('sky blue'),...
    'edgecolor',rgb_x('grey'));
hold on
h2=histogram(bottomDO.oxygen(Cmsk2),0:10:350,'facecolor',rgb_x('salmon'),...
    'edgecolor',rgb_x('black'));
h1.FaceAlpha=0.4;
h2.FaceAlpha=0.95;
xlim([0 350]); grid on; xticks([0 61 100:50:350]); yticks([0:50:200]);

uistack(h1,'top'); yl=ylim;
line([mean(bottomDO.oxygen(Cmsk1),'omitnan') mean(bottomDO.oxygen(Cmsk1),'omitnan')],...
    yl,'color',rgb_x('azure'),'linewidth',1,'linestyle','-');
line([mean(bottomDO.oxygen(Cmsk2),'omitnan') mean(bottomDO.oxygen(Cmsk2),'omitnan')],...
    yl,'color',rgb_x('deep pink'),'linewidth',1,'linestyle','--');
xlabel('$O_\mathrm{b}$ ($\mu$mol kg$^{-1}$)','Interpreter','latex')
ylabel('\it{n}');
box off
set(gca,'tickdir','out');
text(0.025, 0.95, 'd)', 'Units', 'normalized');

%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%
axes(ax3)
for i=1:size(cmb,1)
    scatter(2003.5,130-i*4,15,'s','markerfacecolor',cmb(i,:),'markeredgecolor','none');
    if i==1
        text(2004,130-i*4,'Jan.','fontsize',8);
    elseif i==6
        text(2004,130-i*4,'Jun.','fontsize',8);
    elseif i==12
        text(2004,130-i*4,'Dec.','fontsize',8);
    end
end
% scatter(2003.5,130-9*3,10,'s','markerfacecolor',cmb(9,:),'markeredgecolor','none');

set(findall(gcf,'-property','fontsize'),'fontsize',11);
set(findall(gcf,'-property','fontname'),'fontname','Latin Modern Roman');
% axes(ax2);
% ax2.XTickLabelRotation=45;

%%
export_fig /Users/samuelstevens/Dropbox/Hakai/gliders/oxygen/figures/historical_R1.pdf -dpdf -nofontswap
% exportgraphics(gcf, '/Users/samuelstevens/Dropbox/Hakai/gliders/oxygen/figures/historicalHR.pdf', ...
%     'ContentType', 'vector', 'Resolution', 300);

% export_fig /Users/samuelstevens/Dropbox/Hakai/gliders/oxygen/figures/poster/history.png -m6 -nofontswap -transparent

%% Compare pre and post 2022 periods
clc

% Prepare data for the first period (2003-2020)
data1 = bottomDO.oxygen(Cmsk1);
data1 = data1(~isnan(data1)); % Remove NaNs
n1 = length(data1);

if n1 > 1
    nBoot = 1000; % Number of bootstrap samples
    boot_means1 = bootstrp(nBoot, @mean, data1);
    ci1 = prctile(boot_means1, [2.5 97.5]);
    uncertainty1 = (ci1(2) - ci1(1)) / 2;
else
    ci1 = [NaN NaN];
    uncertainty1 = NaN;
end

mean1 = mean(data1);
% uncertainty1=nanstd(data1)*2;

% Calculate the percentage of hypoxic observations for the first period
percent1 = (count1 / length(Cmsk1)) * 100;

% Prepare data for the second period (after 2021)
data2 = bottomDO.oxygen(Cmsk2);
data2 = data2(~isnan(data2)); % Remove NaNs
n2 = length(data2);

if n2 > 1
    nBoot = 1000; % Number of bootstrap samples
    boot_means2 = bootstrp(nBoot, @mean, data2);
    ci2 = prctile(boot_means2, [2.5 97.5]);
    uncertainty2 = (ci2(2) - ci2(1)) / 2;
else
    ci2 = [NaN NaN];
    uncertainty2 = NaN;
end

mean2 = mean(data2);
% uncertainty2=nanstd(data2)*2;

% Calculate the percentage of hypoxic observations for the second period
percent2 = (count3 / length(Cmsk2)) * 100;

% Print results
fprintf('%3.1f%% hypoxic O_b 2003-2020 (n=%3.0f, mean=%4.1f ± %4.1f µmol/kg)\n', ...
    percent1, n1, mean1, uncertainty1);
fprintf('%3.1f%% hypoxic O_b after 2021 (n=%3.0f, mean=%4.1f ± %4.1f µmol/kg)\n', ...
    percent2, n2, mean2, uncertainty2);

%% check other periods 
wnd=2; p=[];
for i=2003:1:2024-wnd
    Cmsk3=find(bottomDO.time>datenum(i,01,01) & bottomDO.time<datenum(i+wnd,01,01) &...
        bottomDO.z<=350 & bottomDO.latitude<52.5 & mnth>=5 & mnth<=10 &...
        (bottomDO.latitude<52 & bottomDO.longitude>-130.263));
    Cmsk4=find(bottomDO.oxygen(Cmsk3)<1.4*43.7);

    p(i-2002)=(length(Cmsk4)/length(Cmsk3))*100;
end


%% Check single month (September)

% check other periods 
wnd=1; p=[];
for i=2003:1:2024-wnd
    Cmsk3=find(bottomDO.time>=datenum(i,08,01) & bottomDO.time<datenum(i+wnd,10,01) &...
        bottomDO.z<=350 & bottomDO.latitude<52.5 & mnth>=9 & mnth<=10 &...
        (bottomDO.latitude<52 & bottomDO.longitude>-130.263));
    Cmsk4=find(bottomDO.oxygen(Cmsk3)<1.4*43.7);

    p(i-2002)=(length(Cmsk4)/length(Cmsk3))*100;
end


Smsk1=find(bottomDO.time>=datenum(2003,01,01) & bottomDO.time<datenum(2022,01,01) &...
        bottomDO.z<=350 & bottomDO.latitude<52.5 & mnth>=9 & mnth<=10 &...
        (bottomDO.latitude<52 & bottomDO.longitude>-130.263));

Smsk2=find(bottomDO.time>=datenum(2022,01,01) & bottomDO.time<datenum(2024,01,01) &...
        bottomDO.z<=350 & bottomDO.latitude<52.5 & mnth>=9 & mnth<=10 &...
        (bottomDO.latitude<52 & bottomDO.longitude>-130.263));
clc
fprintf('2003-2021 (n=%3.0f, mean=%4.1f+-%4.1f)\n 2022-23 (n=%3.0f, mean=%4.1f+-%4.1f)\n',...
    length(Smsk1),mean(bottomDO.oxygen(Smsk1),'omitnan'),(2*std(bottomDO.oxygen(Smsk1),'omitnan'))/sqrt(length(Smsk1)),...
    length(Smsk2),mean(bottomDO.oxygen(Smsk2),'omitnan'),(2*std(bottomDO.oxygen(Smsk2),'omitnan'))/sqrt(length(Smsk2)));


%% Check single month at CS09 (September)

load CS09_ctd.mat

mnth=month(CS09.mtime);
ob=min(CS09.ox(174:end,:));

Cmsk1=find(CS09.mtime>datenum(2003,1,1) & CS09.mtime<datenum(2022,1,1) & mnth==9);
Cmsk2=find(CS09.mtime>datenum(2022,1,1) & CS09.mtime<datenum(2024,1,1) & mnth==9);

count1=sum(ob(Cmsk1)<61);
count3=sum(ob(Cmsk2)<61);

clc

% Prepare data for the first period (2003-2020)
data1 = ob(Cmsk1);
data1 = data1(~isnan(data1)); % Remove NaNs
n1 = length(data1);

if n1 > 1
    nBoot = 1000; % Number of bootstrap samples
    boot_means1 = bootstrp(nBoot, @mean, data1);
    ci1 = prctile(boot_means1, [2.5 97.5]);
    uncertainty1 = (ci1(2) - ci1(1)) / 2;
else
    ci1 = [NaN NaN];
    uncertainty1 = NaN;
end

mean1 = mean(data1);

% Calculate the percentage of hypoxic observations for the first period
percent1 = (count1 / length(Cmsk1)) * 100;

% Prepare data for the second period (after 2021)
data2 = ob(Cmsk2);
data2 = data2(~isnan(data2)); % Remove NaNs
n2 = length(data2);

if n2 > 1
    nBoot = 1000; % Number of bootstrap samples
    boot_means2 = bootstrp(nBoot, @mean, data2);
    ci2 = prctile(boot_means2, [2.5 97.5]);
    uncertainty2 = (ci2(2) - ci2(1)) / 2;
else
    ci2 = [NaN NaN];
    uncertainty2 = NaN;
end

mean2 = mean(data2);

% Calculate the percentage of hypoxic observations for the second period
percent2 = (count3 / length(Cmsk2)) * 100;

% Print results
fprintf('%3.1f%% hypoxic O_b 2003-2020 (n=%3.0f, mean=%4.1f ± %4.1f µmol/kg)\n', ...
    percent1, n1, mean1, uncertainty1);
fprintf('%3.1f%% hypoxic O_b after 2021 (n=%3.0f, mean=%4.1f ± %4.1f µmol/kg)\n', ...
    percent2, n2, mean2, uncertainty2);


%% How much lower is 23 than other years
yr=year(M.Hak1.time);
yrs=unique(yr);

for i=1:length(yrs)
    msk=yr==yrs(i);
    minO(i)=min(M.Hak1.CTD133m.oxygen(msk));
end

%% CUI on July 13th 
dday=day(datetime(2023,07,19),'dayofyear');

[~,rnk]=sort(CUI(2:end,dday),'descend');
idx=find(rnk==length(1968:2023))


%% Look at along isopycnal AOU gradient

allProfs.AOU=aou(allProfs.salinity,allProfs.potential_temperature,...
    allProfs.oxygen_corrected);
allProfs.isoAOU=NaN(3,length(allProfs.AOU));
allProfs.isoDO=NaN(3,length(allProfs.AOU));
allProfs.isoSpice=NaN(3,length(allProfs.AOU));

sig_thresh=[26 26.5 26.7];

for i=1:size(allProfs.AOU,2)
        for ii=1:length(sig_thresh)
            msk=~isnan(allProfs.AOU(:,i));
            if sum(msk)>1
                allProfs.isoAOU(ii,i)=interp1(inpaint_nans(allProfs.s_t(msk,i)),allProfs.AOU(msk,i),...
                    sig_thresh(ii));
                allProfs.isoDO(ii,i)=interp1(inpaint_nans(allProfs.s_t(msk,i)),allProfs.oxygen_corrected(msk,i),...
                    sig_thresh(ii));
                allProfs.isoSpice(ii,i)=interp1(inpaint_nans(allProfs.s_t(msk,i)),allProfs.spice(msk,i),...
                    sig_thresh(ii));
            end
        end
end


%% make along canyon average

xGrid=0.25:0.25:159.75;
isoAOUmn=NaN(3,length(xGrid));
isoAOUsem=NaN(3,length(xGrid));
ctdAOUmn=NaN(3,length(xGrid));
ctdAOUsem=NaN(3,length(xGrid));

load isoDO_DFOandHakai.mat
isoDO.AOU(isoDO.AOU>275)=NaN;

count=0;
for i=0:0.25:159.5
    count=count+1;
    msk=allProfs.canyonX>=i & allProfs.canyonX<=i+0.25;
    
    isoAOUmn(:,count)=mean(allProfs.isoAOU(:,msk),2,'omitnan');
    isoAOUsem(:,count)=2.*std(allProfs.isoAOU(:,msk),0,2,'omitnan');%./...
        % sqrt(sum(~isnan(allProfs.isoAOU(:,msk)),2,'omitnan'));

    msk=isoDO.OSdist/1e3>=i & isoDO.OSdist/1e3<=i+0.25;
    ctdAOUmn(:,count)=mean(isoDO.AOU(:,msk),2,'omitnan');
    ctdAOUsem(:,count)=2.*std(isoDO.AOU(:,msk),0,2,'omitnan');%./...
        % sqrt(sum(~isnan(isoDO.AOU(:,msk)),2,'omitnan'));

end
%%
xGridS=(xGrid-120).*-1;

figure('units','centimeters','outerposition',[0 0 12 10],'color','w');
t = tiledlayout(2,2,'TileSpacing','Compact','Padding','Compact');
nexttile(t,1:2);
load('BCcoast');
lat_lim=[50.75 52]; lon_lim=[-130.5 -127.5];

fname='/Users/samst/Dropbox/UBC/Misc/british_columbia_3_msl_2013.nc';
lat=ncread(fname,'lat');    
lon=ncread(fname,'lon');
ilon=lon>=lon_lim(1) & lon<=lon_lim(2);
ilat=lat>=lat_lim(1) & lat<=lat_lim(2);
Z=ncread(fname,...
    'Band1',[ find((ilon),1,'first') find((ilat),1,'first')],...
    [ sum(ilon) sum(ilat)],[1 1]);

col=[255 214 140]/255; % YELLOW!
m_proj('equidistant','lon',lon_lim,'lat',lat_lim);   % Projection
hold on
deep=-4000:1000:-1000; shelf=[-600:50:-100 0];
[CS,CH]=m_contourf(lon(ilon),lat(ilat),Z',[deep shelf],'linestyle','none');

g=cmocean('gray',600);g(1:390,:)=[];g(end-20:end,:)=[];
cm=[g;m_colmap('blues',50)];
colormap(cm);

ax = colorbar('Position', [0.9 0.8 0.01 0.1]); 
set(ax, 'ylim', [-375 0]);
ax.Ticks=[-300  -200  -100     0];
ax.TickLabels={'300';'200';'100';'0'};
title(ax,'Depth (m)','Interpreter','latex')
CH.FaceAlpha=0.9;

m_gshhs_f('patch',col,'edgecolor',[0.3 0.3 0.3],'linewidth',0.3);
m_grid('linestyle','none','linewidth',0.5,'tickdir','out','xtick',3);

m_plot(canyonAX(:,1),canyonAX(:,2),'k','linewidth',3);

% Build axis
lnes=linspecer(length(0:25:max(canyonX)));
count=0;

for i=0:30:max(canyonX)
    [~,idx]=min(abs(canyonX-i));
    count=count+1;
    m_scatter(canyonAX(idx,1),canyonAX(idx,2),30,'filled','markerfacecolor','k',...
        'markeredgecolor','w','marker','s');
        m_text(canyonAX(idx,1)-0.4,canyonAX(idx,2),sprintf('%3.0f km',canyonX(idx)-120));

end
text(0.05, 0.9, 'a)', 'Units', 'normalized');

nexttile(1:2);
hold on;
lnes = linspecer(5);lnes(3:4,:)=[]; lnes(1,:)=[0.5 0.5 0.5]; lnes=flipud(lnes);
coeffs=[]; 
coeffs23=[];
for i=1:length(sig_thresh)
    fittedy=[];
    msk=~isnan(isoAOUmn(i,:));
    % errorbar(xGridS,ctdAOUmn(i,:), ctdAOUsem(i,:), 'LineWidth', 0.25, ...
    %     'CapSize', 2, 'Color', [lnes(i,:)-[0.2 0.2 0.2] 0.025], 'LineStyle', 'none');
    fill([xGridS(msk) fliplr(xGridS(msk))]',...
        [isoAOUmn(i,msk)-isoAOUsem(i,msk) fliplr(isoAOUmn(i,msk)+isoAOUsem(i,msk))]',...
        lnes(i,:),'edgecolor','none','facealpha',0.1);

    scatter(xGridS,isoAOUmn(i,:),5,'filled','markerfacecolor',lnes(i,:),...
        'markeredgecolor','none','markerfacealpha',0.3);
    
    scatter(xGridS,ctdAOUmn(i,:),5,'filled','markerfacecolor',lnes(i,:)-[0.1 0.1 0.1],...
        'markeredgecolor','none','marker','^','markerfacealpha',0.7);

    trendMsk=~isnan(isoAOUmn(i,:)) & xGrid<=120;
    [coeffs(i,1),coeffs(i,2)]=tsreg(xGridS(trendMsk), isoAOUmn(i,trendMsk));
    % [p,yhat,CS09.ci1]=polypredci(yrGrid(~isnan(CS09.Obmn)),...
    %     CS09.Obmn(~isnan(CS09.Obmn)),1);
    fittedy(i,:)=polyval(coeffs(i,:),xGridS(trendMsk));

    plot(xGridS(trendMsk),fittedy(i,:),'color','w','LineWidth',3);
    sp(i)=plot(xGridS(trendMsk),fittedy(i,:),'color',lnes(i,:),'LineWidth',1.5);

    fittedy=[];
    trendMsk=~isnan(isoAOUmn(i,:)) & xGrid>120;
    [coeffs(i,1),coeffs(i,2)]=tsreg(xGridS(trendMsk), isoAOUmn(i,trendMsk));
    fittedy(i,:)=polyval(coeffs(i,:),xGridS(trendMsk));

    % plot(xGridS(trendMsk),fittedy(i,:),'color','w','LineWidth',2);
    % sp(i)=plot(xGridS(trendMsk),fittedy(i,:),'--','color',lnes(i,:));
    % 
    % msk=~isnan(isoAOU23mn(i,:));
    % [coeffs23(i,1),coeffs23(i,2)]=tsreg(xGrid(msk), isoAOU23mn(i,msk));
    % fittedy(i,:)=polyval(coeffs23(i,:),xGrid);
    % 
    % plot(xGrid,fittedy(i,:),'color','w','LineWidth',2);
    % sp(i)=plot(xGrid,fittedy(i,:),'--','color',lnes(i,:));
    
end

%set(gca,'xdir','reverse');
xlabel('Distance from 300 m isobath (km)');
ylabel('AOU ($\mu$mol kg$^{-1}$)','Interpreter','latex');
grid on; axis tight; ylim([80 270]);
text(0.05, 0.9, 'b)', 'Units', 'normalized');

set(findall(gcf,'-property','fontsize'),'fontsize',10);
% legend(sp,'26.00 kg m^{-3}','26.50 kg m^{-3}','26.75 kg m^{-3}','fontsize',...
%     6,'location','best');

set(findall(gcf,'-property','fontname'),'fontname','Latin Modern Roman');

% Rough transit time estimate
% 
load SWSadcpS2.mat

msk1=month(adcp{1, 1}.Time)>=5 & month(adcp{1, 1}.Time)<=10;
msk2=adcp{1, 1}.bin_depth>100;
meanAcrossU=mean(adcp{1,1}.velcs(msk1,msk2),'omitnan');
% TT(1)=(120e3/max(meanAcrossU))/86400; %days
TT=(120e3/mean(meanAcrossU,'omitnan'))/86400; %days
% TT(2)=(120e3/min(meanAcrossU))/86400; %days

OUR=(coeffs(1,1).*120)./TT
OUR=(coeffs(3,1).*120)./TT


%%
% exportgraphics(gcf, '/Users/samuelstevens/Dropbox/Hakai/gliders/oxygen/figures/AOUvsXv5.pdf', ...
%     'ContentType', 'vector', 'Resolution', 100);
% export_fig /Users/samuelstevens/Dropbox/Hakai/gliders/oxygen/figures/AOUvsXv5.pdf -dpdf -nofontswap
%% Cumulative across shelf-break flow (estuarine)

msk2=adcp{1, 1}.bin_depth>200;
meanAcrossU=mean(adcp{1,1}.velcs(:,msk2),2,'omitnan');
 
% figure;
% plot(adcp{1,1}.Time,smoothdata(meanAcrossU,'movmean',24*7*2));

EI=NaN(7,365);
dday=day(adcp{1,1}.Time,'dayofyear');
yr=year(adcp{1,1}.Time);

for i=2017:2023
    for ii=1:365
        msk=dday==ii & yr==i;
        EI(i-2016,ii)=mean(meanAcrossU(msk),'omitnan');
    end
end

CEI=cumsum(EI,2,'omitnan');
l=[];
cm=linspecer(7);
figure; hold on;
for i=1:7
    l(i)=plot(1:365,CEI(i,:),'color',cm(i,:),'linewidth',2);
end
legend(l,'2017','2018','2019','2020','2021','2022','2023');

%% Across shelf currents

figure('color','w'); 
subplot(2,1,1);
hold on;
plot(adcp{1, 1}.Time,meanAcrossU,'color',[0.9 0.9 0.9]);
plot(adcp{1, 1}.Time,smoothdata(meanAcrossU,'movmean',24*14),'r');
axis tight; grid on;
ylim([-0.2 0.2]);
% ylabel('Onshore current > 200 m (cm s^{-1})')

subplot(2,1,2); hold on;
lnes=linspecer(6);
for i=1:6
    plot(1:365,EI(i,:),'color',lnes(i,:))
end
plot(1:365,smoothdata(mean(EI,1,'omitnan'),'movmean',14),'k','linewidth',2);
axis tight; grid on;
ylabel('Onshore current > 200 m (cm s^{-1})')
xlabel('Julian Day');

set(findall(gcf, '-property', 'fontsize'), 'fontsize', 10);
set(findall(gcf, '-property', 'fontname'), 'fontname', 'Latin Modern Roman');

%%
% export_fig /Users/samuelstevens/Dropbox/Hakai/gliders/oxygen/figures/AOUvsXv3.pdf -dpdf -nofontswap

%% Does CS09 mean oxygen link with winter upwelling?

yrs=year(CS09.mtime);
yr=unique(yrs);
mnth=month(CS09.mtime);
mnCS=NaN(length(yr),1);
Byr=1967:2023;
cmb=NaN(length(yr),2);

for i=1:length(yr)
    msk=yrs==yr(i) & mnth==9;
    if sum(msk)>0
        mnCS(i)=mean(CS09.ox(:,msk),'omitnan');
    end
    cmb(i,:)=[i find(Byr==yr(i))];
end


figure;
scatter(CUIw(cmb(:,2)),mnCS)
hold on
scatter(CUIw(cmb(end,2)),mnCS(end),10,'filled','r')

cmb=NaN(length(yr),2);
mnCS=NaN(length(yr),1);
for i=1:length(yr)
    msk=yrs==yr(i) & mnth==7;
    if sum(msk)>0
        mnCS(i)=mean(mean(CS09.ox(:,msk),'omitnan'));
    end
    cmb(i,:)=[i find(Byr==yr(i))];
end

scatter(CUIw(cmb(:,2)),mnCS)
hold on
scatter(CUIw(cmb(end,2)),mnCS(end),10,'filled','r')

%%
KC=HakaiWaterPropertiesInstrumentP;
KC.mtime=datenum(1970,01,01,00,00,00)+(KC.time/86400);
KC.pressure=double(KC.pressure);
KC.dissolved_oxygen_ml_l=KC.dissolved_oxygen_ml_l.*43.7;
sa=gsw_SA_from_SP(KC.salinity,KC.pressure,KC.longitude,KC.latitude);
pt=gsw_pt0_from_t(sa,KC.temperature,KC.pressure);
ct=gsw_CT_from_t(sa,KC.temperature,KC.pressure);
KC.s_ti=gsw_sigma0(sa,ct);
KC.AOUi=aou(KC.salinity,pt,KC.dissolved_oxygen_ml_l);

t=unique(KC.mtime);
KC.ox=NaN(371,length(t));
KC.AOU=NaN(371,length(t));

KC.depth=KC.depth+rand(size(KC.depth))/1e5;
for i=1:length(t)
    msk=KC.mtime==t(i);
    KC.ox(:,i)=interp1(KC.depth(msk),KC.dissolved_oxygen_ml_l(msk),0:370)';
    KC.AOU(:,i)=interp1(KC.depth(msk),KC.AOUi(msk),0:370)';
    KC.s_t(:,i)=interp1(KC.depth(msk),KC.s_ti(msk),0:370)';

end
t(KC.ox(1,:)>550)=[];
KC.ox(:,KC.ox(10,:)>550)=[];
KC.AOU(:,KC.AOU(10,:)<-550)=[];

yrs=year(t);
yr=unique(yrs);
mnth=month(t);
mnKC=NaN(length(yr),1);
Byr=1967:2023;
cmbK=NaN(length(yr),2);

for i=1:length(yr)
    msk=yrs==yr(i) & mnth==9;
    if sum(msk)>0
        mnKC(i)=mean(mean(KC.ox(100:end,msk),'omitnan'));
    end
    cmbK(i,:)=[i find(Byr==yr(i))];
end

figure;
scatter(CUIw(cmb(:,2)),mnCS)
hold on
scatter(CUIw(cmb(end,2)),mnCS(end),10,'filled','r')
scatter(CUIw(cmbK(:,2)),mnKC)
hold on
scatter(CUIw(cmb(end,2)),mnKC(end-1),20,'filled','r')


%% look at April spiciness
% 
% msk=allProfs.missionIdx==1 & allProfs.canyonX'<50;
% figure;
% scatter(allProfs.mtime(msk),allProfs.isoSpice(:,msk),10,'filled')
% axdate; grid on;
% 
% mnth=month(allProfs.mtime);
% 
% % for i=200


% figure; 
% scatter(M.Scott2.time,M.Scott2.CTD280m.spice,10,'filled')
% hold on
% scatter(M.Hak1.time,M.Hak1.CTD133m.spice,10,'filled')
% axdate(24)


msk=allProfs.missionIdx==1 & allProfs.canyonX'<50 & allProfs.canyonX'>20 &...
    allProfs.mtime<datenum(2023,04,01);

[N, edgesX, edgesY] = histcounts2(allProfs.isoSpice(2,msk),allProfs.isoDO(2,msk), ...
    -0.8:0.05:0.25,30:5:350);

msk=allProfs.missionIdx==1 & allProfs.canyonX'<50 & allProfs.canyonX'>20 &...
    allProfs.mtime>datenum(2023,04,01) & allProfs.mtime<datenum(2023,05,01);

figure;
hold on
X = edgesX(1:end-1) + diff(edgesX)/2;
Y = edgesY(1:end-1) + diff(edgesY)/2;
p=pcolor(X, Y, N');p.FaceAlpha=0.8;shading flat; %colormap(gca,cm1);
% xlabel('$O$ ($\mu$mol kg$^{-1}$)','Interpreter','latex');
% ylabel('\sigma_\theta (kg m^{-3})');
axis tight;
scatter(allProfs.isoSpice(2,msk),allProfs.isoDO(2,msk),10,'r','filled')

%%
% export_fig /Users/samuelstevens/Dropbox/Hakai/gliders/oxygen/figures/AOUvsXV2.pdf -dpdf -nofontswap
 
 %% Old code
% % CS line CTD %
% axCS=axes('position',[0.1 0.15 0.5 0.2]);
% CSlat=CS_ctd.lat(:,CSmsk);
% CSlon=CS_ctd.lon(:,CSmsk);
% [~,msk]=max(CSlat);
% CSx=gsw_distance([ones(size(CSlon)).*CSlon(msk);CSlon]',...
%     [ones(size(CSlat)).*CSlat(msk);CSlat]')/1000;
% [CSx,msk]=sort(CSx);
% CSo=CS_ctd.ox(:,CSmsk);CSo=CS_ctd.ox(:,msk);
% CSz=CS_ctd.pr(:,CSmsk);
% 
% [xx,yy]=meshgrid(CSx,CSz(:,1));
% contourf(xx,yy,CSo,0:20:400,'LineStyle','none');
% clim([0 400]);
% colormap(axCS,cmocean('balance','pivotpoint',1.4*43.7));
% axis ij; ylim([0 500]);
% title(['Glider transect (',G.S33.t,')']);
% hold on
% [CS,CH]=contour(xx,yy,CSo,[1.4*43.7 1.4*43.7],'linewidth',0.75,'linecolor',[0.1 0.1 0.1]);
% % l=area(canyonX,canyonZ*-1,1000,'facecolor',[0.9 0.9 0.9],'facealpha',1,...
% %     'edgecolor',[0.4 0.4 0.4]);
% yticks(0:150:450); xticks(-30:30:120);
% xlabel('Along-canyon distance (km)');
% ylabel('Depth (m)');
% CB=m_contfbar([.1 .4],0.35,CS,CH,'axfrac',0.075,'edgecolor','none');
% yl=CB.YLim; axes(CB);
% line([1.4*43.7 1.4*43.7],yl,'color','k','linewidth',1);
% CB.XTick=[0 61.2 100:50:250];
% xlabel(CB,'D.O. (\mumol kg^{-1})');

% 
% % DENSITY zoom 1 %
% axMr=axes('position',[0.3 0.675 0.3 0.05]);
% p(2)=plot(M.Hak1.time,M.Hak1.CTD133m.sigma_theta,'color',lnes(2,:));
% hold on
% p(2).Color=[p(2).Color 0.075];
% s(2)=plot(M.Hak1.time,smoothdata(M.Hak1.CTD133m.sigma_theta,2,'movmean',wind,'omitnan'),...
%     'color',lnes(2,:));
% p(1)=plot(M.Scott2.time,M.Scott2.CTD280m.sigma_theta,'color',lnes(1,:));
% p(1).Color=[p(1).Color 0.075];
% s(1)=plot(M.Scott2.time,smoothdata(M.Scott2.CTD280m.sigma_theta,2,'movmean',wind,'omitnan'),...
%    'color',lnes(1,:));
% xlim([datenum(2022,04,01) datenum(2022,12,01)]);
% ylim([26.6 26.8]);
% axdate(10);xticklabels([]);
% grid on; box off;
% ylabel('Bottom \sigma_\theta (kg m^{-3})');
% 
% 
% % DENSITY plot %
% axMr=axes('position',[0.1 0.56 0.5 0.1]);
% % ll=line([datenum(2022,08,09) datenum(2022,08,09)],[25 175], ...
% %     'color','k','linewidth',0.75,'linestyle',':');
% hold on
% rectangle('Position',[datenum(2022,04,01) 26.6...
%     datenum(2022,12,01)-datenum(2022,04,01) 0.35],'edgecolor',[0.6 0.6 0.6]);
% p(1)=plot(M.Scott2.time,M.Scott2.CTD280m.sigma_theta,'color',lnes(1,:));
% p(1).Color=[p(1).Color 0.075];
% s(1)=plot(M.Scott2.time,smoothdata(M.Scott2.CTD280m.sigma_theta,2,'movmean',wind,'omitnan'),...
%    'color',lnes(1,:));
% p(2)=plot(M.Hak1.time,M.Hak1.CTD133m.sigma_theta,'color',lnes(2,:));
% hold on
% p(2).Color=[p(2).Color 0.075];
% s(2)=plot(M.Hak1.time,smoothdata(M.Hak1.CTD133m.sigma_theta,2,'movmean',wind,'omitnan'),...
%     'color',lnes(2,:));
% xlim([datenum(2022,04,01) datenum(2022,12,01)]);
% axdate(10);
% grid on; box off



% Density TS
% axMr=axes('position',[0.1 0.665 0.5 0.075]);
% % ll=line([datenum(2022,08,09) datenum(2022,08,09)],[25 175], ...
% %     'color','k','linewidth',0.75,'linestyle',':');
% hold on
% p(1)=plot(M.Scott2.time,M.Scott2.CTD280m.sigma_theta,'color',lnes(1,:));
% p(1).Color=[p(1).Color 0.075];
% s(1)=plot(M.Scott2.time,smoothdata(M.Scott2.CTD280m.sigma_theta,2,'movmean',wind,'omitnan'),...
%    'color',lnes(1,:));
% 
% xlim([datenum(2022,04,01) datenum(2022,12,01)]);
% axdate(10);
% xticklabels([]);
% grid on; box off;
% 
% axMr=axes('position',[0.1 0.575 0.5 0.075]);
% p(2)=plot(M.Hak1.time,M.Hak1.CTD133m.sigma_theta,'color',lnes(2,:));
% hold on
% p(2).Color=[p(2).Color 0.075];
% s(2)=plot(M.Hak1.time,smoothdata(M.Hak1.CTD133m.sigma_theta,2,'movmean',wind,'omitnan'),...
%     'color',lnes(2,:));
% xlim([datenum(2022,04,01) datenum(2022,12,01)]);
% % ylim([25 175]);
% axdate(10);
% grid on; box off;
% ylabel('Bottom \sigma_\theta (kg m^{-3})');


% plot(NS_ctd.ox(:,NSmsk1),NS_ctd.pr(:,NSmsk1),'linewidth',1.5,...
%     'color',rgb_x('pastel green'));
% plot(CS_ctd.ox(:,CSmsk(1)),CS_ctd.pr(:,CSmsk(1)),'linewidth',1.5,...
%     'color',rgb_x('light purple'));
% 
% axis tight ij;grid on; xlim([40 160]); ylim([50 500]);
% yl=ylim;
% line([1.4*43.7 1.4*43.7],yl,'color','k','linestyle','--');

% 
% %% What's the mean properties at S2 compared to offshore? 
% 
% mnth=month(M.Scott2.time);
% msk=mnth>=6 & mnth<=9;
% 
% mnST=mean(M.Scott2.CTD280m.sigma_theta(msk),'omitnan');
% mnT=mean(M.Scott2.CTD280m.temperature(msk),'omitnan');
% mnO=mean(M.Scott2.CTD280m.oxygen(msk),'omitnan')*43.7;
% mnS=mean(M.Scott2.CTD280m.salinity(msk),'omitnan');
% SA=gsw_SA_from_SP(M.Scott2.CTD280m.salinity,...
%     ones(size(M.Scott2.CTD280m.salinity)).*280,...
%     M.Scott2.longitude,M.Scott2.latitude);
% CT=gsw_CT_from_t(SA,M.Scott2.CTD280m.temperature,...
%     ones(size(M.Scott2.CTD280m.salinity)).*280);
% M.Scott2.CTD280m.rho=gsw_rho(SA,CT,ones(size(M.Scott2.CTD280m.salinity)).*280);
% mnR=mean(M.Scott2.CTD280m.rho(msk),'omitnan');
% 
% load CS02OS.mat
% tmpZ=mean(gsw_z_from_p(OSctd.pr,OSctd.lat),2,'omitnan')*-1;
% 
% tmpT=mean(OSctd.temp,2,'omitnan');
% tmpS=mean(OSctd.sal,2,'omitnan');
% tmpST=mean(OSctd.s_t,2,'omitnan');
% tmpDO=mean(OSctd.ox,2,'omitnan');
% tmpR=mean(OSctd.rho,2,'omitnan');
% 
% osT=interp1(tmpST(2:end),tmpT(2:end),mnST);
% osS=interp1(tmpST(2:end),tmpS(2:end),mnST);
% osO=interp1(tmpST(2:end),tmpDO(2:end),mnST);
% osZ=interp1(tmpST(2:end),tmpZ(2:end),mnST);

% osT=interp1(tmpR(2:end),tmpT(2:end),mnR);
% osS=interp1(tmpR(2:end),tmpS(2:end),mnR);
% % osO=interp1(tmpR(2:end),tmpDO(2:end),mnR);
% % osZ=interp1(tmpR(2:end),tmpZ(2:end),mnR);
% 
% %% Crawford style scatter
% 
% % for i=1:length(Cmsk)
% %     if bottomDO.oxygen(Cmsk(i))<0.5*43.7
% %         sm(1)=m_scatter(bottomDO.longitude(Cmsk(i)),bottomDO.latitude(Cmsk(i)),...
% %             20,'+','linewidth',1,'markeredgecolor','k');
% % 
% %     elseif bottomDO.oxygen(Cmsk(i))>=0.5*43.7 && bottomDO.oxygen(Cmsk(i))<=1*43.7
% %         sm(2)=m_scatter(bottomDO.longitude(Cmsk(i)),bottomDO.latitude(Cmsk(i)),...
% %             20,'+','linewidth',1,'markeredgecolor','r');
% %         m_scatter(bottomDO.longitude(Cmsk(i)),bottomDO.latitude(Cmsk(i)),...
% %             10,'.','linewidth',1,'markeredgecolor','r');
% % 
% %     elseif bottomDO.oxygen(Cmsk(i))>1*43.7 && bottomDO.oxygen(Cmsk(i))<=1.4*43.7
% %         sm(3)=m_scatter(bottomDO.longitude(Cmsk(i)),bottomDO.latitude(Cmsk(i)),...
% %             20,'filled','linewidth',1,'markeredgecolor','b');
% % 
% %     elseif bottomDO.oxygen(Cmsk(i))>1.4*43.7 && bottomDO.oxygen(Cmsk(i))<=3*43.7
% %         sm(4)=m_scatter(bottomDO.longitude(Cmsk(i)),bottomDO.latitude(Cmsk(i)),...
% %             20,'x','linewidth',1,'markeredgecolor','g');
% % 
% %     elseif bottomDO.oxygen(Cmsk(i))>3*43.7
% %         sm(5)=m_scatter(bottomDO.longitude(Cmsk(i)),bottomDO.latitude(Cmsk(i)),...
% %             20,'x','linewidth',1,'markeredgecolor',rgb_x('orange'));
% % 
% %     end
% % end
% % uistack(sm(5),'top'); uistack(sm(4),'top'); uistack(sm(3),'top'); uistack(sm(2),'top'); uistack(sm(1),'top');
% 
% 
%  %% Check individual deployments
% 
% filelist={'dfo-k999-20230915_grid_delayed.nc'};
% 
% fileInfo = ncinfo(filelist{1});
% a=struct();
% for i=1:length(fileInfo.Variables)
%     a.(fileInfo.Variables(i).Name)=ncread(filelist{1},...
%         fileInfo.Variables(i).Name);
% end
% 
% msk=a.longitude>-129.5;
% figure;
% plot(a.oxygen_corrected(msk),a.depth(msk));
% axis ij; ylim([0 500]);
% line([60 60],[0 400],'color','k','linewidth',2);
% 
% 
% filelist={'dfo-eva035-20230620.nc'};
% 
% fileInfo = ncinfo(filelist{1});
% a=struct();
% for i=1:length(fileInfo.Variables)
%     a.(fileInfo.Variables(i).Name)=ncread(filelist{1},...
%         fileInfo.Variables(i).Name);
% end
% 
% msk=a.longitude>-129.5;
% figure;
% plot(a.oxygen_corrected(msk),a.depth(msk));
% axis ij; ylim([0 500]);
% line([60 60],[0 400],'color','k','linewidth',2);
% 
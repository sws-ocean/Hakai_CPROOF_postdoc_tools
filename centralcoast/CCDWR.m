%% Mean properties and seasonality (Fitzhugh)

addpath(genpath('/Users/samst/Dropbox/Hakai/'));
addpath(genpath('/Users/samst/Dropbox/UBC/GSW/'));
addpath(genpath('/Users/samst/Dropbox/UBC/m_map/'));
addpath(genpath('/Users/samst/Dropbox/UBC/Misc/'));

clear 
load allCalvert_wStations.mat
H = HakaiWaterPropertiesInstrumentP;
clearvars HakaiWaterPropertiesInstrumentP
H.pressure=double(H.pressure);


H.SA=gsw_SA_from_SP(H.salinity,H.pressure,H.longitude,H.latitude);
H.CT=gsw_CT_from_t(H.SA,H.temperature,H.pressure);
H.s_t=gsw_sigma0(H.SA,H.CT);
H.spice=gsw_spiciness0(H.SA,H.CT);
H.dissolved_oxygen_ml_l=H.dissolved_oxygen_ml_l.*43.7;
H.time = datenum(1970,01,01,00,00,00) + (H.time / 86400);

load H1.mat;
load S2.mat

%% work through Hakai datasets
sig_thresh = [25.5];
H.station = string(H.station);
allStations=unique(H.station);
allHak=struct();

sanitizedStations = matlab.lang.makeValidName(allStations);


for i=1:length(allStations)
    idx=strcmp(H.station,allStations(i));
    
    zs=H.depth(idx);
    ts=H.temperature(idx);
    ss=H.salinity(idx);
    os=H.dissolved_oxygen_ml_l(idx);
    sts=H.s_t(idx);
    tts=H.time(idx);
    spi=H.spice(idx);

    profStop=[find(sign(diff(zs))==-1);length(zs)];
    
    if length(profStop)>10

        allHak.(sanitizedStations(i)).lon=mean(H.longitude(idx),'omitnan');
        allHak.(sanitizedStations(i)).lat=mean(H.latitude(idx),'omitnan');

        for ii=1:length(profStop)-1
            
            if ii==1
                profRng=profStop(1):profStop(2);
            else
                profRng=profStop(ii)+1:profStop(ii+1);
            end

            if length(profRng)>1

                allHak.(sanitizedStations(i)).time(ii)=mean(tts(profRng),'omitnan');
                allHak.(sanitizedStations(i)).dday(ii)=...
                    day(datetime(datestr(allHak.(sanitizedStations(i)).time(ii))),'dayofyear');


                fuzz=rand(length(profRng),1)/1000;
                allHak.(sanitizedStations(i)).depth(:,ii)=[1:500]';
                allHak.(sanitizedStations(i)).t(:,ii)=...
                    interp1(zs(profRng)+fuzz,ts(profRng),1:500);
                allHak.(sanitizedStations(i)).s(:,ii)=...
                    interp1(zs(profRng)+fuzz,ss(profRng),1:500);
                allHak.(sanitizedStations(i)).o(:,ii)=...
                    interp1(zs(profRng)+fuzz,os(profRng),1:500);
                allHak.(sanitizedStations(i)).s_t(:,ii)=...
                    interp1(zs(profRng)+fuzz,sts(profRng),1:500);
                allHak.(sanitizedStations(i)).spice(:,ii)=...
                    interp1(zs(profRng)+fuzz,spi(profRng),1:500);
                
                try
                    % calculate isopycnal temperature
                    allHak.(sanitizedStations(i)).IPtemp(:,ii)=...
                        interp1(sts(profRng)+fuzz,ts(profRng),sig_thresh);

                    % calculate isopycnal oxygen
                    allHak.(sanitizedStations(i)).IPox(:,ii)=...
                        interp1(sts(profRng)+fuzz,os(profRng),sig_thresh);
                    
                    % calculate isopycnal salinity
                    allHak.(sanitizedStations(i)).IPsal(:,ii)=...
                        interp1(sts(profRng)+fuzz,ss(profRng),sig_thresh);

                    % calculate isopycnal spice
                    allHak.(sanitizedStations(i)).IPspice(:,ii)=...
                        interp1(sts(profRng)+fuzz,spi(profRng),sig_thresh);

                    % calculate isopycnal depth
                    allHak.(sanitizedStations(i)).IPz(:,ii)=...
                        interp1(sts(profRng)+fuzz,zs(profRng),sig_thresh);
    
                     allHak.(sanitizedStations(i)).Ob(:,ii)=os(profRng(end));

                end
            end
        end
    end
end


%% Pruth
load HakaiPruthDockProvisional.mat
PD=HakaiPruthDockProvisional;

load BAKmay2024.mat  

PD.mtime=datenum(1970,01,01,00,00,00)+(PD.time/86400);
PD.dailyT=floor(min(PD.mtime)):floor(max(PD.mtime));
PD.daily=NaN(1,length(PD.dailyT));

count=0;
for i=PD.dailyT(1):PD.dailyT(end)
    count=count+1;

    msk=floor(PD.mtime)==i;
    PD.daily(count)=mean(PD.tideheight_avg(msk),'omitnan');
end

PD.daily=PD.daily-mean(PD.daily,'omitnan');

%% Find density increases at KC10

tgrid=floor(allHak.KC10.time(1)):floor(allHak.KC10.time(end));
KCi=interp1(allHak.KC10.time,allHak.KC10.s_t(250,:),tgrid);
% KCi=interp1([allHak.KC10.time allHak.FZH08.time],...
%     [allHak.KC10.s_t(250,:) allHak.FZH08.s_t(300,:)],tgrid);

KCi=inpaint_nans(KCi);

UPmsk=diff(KCi)>0;

%%
figure('units','centimeters','outerposition',[0 0 30 30],'color','w');
allHak.KC10.s_t(250,allHak.KC10.s_t(250,:)<24)=NaN;

lnes=linspecer(5);

ax1=subplot(3,1,1); hold on;
line([tgrid(UPmsk)' tgrid(UPmsk)'],[25.2 26.9],'color',[0 0 0 0.02])
le(1)=plot(H1.time,smoothdata(H1.CTD133m.sigma_theta,'movmean',24*7,'omitnan'),...
    'color',[0 0 0 0.3]);
le(2)=plot(S2.time,smoothdata(S2.CTD280m.sigma_theta,'movmean',24*7,'omitnan'),...
    'color',[0 0 0 1]);
le(3)=plot(allHak.KC10.time,allHak.KC10.s_t(250,:),'color',lnes(1,:));
scatter(allHak.KC10.time,allHak.KC10.s_t(250,:),20,'filled','markerfacecolor',lnes(1,:));
le(4)=plot(allHak.FZH08.time,allHak.FZH08.s_t(300,:),'color',lnes(2,:));
scatter(allHak.FZH08.time,allHak.FZH08.s_t(300,:),20,'filled','markerfacecolor',lnes(2,:));
le(5)=plot(allHak.DFO3.time,allHak.DFO3.s_t(275,:),'color',lnes(3,:));
scatter(allHak.DFO3.time,allHak.DFO3.s_t(275,:),20,'filled','markerfacecolor',lnes(3,:));
% plot(allHak.FZH08.time,allHak.FZH08.s_t(300,:),'color',lnes(2,:));
% scatter(allHak.FZH08.time,allHak.FZH08.s_t(300,:),20,'k','filled','markerfacecolor',lnes(2,:));
axis tight;
xlim([datenum(2013,06,01) datenum(2025,01,01)]);
axdate; grid on; 
le(6)=fill(NaN,NaN,rgb_x('light grey'),'edgecolor','none');
ylabel('\sigma_\theta (kg m^{-3})');


ax2=subplot(3,1,2); hold on
line([tgrid(UPmsk)' tgrid(UPmsk)'],[-275 90],'color',[0 0 0 0.02])

plot(BAK.mtime,smoothdata(BAK.B,'movmean',14),'k');
ylabel({'Bakun index @ 51N';'(m^3 s^{-1} 100m^{-1})'});
set(gca,'ydir','reverse');
ylim([-296.2617  103.7383]);

yyaxis right
plot(PD.dailyT,smoothdata(PD.daily,'movmean',14),'color',lnes(end,:));
xlim([datenum(2013,06,01) datenum(2025,01,01)]);
axdate; grid on; 
ylabel('SSHA @ Pruth Bay (m)')
ylim([-0.3411  0.6589])
set(gca,'ycolor',lnes(end,:));

subplot(3,1,3);
load('BCcoast');
lat_lim=[50.5 52]; lon_lim=[-131.5 -127];
col=[255 214 140]/255; % YELLOW!
m_proj('equidistant','lon',lon_lim,'lat',lat_lim);   % Projection
hold on
[CS,CH]=m_etopo2('contour',[-1000 -300 -100],'edgecolor',rgb_x('light grey'));
clabel(CS,CH,'labelspacing',200,'color',rgb_x('grey'),'fontsize',6);
m_gshhs_f('patch',col,'edgecolor',[0.3 0.3 0.3],'linewidth',0.3);
m_grid('linestyle','none','linewidth',1,'tickdir','out');
m_scatter(-128.129,51.655,100,'filled','markerfacecolor',...
    lnes(end,:),'markeredgecolor','k','marker','v');
m_text(-128.129-0.9,51.655,"Pruth Bay",'fontweight','bold');
m_scatter(-131,51,100,'x','markeredgecolor','k','linewidth',2);
m_text(-131-0.3,51+0.25,{"Bakun Index";"(131W, 51N)"},'fontweight','bold');

m_scatter(allHak.KC10.lon,allHak.KC10.lat,100,'filled',...
    'markerfacecolor',lnes(1,:));
m_text(allHak.KC10.lon,allHak.KC10.lat,{'KC10'},'fontweight','bold');

m_scatter(allHak.FZH08.lon,allHak.FZH08.lat,100,'filled',...
    'markerfacecolor',lnes(2,:));
m_text(allHak.FZH08.lon,allHak.FZH08.lat,{'FZH08'},'fontweight','bold');

m_scatter(allHak.DFO3.lon,allHak.DFO3.lat,100,'filled',...
    'markerfacecolor',lnes(3,:));
m_text(allHak.DFO3.lon,allHak.DFO3.lat,{'DFO3'},'fontweight','bold');


set(findall(gcf,'-property','fontsize'),'fontsize',12);
axes(ax1)
legend(le,'Hak1 (133m)','Scott2 (280m)','KC10 (250m)','FZH08 (300m)',...
    'DFO3 (275m)','Deep Renewal Period','fontsize',8);
set(findall(gcf,'-property','fontname'),'fontname','Latin Modern Roman');

%%
export_fig /Users/samuelstevens/Dropbox/Hakai/synthesis/figs/CCDWR.jpg -m5 -nofontswap


%% Temp and O2

figure('units','centimeters','outerposition',[0 0 30 30],'color','w');

allHak.KC10.o(allHak.KC10.o>500)=NaN;
allHak.DFO3.o(allHak.DFO3.o>500)=NaN;
allHak.FZH08.o(allHak.FZH08.o>500)=NaN;
allHak.FZH08.o(300,allHak.FZH08.o(300,:)>250)=NaN;

lnes=linspecer(5);

ax1=subplot(3,1,1); hold on;
line([tgrid(UPmsk)' tgrid(UPmsk)'],[5.25 8.75],'color',[0 0 0 0.02])
le(1)=plot(H1.time,smoothdata(H1.CTD133m.temperature,'movmean',24*7,'omitnan'),...
    'color',[0 0 0 0.3]);
le(2)=plot(S2.time,smoothdata(S2.CTD280m.temperature,'movmean',24*7,'omitnan'),...
    'color',[0 0 0 1]);
le(3)=plot(allHak.KC10.time,allHak.KC10.t(250,:),'color',lnes(1,:));
scatter(allHak.KC10.time,allHak.KC10.t(250,:),20,'filled','markerfacecolor',lnes(1,:));
le(4)=plot(allHak.FZH08.time,allHak.FZH08.t(300,:),'color',lnes(2,:));
scatter(allHak.FZH08.time,allHak.FZH08.t(300,:),20,'filled','markerfacecolor',lnes(2,:));
le(5)=plot(allHak.DFO3.time,allHak.DFO3.t(275,:),'color',lnes(3,:));
scatter(allHak.DFO3.time,allHak.DFO3.t(275,:),20,'filled','markerfacecolor',lnes(3,:));
% plot(allHak.FZH08.time,allHak.FZH08.s_t(300,:),'color',lnes(2,:));
% scatter(allHak.FZH08.time,allHak.FZH08.s_t(300,:),20,'k','filled','markerfacecolor',lnes(2,:));
axis tight;
xlim([datenum(2013,06,01) datenum(2025,01,01)]);
axdate; grid on; 
le(6)=fill(NaN,NaN,rgb_x('light grey'),'edgecolor','none');
ylabel('Temperature (^oC)');

ax2=subplot(3,1,2); hold on;
line([tgrid(UPmsk)' tgrid(UPmsk)'],[25 205],'color',[0 0 0 0.02])
plot(H1.time,smoothdata(H1.CTD133m.oxygen,'movmean',24*7,'omitnan'),...
    'color',[0 0 0 0.3]);
plot(S2.time,smoothdata(S2.CTD280m.oxygen,'movmean',24*7,'omitnan'),...
    'color',[0 0 0 1]);
plot(allHak.KC10.time,allHak.KC10.o(250,:),'color',lnes(1,:));
scatter(allHak.KC10.time,allHak.KC10.o(250,:),20,'filled','markerfacecolor',lnes(1,:));
plot(allHak.FZH08.time,allHak.FZH08.o(300,:),'color',lnes(2,:));
scatter(allHak.FZH08.time,allHak.FZH08.o(300,:),20,'filled','markerfacecolor',lnes(2,:));
plot(allHak.DFO3.time,allHak.DFO3.o(275,:),'color',lnes(3,:));
scatter(allHak.DFO3.time,allHak.DFO3.o(275,:),20,'filled','markerfacecolor',lnes(3,:));
% plot(allHak.FZH08.time,allHak.FZH08.s_t(300,:),'color',lnes(2,:));
% scatter(allHak.FZH08.time,allHak.FZH08.s_t(300,:),20,'k','filled','markerfacecolor',lnes(2,:));
axis tight;
xlim([datenum(2013,06,01) datenum(2025,01,01)]);
axdate; grid on; 
ylabel('O$_2$ ($\mu$mol kg$^{-1}$)', 'Interpreter', 'latex');



set(findall(gcf,'-property','fontsize'),'fontsize',12);
axes(ax1)
legend(le,'Hak1 (133m)','Scott2 (280m)','KC10 (250m)','FZH08 (300m)',...
    'DFO3 (275m)','Deep Renewal Period','fontsize',8);
set(findall(gcf,'-property','fontname'),'fontname','Latin Modern Roman');


%%
export_fig /Users/samuelstevens/Dropbox/Hakai/synthesis/figs/CCDWR_TnO.jpg -m5 -nofontswap


%% Map isopycnal properties

figure('units','centimeters','outerposition',[0 0 30 30],'color','w');
lat_lim2=[51.35 52.5];  lon_lim2=[-128.4 -126.75];
col=[255 214 140]/255; % YELLOW
fname='/Users/samst/Dropbox/UBC/Misc/british_columbia_3_msl_2013.nc';
lat=ncread(fname,'lat');    
lon=ncread(fname,'lon');
ilon=lon>=lon_lim2(1) & lon<=lon_lim2(2);
ilat=lat>=lat_lim2(1) & lat<=lat_lim2(2);
Z=ncread(fname,...
    'Band1',[ find((ilon),1,'first') find((ilat),1,'first')],...
    [ sum(ilon) sum(ilat)],[1 1]);
m_proj('equidistant','lon',lon_lim2,'lat',lat_lim2);   % Projection
hold on
% [CS,CH]=m_contour(lon(ilon),lat(ilat),Z',[-600:100:0]);
m_gshhs_f('patch',col,'edgecolor',[0.9 0.9 0.9],'linewidth',0.1);
m_grid('tickdir','out','linestyle','none','YaxisLocation','left','fontsize',10);

flds=fieldnames(allHak);
tmpO=NaN(length(flds),1);l=NaN(length(flds),1);ll=NaN(length(flds),1);
for i=1:length(flds)
    try
        if length(allHak.(flds{i}).IPox)>12
            tmpO(i)=mean(allHak.(flds{i}).IPox,'omitnan');
            l(i)=mean(allHak.(flds{i}).lon,'omitnan');
            ll(i)=mean(allHak.(flds{i}).lat,'omitnan');
        end
    end
end

m_scatter(l,ll,100,tmpO','filled');
clim([100 175])


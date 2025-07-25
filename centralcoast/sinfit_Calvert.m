%% Is there seasonality in QCS marginal deep water, and can we use it as a tracer? 

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



%% seasonality: isopycnal oxygen
flds = fieldnames(allHak);

a=[];
for i=1:length(flds)
    if isfield(allHak.(flds{i}),'IPtemp')

        % fit sinusoid to seasonal cycle all all isopycnals
        for iii=1:length(sig_thresh)
                        
            msk=~isnan(allHak.(flds{i}).IPox(iii,:)) & ...
                allHak.(flds{i}).IPox(iii,:)~=0 & ...
                allHak.(flds{i}).IPox(iii,:)<500;
            cleanO=allHak.(flds{i}).IPox(iii,msk);
            dday=allHak.(flds{i}).dday(msk);

            [amp,phase,frac,offset,yy]=fit_harmonics(...
                cleanO,dday,1,365,0.01);
            yt=offset;
            for tt=1:365
                fitted_harm(tt)=offset+amp*cos(2*pi*tt/365 + phase);
            end
            
            [~, allHak.(flds{i}).cold_dayO(iii,1)]=max(fitted_harm);
            allHak.(flds{i}).ampO(iii,1)=amp;
            allHak.(flds{i}).Os=[cleanO;dday];

            % if allHak.(flds{i}).cold_dayO(iii,1)~=1
            %     figure('units','centimeters','outerposition',[0 0 9 9],'color','w');
            %     plot(fitted_harm,'k')
            %     hold on
            %     scatter(dday,cleanO,10,'filled','markerfacecolor',rgb_x('rose'));
            % 
            %     grid on;
            %     axis tight; ylim([50 200]);
            %     line([allHak.(flds{i}).cold_dayO(iii,1) allHak.(flds{i}).cold_dayO(iii,1)],...
            %         [0 max(fitted_harm)],'color','k','linestyle','--');
            %     title(flds{i});
            %     text(0.5,0.2,sprintf('z = %0.0f - %0.0f m\nphase = %0.0f days\namp = %0.2f umol/kg',...
            %         min(allHak.(flds{i}).IPz(iii,msk)),...
            %         max(allHak.(flds{i}).IPz(iii,msk)),...
            %         allHak.(flds{i}).cold_dayO(iii,1),amp),...
            %         'units','normalized');
            %     ylabel('O^2 (\mumol kg^{-1})');
            %     xlabel('Julian Day');
            %     set(findall(gcf, '-property', 'fontsize'), 'fontsize', 10);
            %     set(findall(gcf, '-property', 'fontname'), 'fontname', 'Latin Modern Roman');
            %     exportgraphics(gcf, ['/Users/samuelstevens/Dropbox/Hakai/synthesis/figs/sinfit/',flds{i},'.jpg'], ...
            %         'Resolution', 400);
            %     close all;
            % end

            try
                a=[a;min(allHak.(flds{i}).IPz(1,msk)) max(allHak.(flds{i}).IPz(1,msk)) nanmean(allHak.(flds{i}).IPz(1,msk))];
            end
        end

        % allHak.(flds{i}).cold_dayO(allHak.(flds{i}).cold_dayO==1)=NaN;
    end
end

% cold_dayO(cold_dayO<=1)=NaN;


%% make some maps

figure('units','centimeters','outerposition',[0 0 30 30],'color','w');
lat_lim2=[51.35 52];  lon_lim2=[-128.4 -127.25];
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
[CS,CH]=m_contourf(lon(ilon),lat(ilat),Z',[-1000:50:-100 0],'linestyle','none');
colormap(m_colmap('blues'));
clim([-500 0]);
ax=colorbar;ax.YLim=[-500 0];ax.FontSize=12;
title(ax,'(m)','fontweight','normal');
m_gshhs_f('patch',col,'edgecolor',[0.9 0.9 0.9],'linewidth',0.1);
m_grid('tickdir','out','linestyle','none','YaxisLocation','left','fontsize',14,...
    'xticklabels',[],'yticklabels',[]);

for i=1:length(flds)
    if isfield(allHak.(flds{i}),'cold_dayO')
        if allHak.(flds{i}).cold_dayO>1 && length(allHak.(flds{i}).Os)>10
            m_text(allHak.(flds{i}).lon,allHak.(flds{i}).lat,...
                sprintf('%0.0f/%2.0f',allHak.(flds{i}).cold_dayO,...
                allHak.(flds{i}).ampO),'fontsize',7,'color','k','fontweight','bold');
        end
    end
end
set(gca,'clipping','on');

set(findall(gcf, '-property', 'fontname'), 'fontname', 'Latin Modern Roman');

%%
% exportgraphics(gcf, '/Users/samuelstevens/Dropbox/Hakai/synthesis/figs/seasonalityCalvert.jpg', ...
%      'Resolution', 500);

%%
figure('units','centimeters','outerposition',[0 0 30 30],'color','w');
lat_lim2=[51.35 52.5];  lon_lim2=[-128.4 -127.1];
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
[CS,CH]=m_contourf(lon(ilon),lat(ilat),Z',[-1000:50:-100 0],'linestyle','none');
colormap(m_colmap('blues'));
clim([-500 0]);
ax=colorbar;ax.YLim=[-500 0];ax.FontSize=12;
title(ax,'(m)','fontweight','normal');
m_gshhs_f('patch',col,'edgecolor',[0.9 0.9 0.9],'linewidth',0.1);
m_grid('tickdir','out','linestyle','none','YaxisLocation','left','fontsize',14,...
    'xticklabels',[],'yticklabels',[]);

for i=1:length(flds)
    if isfield(allHak.(flds{i}),'cold_dayO')
        if allHak.(flds{i}).cold_dayO>1 && length(allHak.(flds{i}).Os)>10
            m_text(allHak.(flds{i}).lon,allHak.(flds{i}).lat,...
                sprintf('%0.0f/%2.0f',allHak.(flds{i}).cold_dayO,...
                allHak.(flds{i}).ampO),'fontsize',7,'color','k','fontweight','bold');
        end
    end
end
set(gca,'clipping','on');

set(findall(gcf, '-property', 'fontname'), 'fontname', 'Latin Modern Roman');

%%
% exportgraphics(gcf, '/Users/samuelstevens/Dropbox/Hakai/synthesis/figs/seasonalityCalvertBnD.jpg', ...
%      'Resolution', 500);


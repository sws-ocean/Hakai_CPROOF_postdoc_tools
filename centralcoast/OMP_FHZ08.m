%% OMP analysis- Hakai Pass

addpath(genpath('/Users/samst/Dropbox/Hakai/'));
addpath(genpath('/Users/samst/Dropbox/UBC/GSW/'));
addpath(genpath('/Users/samst/Dropbox/UBC/m_map/'));
addpath(genpath('/Users/samst/Dropbox/UBC/Misc/'));

%% Load observations
clear
load allCalvert_wStations.mat
H = HakaiWaterPropertiesInstrumentP;
clearvars HakaiWaterPropertiesInstrumentP
H.pressure=double(H.pressure);

H.SA=gsw_SA_from_SP(H.salinity,H.pressure,H.longitude,H.latitude);
H.CT=gsw_CT_from_t(H.SA,H.temperature,H.pressure);
H.s_t=gsw_sigma0(H.SA,H.CT);
H.pt=gsw_pt_from_CT(H.SA,H.CT);
H.spice=gsw_spiciness0(H.SA,H.CT);
H.dissolved_oxygen_ml_l=H.dissolved_oxygen_ml_l;
H.dissolved_oxygen_ml_l(H.dissolved_oxygen_ml_l>20)=NaN;
H.time = datenum(1970,01,01,00,00,00) + (H.time / 86400);

%% work through Hakai datasets
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
    pt=H.pt(idx);
    cd=H.cdom_ppb(idx);

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
                allHak.(sanitizedStations(i)).pt(:,ii)=...
                    interp1(zs(profRng)+fuzz,pt(profRng),1:500);
                allHak.(sanitizedStations(i)).s(:,ii)=...
                    interp1(zs(profRng)+fuzz,ss(profRng),1:500);
                allHak.(sanitizedStations(i)).o(:,ii)=...
                    interp1(zs(profRng)+fuzz,os(profRng),1:500);
                allHak.(sanitizedStations(i)).s_t(:,ii)=...
                    interp1(zs(profRng)+fuzz,sts(profRng),1:500);
                allHak.(sanitizedStations(i)).spice(:,ii)=...
                    interp1(zs(profRng)+fuzz,spi(profRng),1:500);
                allHak.(sanitizedStations(i)).cdom_ppb(:,ii)=...
                    interp1(zs(profRng)+fuzz,cd(profRng),1:500);
            end
        end
    end
end

save('/Users/samuelstevens/Dropbox/Hakai/ctd/data/allHakSanitized','allHak');

%% prep FHZ08 data and cononfigure observations for OMP

press=allHak.FZH08.depth(:);       % [dbar]      
ptemp=allHak.FZH08.pt(:);   % [°C]        
sal=allHak.FZH08.s(:) ;%[PSU]       
oxy=allHak.FZH08.o(:);      % [ml/l]
long=allHak.FZH08.lon.*ones(length(press),1);
lat=allHak.FZH08.lat.*ones(length(press),1);

time=[];
for i=1:length(allHak.FZH08.time);
    time=[time;allHak.FZH08.time(i).*ones(500,1)];
end

nanmsk=isnan(press) | isnan(ptemp) |isnan(sal) |isnan(oxy);

press(nanmsk)=[]; % Can use -9 as NaNs
ptemp(nanmsk)=[];
sal(nanmsk)=[];
oxy(nanmsk)=[];
lat(nanmsk)=[];
long(nanmsk)=[];
time(nanmsk)=[];
pdens = sw_pden(sal, ptemp, press, 0) - 1000;

save('/Users/samuelstevens/Dropbox/Hakai/synthesis/data/FZH08_ompdata.mat',...
    'lat', 'long', 'press', 'ptemp', 'sal','oxy','pdens');

%% Quickly visualize the selected data
f1=figure('units','centimeters','outerposition',...
    [0.01 0.01 19 8],'color','w');

tmpT=allHak.FZH08.pt;
tmpS=allHak.FZH08.s;
tmpO=allHak.FZH08.o;
tmpP=allHak.FZH08.depth;
tmpX=ones(size(tmpP)).*allHak.FZH08.time;

subplot(3,1,1);
contourf(tmpX,tmpP,tmpT,...
    5:0.1:12,'linecolor','none');
% clim([5 10]);
axis ij;
colormap(gca,cmocean('thermal'));
c=colorbar('location','eastoutside');
axdate(10);
ylim([0 400]);
ylabel('Depth (m)');
ylabel(c,'Temprature (^oC)');

subplot(3,1,2);
contourf(tmpX,tmpP,tmpS,...
    24:0.05:34,'linecolor','none');
clim([31 34]);
axis ij;
colormap(gca,cmocean('haline'));
c=colorbar('location','eastoutside');
title('FZH08')
axdate(10);
ylim([0 400]);
ylabel('Depth (m)')
ylabel(c,'Salinity (PSU)');

subplot(3,1,3);
contourf(tmpX,tmpP,tmpO,...
    1:0.1:7,'linecolor','none');
% clim([6 9]);
axis ij;
colormap(gca,cmocean('amp'));
c=colorbar('location','eastoutside');
hold on;
axdate(10);
ylim([0 400]);
ylabel('Depth (m)')
ylabel(c,'DO (ml/l)');


set(findall(gcf, '-property', 'fontname'), 'fontname', 'Latin Modern Roman');
set(findall(gcf, '-property', 'fontsize'), 'fontsize', 12);


%%
export_fig /Users/samuelstevens/Dropbox/Hakai/synthesis/figs/OMP_input.jpg -m5 -nofontswap

%% Calculate weights 
% Use inverse variance between water masses (Masson 2006)

qwt = qwt_shelfMn([1 2 3]); 
tracer_std = std(qwt); 
weights = 1 ./ tracer_std;

% Optional: Normalize so largest weight = 1
% weights = weights / max(weights);

disp('Weights (normalized):');
disp(weights);

load testwght.mat

Wx(1,1)=weights(1);
Wx(2,2)=weights(2);
Wx(3,3)=weights(3);

% Wx(:,5:8)=0.;

save('/Users/samuelstevens/Dropbox/UBC/Misc/omp2/omp2/shelfMnwght.mat','ratio','Wx');

% save('/Users/samuelstevens/Dropbox/Hakai/synthesis/data/HKPweights.mat','weights');

%% save workspace
clearvars allHak
close all
save('/Users/samuelstevens/Dropbox/Hakai/synthesis/data/FZH08_ompWS.mat');

%% Run OMP
omp2gui;

%% Plot some results
clear
load FZH08_ompWS
load test.mat

%% Do some plotting
f1=figure('units','centimeters','outerposition',...cl
    [0.01 0.01 20 30],'color','w');
old_data=load('FZH08_ompdata.mat');
time(isnan(old_data.ptemp))=[];

for jj=1:size(A,1)
    para=A(jj,:)*100;

    % create regular grid:
    XI=linspace(min(time),max(time),200)';
    YI=linspace(min(press),max(press),400);
    
    % interpolate to regular grid:
    para2=griddata(time(1:10:end),press(1:10:end)',para(1:10:end),XI,YI,'v4');
    para3=griddata(time(1:10:end),press(1:10:end)',pdens(1:10:end),XI,YI,'v4');

    subplot(size(A,1),1,jj);

    contourf(XI,YI,para2,[0:10:100],'linestyle',...
        'none');
    hold on
    contour(XI,YI,para3,[20:1:27],'linecolor',[0.5 0.5 0.5]);

    if jj==1
        % colormap(gca,cmocean('amp'));
        title(tit_index(1:5))
    elseif jj==2
        % colormap(gca,flipud(cmocean('ice')));
        title(tit_index(6:10));
    elseif jj==3
        % colormap(gca,cmocean('grey'));
        title(tit_index(11:15));
    else
        colormap(cmocean('amp'));
    end

    colormap(cmocean('amp'));
    clim([0 100]);
    c=colorbar;
    hold on
    drawnow
    axis ij
    % xlim([0 38]);
    ylim([0 385]);
    ylabel('Depth (m)');
    axdate(10);
    title(c,'Water Mass Content (%)')
end
set(findall(gcf, '-property', 'fontname'), 'fontname', 'Latin Modern Roman');
set(findall(gcf, '-property', 'fontsize'), 'fontsize', 10);

%%
export_fig /Users/samuelstevens/Dropbox/Hakai/synthesis/figs/OMP_SW123_output.jpg -m5 -nofontswap

%% Do some plotting-density space
f1=figure('units','centimeters','outerposition',...
    [0.01 0.01 8 30],'color','w');
for jj=1:size(A,1)
    para=A(jj,:)*100;

    % create regular grid:
    XI=linspace(min(time),max(time),200)';
    YI=linspace(min(pdens),max(pdens),400);
    
    % interpolate to regular grid:
    para2=griddata(time(1:10:end),pdens(1:10:end)',para(1:10:end),XI,YI,'v4');

    subplot(size(A,1),1,jj);

    contourf(XI,YI,para2,[0:10:100],'linestyle',...
        'none');

    if jj==1
        % colormap(gca,cmocean('amp'));
        title(tit_index(1:5))
    elseif jj==2
        % colormap(gca,flipud(cmocean('ice')));
        title(tit_index(6:10));
    elseif jj==3
        % colormap(gca,cmocean('grey'));
        title(tit_index(11:15));
    else
        colormap(cmocean('amp'));
    end

    colormap(cmocean('amp'));
    clim([0 100]);
    c=colorbar;
    hold on
    drawnow
    axis ij
    % xlim([0 38]);
    ylim([24 26.5]);
    ylabel('Depth (m)');
    axdate(10);
    title(c,'Water Mass Content (%)')
end
set(findall(gcf, '-property', 'fontname'), 'fontname', 'Latin Modern Roman');
set(findall(gcf, '-property', 'fontsize'), 'fontsize', 10);
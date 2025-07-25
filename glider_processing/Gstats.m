%% Glider stats
addpath(genpath('/Users/samst/Dropbox/Hakai/'));
addpath(genpath('/Users/samst/Dropbox/UBC/GSW/'));
addpath(genpath('/Users/samst/Dropbox/UBC/m_map/'));
addpath(genpath('/Users/samst/Dropbox/UBC/Misc/'));

clear
load glider_allProfs.mat

%% Get some numbers
clc
missions=max(unique(allProfs.missionNum(:)));

fprintf('%5i glider surveys\n',missions);
fprintf('%5i glider profiles\n',size(allProfs.temperature,2));

% how far apart are the profiles? How often?
msk=canyonX(allProfs.canyonIDX)'~=1 & canyonX(allProfs.canyonIDX)'~=max(canyonX);
X=gsw_distance(allProfs.longitude(msk),allProfs.latitude(msk))/1000;
T=diff(allProfs.mtime(msk));

fprintf('%5.1f km average spacing\n',mean(X(X<10),'omitnan'));
fprintf('%5.0f minutes average frequency\n',mean(T(T<1),'omitnan')*3600);

figure('units','centimeters','outerposition',[0 0 15 10],'color','w');
subplot(1,2,1)
histogram(T(T<0.2)*24)
xlabel('Time between profiles (hrs)');
grid on; axis tight;
subplot(1,2,2)
histogram(X(X<2));
xlabel('Distance between profiles (kms)');
grid on; axis tight;

% export_fig /Users/samuelstevens/Dropbox/Hakai/gliders/figs/resolution.png -m3 -nofontswap

%% time vs canyon axis plot

figure('units','normalized','position',[0 0 0.75 0.75],'color','w');

hold on
lnes=lines(missions);

for i=1:6
    subplot(6,1,i)
    hold on
    for ii=1:missions
        msk=allProfs.missionNum==ii & canyonX(allProfs.canyonIDX)'~=1 ...
            & canyonX(allProfs.canyonIDX)'~=max(canyonX);
        scatter(allProfs.mtime(msk),canyonX(allProfs.canyonIDX(msk)),...
            40,'filled','MarkerFaceColor',lnes(ii,:));

        msk=allProfs.missionNum==ii & canyonX(allProfs.canyonIDX)'==1 ...
            | canyonX(allProfs.canyonIDX)'==max(canyonX);
        scatter(allProfs.mtime(msk),canyonX(allProfs.canyonIDX(msk)),...
            20,'filled','MarkerFaceColor',[0.8 0.8 0.8]);
    end
    xlim([datenum(2018+i,01,01) datenum(2019+i,01,01)])
    axdate; grid on;

end
ylabel('Distance along canyon (km)');
set(findall(gcf,'-property','fontsize'),'fontsize',16);
set(findall(gcf,'-property','fontname'),'fontname','Latin Modern Roman');
%%
% export_fig /Users/samuelstevens/Dropbox/Hakai/gliders/figs/timeSpace.png -m3 -nofontswap
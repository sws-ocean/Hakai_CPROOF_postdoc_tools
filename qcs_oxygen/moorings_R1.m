%% Mooring climatologies and PPs
% 
% addpath(genpath('/Users/samst/Dropbox/Hakai/'));
% addpath(genpath('/Users/samst/Dropbox/UBC/GSW/'));
% addpath(genpath('/Users/samst/Dropbox/UBC/m_map/'));
% addpath(genpath('/Users/samst/Dropbox/UBC/Misc/'));
clear 

load S2.mat
load H1.mat

M=struct();
M.S2=S2; M.H1=H1;
clearvars -except M
load S2_17.mat
M.S2_17=S2;
clearvars -except M


M.H1.CTD133m.sigma_theta(32093:35236)=NaN;

%%
syr=year(M.S2.time); smnth=month(M.S2.time);
hyr=year(M.H1.time); hmnth=month(M.H1.time);
ddayH=day(datetime(datestr(M.H1.time)),'dayofyear');
ddayS=day(datetime(datestr(M.S2.time)),'dayofyear');

pr=ones(size(M.S2.CTD280m.salinity)).*(280*1.0051);
SA=gsw_SA_from_SP(M.S2.CTD280m.salinity,...
    pr,M.S2.longitude,M.S2.latitude);
PT=gsw_pt0_from_t(SA,M.S2.CTD280m.temperature,pr);
M.S2.CTD280m.AOU=aou(M.S2.CTD280m.salinity,PT,M.S2.CTD280m.oxygen);

pr=ones(size(M.H1.CTD133m.salinity)).*(133*1.0051);
SA=gsw_SA_from_SP(M.H1.CTD133m.salinity,...
    pr,M.H1.longitude,M.H1.latitude);
PT=gsw_pt0_from_t(SA,M.H1.CTD133m.temperature,pr);
M.H1.CTD133m.AOU=aou(M.H1.CTD133m.salinity,PT,M.H1.CTD133m.oxygen);

SclimWM=NaN(1,365);
HclimWM=NaN(1,365);
SclimWA=NaN(1,365);
HclimWA=NaN(1,365);
HrC=NaN(1,365);
SrC=NaN(1,365);
SclimWE=NaN(1,365);
HclimWE=NaN(1,365);
HrCE=NaN(1,365);
SrCE=NaN(1,365);
SclimAE=NaN(1,365);
HclimAE=NaN(1,365);

for i=1:365
    SclimWM(i)=mean(M.S2.CTD280m.oxygen(ddayS==i),'omitnan');
    HclimWM(i)=mean(M.H1.CTD133m.oxygen(ddayH==i),'omitnan');

    SclimWE(i)=std(M.S2.CTD280m.oxygen(ddayS==i),0,'omitnan');
    HclimWE(i)=std(M.H1.CTD133m.oxygen(ddayH==i),0,'omitnan');

    HrC(i)=mean(M.H1.CTD133m.sigma_theta(ddayH==i),'omitnan');
    SrC(i)=mean(M.S2.CTD280m.sigma_theta(ddayH==i),'omitnan');

    HrCE(i)=std(M.H1.CTD133m.sigma_theta(ddayH==i),0,'omitnan');
    SrCE(i)=std(M.S2.CTD280m.sigma_theta(ddayH==i),0,'omitnan');

    SclimWA(i)=mean(M.S2.CTD280m.AOU(ddayS==i),'omitnan');
    HclimWA(i)=mean(M.H1.CTD133m.AOU(ddayH==i),'omitnan');
    
    SclimAE(i)=std(M.S2.CTD280m.AOU(ddayS==i),0,'omitnan');
    HclimAE(i)=std(M.H1.CTD133m.AOU(ddayH==i),0,'omitnan');
end

% Create dummy timeseries
Hdummy=[repmat(HclimWM,1,7);repmat(HclimWM-1*HclimWE,1,7);repmat(HclimWM+1*HclimWE,1,7);...
    repmat(HclimWM-2*HclimWE,1,7);repmat(HclimWM+2*HclimWE,1,7)];
Sdummy=[repmat(SclimWM,1,7);repmat(SclimWM-1*SclimWE,1,7);repmat(SclimWM+1*SclimWE,1,7);
    repmat(SclimWM-2*SclimWE,1,7);repmat(SclimWM+2*SclimWE,1,7)];

HSTdummy=[repmat(HrCE,1,7);repmat(HrC-1*HrCE,1,7);repmat(HrC+1*HrCE,1,7);...
    repmat(HrC-2*HrCE,1,7);repmat(HrC+2*HrCE,1,7)];
SSTdummy=[repmat(SrC,1,7);repmat(SrC-1*SrCE,1,7);repmat(SrC+1*SrCE,1,7);
    repmat(SrC-2*SrCE,1,7);repmat(SrC+2*SrCE,1,7)];

dummyT=datenum(2017,01,01):datenum(2024,01,01)-2;



%% Creat 2d histograms
tmp=[M.S2_17.CTD280m.oxygen M.S2.CTD280m.oxygen];
timeEdges = floor(min(M.S2_17.time)):7:ceil(max(M.S2.time));
oxygenEdges = min(M.S2.CTD280m.oxygen):1:max(M.S2.CTD280m.oxygen);
[M.S2.N, M.S2.timeBinEdges, M.S2.oxygenBinEdges] =...
    histcounts2([M.S2_17.time M.S2.time], tmp, timeEdges, oxygenEdges);
M.S2.N(M.S2.N==0)=NaN;

tmp=[M.S2_17.CTD280m.sigma_theta M.S2.CTD280m.sigma_theta];
sigma_thetaEdges = min(M.S2.CTD280m.sigma_theta):0.025:max(M.S2.CTD280m.sigma_theta);
[M.S2.Ns, M.S2.timeBinEdges, M.S2.sigma_thetaBinEdges] =...
    histcounts2([M.S2_17.time M.S2.time], tmp, timeEdges, sigma_thetaEdges);
M.S2.Ns(M.S2.Ns==0)=NaN;

timeEdges = floor(min(M.H1.time)):7:ceil(max(M.H1.time));
oxygenEdges = min(M.H1.CTD133m.oxygen):1:max(M.H1.CTD133m.oxygen);
[M.H1.N, M.H1.timeBinEdges, M.H1.oxygenBinEdges] =...
    histcounts2(M.H1.time, M.H1.CTD133m.oxygen, timeEdges, oxygenEdges);
M.H1.N(M.H1.N==0)=NaN;

sigma_thetaEdges = min(M.H1.CTD133m.sigma_theta):0.025:max(M.H1.CTD133m.sigma_theta);
[M.H1.Ns, M.H1.timeBinEdges, M.H1.sigma_thetaBinEdges] =...
    histcounts2(M.H1.time, M.H1.CTD133m.sigma_theta, timeEdges, sigma_thetaEdges);
M.H1.Ns(M.H1.Ns==0)=NaN;

%% Mooring TS

figure('units','centimeters','outerposition',[0 0 18 21],'color','w');
cm1=cmocean('algae');cm2=cmocean('turbid'); 
cm1(1:3,:)=[0.95 0.95 0.95;0.95 0.95 0.95;0.95 0.95 0.95];
cm2(1:3,:)=[0.95 0.95 0.95;0.95 0.95 0.95;0.95 0.95 0.95];

%%%%%%% S2
ax1=axes('position',[0.1 0.75+0.035 0.5 0.15]); hold on
pcolor(M.S2.timeBinEdges(1:end-1), M.S2.oxygenBinEdges(1:end-1), M.S2.N');
shading flat; cm=cmocean('speed');colormap(gca,cm1);
f=fill([dummyT fliplr(dummyT)]',[Sdummy(4,:) fliplr(Sdummy(5,:))]',...
    [0.9 0.9 0.9],'edgecolor','none');
uistack(f,'bottom');

tmp=[M.S2_17.CTD280m.oxygen M.S2.CTD280m.oxygen];
plot([M.S2_17.time M.S2.time],smoothdata(tmp,'movmean',24*7,'omitnan'),...
        'color',rgb_x('black'),'LineWidth',0.5);
line([dummyT(1) dummyT(end)],[61 61],'linestyle','--','color','k');

tmp=smoothdata(M.S2.CTD280m.oxygen,'movmean',24*7,'omitnan');
msk2=find(tmp<interp1(dummyT,Sdummy(4,:),M.S2.time) & tmp<61);
scatter(M.S2.time(msk2),tmp(msk2),10,'filled', 'markeredgecolor', [0.2 0 0],...
    'MarkerFaceColor', [0.1 0 0],'markerfacealpha',0.3);  

grid on; axis tight; xticklabels([]);
ylabel('O$_2$ ($\mu$mol kg$^{-1}$)','Interpreter','latex');
title ('Scott2 (\it{z}\rm{=280 m)}','fontweight','normal'); yticks(50:25:175);
box on
axdate(7); xticklabels([]);xlim([datenum(2017,01,01) datenum(2023,10,01)]);
yl=ylim;
zoomD=[datenum(2023,05,01) datenum(2023,09,01)];
rectangle('position',[zoomD(1) yl(1) diff(zoomD) diff(yl)],'edgecolor',...
    rgb_x('red'));
text(0.01,0.95,'a)','units','normalized');


%%%%%% S2 O ZOOM
ax1m=axes('position',[0.615 0.75+0.035 0.13 0.15]); hold on
pcolor(M.S2.timeBinEdges(1:end-1), M.S2.oxygenBinEdges(1:end-1), M.S2.N');
shading flat; cm=cmocean('speed');colormap(gca,cm1);
f=fill([dummyT fliplr(dummyT)]',[Sdummy(4,:) fliplr(Sdummy(5,:))]',...
    [0.9 0.9 0.9],'edgecolor','none');
uistack(f,'bottom');
tmp=[M.S2_17.CTD280m.oxygen M.S2.CTD280m.oxygen];
plot([M.S2_17.time M.S2.time],smoothdata(tmp,'movmean',24*7,'omitnan'),...
        'color',rgb_x('black'),'LineWidth',0.5);
line([dummyT(1) dummyT(end)],[61 61],'linestyle','--','color','k');
grid on; 
axis tight;
xlim(zoomD);ylim(yl);yticks(50:25:175);
box on;  set(gca,'xcolor','r','ycolor','r')
axdate(3); yticklabels([]);xticklabels([]);
set(gca, 'GridLineStyle', '-', ...
         'GridColor', [0.15, 0.15, 0.15],  ...
         'GridAlpha',0.15);      
cbar = colorbar(ax1m,'location','east outside'); 
% cbar.Position = [0.1 0.75 0.1 0.05];
title(cbar,'\it{n}','fontweight','normal');


%%%%%% S2 density 
axb=axes('position',[0.1 0.60+0.035 0.5 0.14]); hold on

pcolor(M.S2.timeBinEdges(1:end-1), M.S2.sigma_thetaBinEdges(1:end-1), M.S2.Ns');
shading flat; cm=cmocean('speed');colormap(gca,cm);
f=fill([dummyT fliplr(dummyT)]',[SSTdummy(4,:) fliplr(SSTdummy(5,:))]',...
    [0.9 0.9 0.9],'edgecolor','none');
uistack(f,'bottom');
tmp=[M.S2_17.CTD280m.sigma_theta M.S2.CTD280m.sigma_theta];
plot([M.S2_17.time M.S2.time],smoothdata(tmp,'movmean',24*7,'omitnan'),...
        'color',rgb_x('black'),'LineWidth',0.5);

grid on; axis tight; xticklabels([]);
ylabel('$\sigma_\theta$ (kg m$^{-3}$)','Interpreter','latex');
box on
axdate(7);xlim([datenum(2017,01,01) datenum(2023,10,01)]);
ylim([26.3 26.9]); yl=ylim; yticks(26.4:0.2:26.8);
zoomD=[datenum(2023,05,01) datenum(2023,09,01)];
rectangle('position',[zoomD(1) yl(1) diff(zoomD) diff(yl)],'edgecolor',...
    rgb_x('red'));
text(0.01,0.95,'b)','units','normalized');
xlabel('Year');

% S2 DENSITY ZOOM
ax1m=axes('position',[0.615 0.60+0.035 0.13 0.14]); hold on
pcolor(M.S2.timeBinEdges(1:end-1), M.S2.sigma_thetaBinEdges(1:end-1), M.S2.Ns');
shading flat; cm=cmocean('speed');colormap(gca,cm);
f=fill([dummyT fliplr(dummyT)]',[SSTdummy(4,:) fliplr(SSTdummy(5,:))]',...
    [0.9 0.9 0.9],'edgecolor','none');
uistack(f,'bottom');
tmp=[M.S2_17.CTD280m.sigma_theta M.S2.CTD280m.sigma_theta];
plot([M.S2_17.time M.S2.time],smoothdata(tmp,'movmean',24*7,'omitnan'),...
        'color',rgb_x('black'),'LineWidth',0.5);
grid on; ylim(yl);
axis tight;
xlim(zoomD);
box on;  set(gca,'xcolor','r','ycolor','r')
axdate(3); yticklabels([]);
set(gca, 'GridLineStyle', '-', ...
         'GridColor', [0.15, 0.15, 0.15],  ...
         'GridAlpha', 0.15);
ylim(yl); yticks(26.4:0.2:26.8);

cbar = colorbar(ax1m,'location','east outside'); 
title(cbar,'\it{n}','fontweight','normal');
% cbar.Position = [0.1 0.75 0.1 0.05];

% Set the color of the x-axis tick labels
ax1m.XColor = 'k';  % Tick labels color
%%
% Set the color of the x-axis ruler (the line and tick marks)
ax1m.XRuler.Axle.ColorData = uint8([255; 0; 0; 255]);  % Ruler (line and tick marks) color

%%%%%%%%%%%%%% H1
axc=axes('position',[0.1 0.4 0.5 0.15]); hold on
pcolor(M.H1.timeBinEdges(1:end-1), M.H1.oxygenBinEdges(1:end-1), M.H1.N');
shading flat; cm=cmocean('speed');colormap(gca,cm1);
f=fill([dummyT fliplr(dummyT)]',[Hdummy(4,:) fliplr(Hdummy(5,:))]',...
    [0.9 0.9 0.9],'edgecolor','none');
ylim([25 225]);yl=ylim;yticks(50:50:200);

tmp=smoothdata(M.H1.CTD133m.oxygen,'movmean',24*7,'omitnan');
msk2=find(tmp<interp1(dummyT,Hdummy(4,:),M.H1.time) & tmp<61);
l=plot([M.H1.time(msk2); M.H1.time(msk2)],...
    [min(yl).*ones(size(M.H1.time(msk2)));max(yl).*ones(size(M.H1.time(msk2)))],...
    'color',[1 0 0 0.005],'linewidth',2);
uistack(l,'bottom');
uistack(f,'bottom');

plot(M.H1.time,smoothdata(M.H1.CTD133m.oxygen,'movmean',24*7,'omitnan'),...
        'color',rgb_x('black'),'LineWidth',0.5);

line([dummyT(1) dummyT(end)],[61 61],'linestyle','--','color','k');

tmp=smoothdata(M.H1.CTD133m.oxygen,'movmean',24*7,'omitnan');
msk2=find(tmp<interp1(dummyT,Hdummy(4,:),M.H1.time) & tmp<61);
% s=scatter(M.H1.time(msk2),tmp(msk2),10,'s', 'markeredgecolor', [1 0 0],...
%     'MarkerFaceColor', [1 0 0],'markerfacealpha',0.3);  
% uistack(s,'bottom');

grid on; axis tight;
title ('Hak1 (\it{z}\rm{=133 m)}','fontweight','normal'); ylim([25 225]);
ylabel('O$_2$ ($\mu$mol kg$^{-1}$)','Interpreter','latex');

Hso=M.H1.CTD133m.oxygen(msk2);
Hst=M.H1.CTD133m.sigma_theta(msk2);
Ha=M.H1.CTD133m.AOU(msk2);
Ht=M.H1.CTD133m.temperature(msk2);
Hs=M.H1.CTD133m.salinity(msk2);


box on;
axdate(6); yl=ylim;xlim([datenum(2017,01,01) datenum(2023,10,01)]);
rectangle('position',[zoomD(1) yl(1) diff(zoomD) diff(yl)],'edgecolor',...
    rgb_x('red'));
text(0.01,0.95,'c)','units','normalized');
axdate(7); xticklabels([]);xlim([datenum(2017,01,01) datenum(2023,10,01)]);


%%%%% ZOOM 
ax1m=axes('position',[0.615 0.4 0.13 0.15]); hold on
pcolor(M.H1.timeBinEdges(1:end-1), M.H1.oxygenBinEdges(1:end-1), M.H1.N');
shading flat; cm=cmocean('speed');colormap(gca,cm1);
f=fill([dummyT fliplr(dummyT)]',[Hdummy(4,:) fliplr(Hdummy(5,:))]',...
    [0.9 0.9 0.9],'edgecolor','none');
l=plot([M.H1.time(msk2); M.H1.time(msk2)],...
    [min(yl).*ones(size(M.H1.time(msk2)));max(yl).*ones(size(M.H1.time(msk2)))],...
    'color',[1 0 0 0.005],'linewidth',2);
uistack(l,'bottom');
uistack(f,'bottom');


plot(M.H1.time,smoothdata(M.H1.CTD133m.oxygen,'movmean',24*7,'omitnan'),...
        'color',rgb_x('black'),'LineWidth',0.5);

line([dummyT(1) dummyT(end)],[61 61],'linestyle','--','color','k');


% plot(M.H1.time(msk2),tmp(msk2),'r.'); 
grid on; 
axis tight;
xlim(zoomD);
box on;axdate(3); yticklabels([]);xticklabels([]);
set(gca,'xcolor','r','ycolor','r');
set(gca,'xticklabelrotation',0);
ylim(yl);yticks(50:50:200);

% Set the color of the x-axis tick labels
ax1m.XColor = 'k';  % Tick labels color
cbar = colorbar(ax1m,'location','east outside'); 
title(cbar,'\it{n}','fontweight','normal');

%%
% Set the color of the x-axis ruler (the line and tick marks)
ax1m.XRuler.Axle.ColorData = uint8([255; 0; 0; 255]);  % Ruler (line and tick marks) color
set(gca, 'GridLineStyle', '-', ...
         'GridColor', [0.15, 0.15, 0.15],  ...
         'GridAlpha', 0.15);

%%%%%% H1 density 
axd=axes('position',[0.1 0.4-0.15 0.5 0.14]); hold on

pcolor(M.H1.timeBinEdges(1:end-1), M.H1.sigma_thetaBinEdges(1:end-1), M.H1.Ns');
shading flat; cm=cmocean('speed');colormap(gca,cm);
f=fill([dummyT fliplr(dummyT)]',[HSTdummy(4,:) fliplr(HSTdummy(5,:))]',...
    [0.9 0.9 0.9],'edgecolor','none');
uistack(f,'bottom');
tmp=smoothdata(M.H1.CTD133m.sigma_theta,'movmean',24*7);
tmp(32093:35236)=NaN;
plot(M.H1.time,tmp,'-','color',rgb_x('black'),'LineWidth',0.5);

grid on; axis tight; xticklabels([]);
ylabel('$\sigma_\theta$ (kg m$^{-3}$)','Interpreter','latex');
box on
axdate(7);xlim([datenum(2017,01,01) datenum(2023,10,01)]);
zoomD=[datenum(2023,05,01) datenum(2023,09,01)];
text(0.01,0.95,'d)','units','normalized');
yl=[24.8 27];ylim(yl);yticks(25:27);
rectangle('position',[zoomD(1) yl(1) diff(zoomD) diff(yl)],'edgecolor',...
    rgb_x('red'));
xlabel('Year');

% H1 DENSITY ZOOM
ax1m=axes('position',[0.615 0.4-0.15 0.13 0.14]); hold on
pcolor(M.H1.timeBinEdges(1:end-1), M.H1.sigma_thetaBinEdges(1:end-1), M.H1.Ns');
shading flat; cm=cmocean('speed');colormap(gca,cm);
f=fill([dummyT fliplr(dummyT)]',[HSTdummy(4,:) fliplr(HSTdummy(5,:))]',...
    [0.9 0.9 0.9],'edgecolor','none');
uistack(f,'bottom');
plot(M.H1.time,smoothdata(M.H1.CTD133m.sigma_theta,'movmean',24*7,'omitnan'),...
        'color',rgb_x('black'),'LineWidth',0.5);
grid on; 
axis tight;
xlim(zoomD);
box on;  set(gca,'xcolor','r','ycolor','r')
axdate(3); yticklabels([]);
set(gca, 'GridLineStyle', '-', ...
         'GridColor', [0.15, 0.15, 0.15],  ...
         'GridAlpha', 0.15);
ylim(yl);yticks(25:27);

% Set the color of the x-axis tick labels
ax1m.XColor = 'k';  % Tick labels color
cbar = colorbar(ax1m,'location','east outside'); 
title(cbar,'\it{n}','fontweight','normal');
%%
% Set the color of the x-axis ruler (the line and tick marks)
ax1m.XRuler.Axle.ColorData = uint8([255; 0; 0; 255]);  % Ruler (line and tick marks) color

set(findall(gcf, '-property', 'fontsize'), 'fontsize', 9);

% Mooring time-series hists
mnth=month(M.S2.time);
yr=year(M.S2.time);
yrs=2017:2023;
lnes=[rgb_x('steel blue');rgb_x('rose')];

pS=[];mS=[];
for i=1:length(yrs)
    
    axes('position', [0.7+0.04 0.87-((i-1)*0.1) 0.095 0.08],'TickDir','out','fontsize',8);
    hold on


    msk=yr==yrs(i) & mnth>=5 & mnth<=10;
    msk2=M.S2.CTD280m.oxygen(msk)<61;
    pS(i)=sum(msk2)/sum(msk);
    mS(i)=median(M.S2.CTD280m.oxygen(msk),'omitnan');
    tmp=M.S2.CTD280m.oxygen(msk);
    
    D{i,1}=M.S2.CTD280m.oxygen(msk);
    lqS(i)=quantile(tmp,0.25);
    uqS(i)=quantile(tmp,0.75);

    line([mS(i) mS(i)],[0 2000],'color',lnes(1,:),'linewidth',1.75)
    line([quantile(tmp,0.25) quantile(tmp,0.25)],[0 2000],'linestyle',...
        '-','color',lnes(1,:),'linewidth',0.75)
    line([quantile(tmp,0.75) quantile(tmp,0.75)],[0 2000],'linestyle',...
        '-','color',lnes(1,:),'linewidth',0.75)


    % subplot(length(yrs),2,i*2-1)
    histogram(tmp,40:2.5:160,'facecolor',lnes(1,:),...
        'edgecolor','none');
    grid on;
    ylim([0 1220]);
    yticks([0 500 1000])
    yticklabels([]);
    set(gca,'xcolor','none','ycolor','none','fontsize',8);
    box off
    xlim([40 90]);


    if i<length(yrs)
        xticklabels([]);
    end

    if i==1
        title('Scott2','fontweight','normal')
        text(0.025,0.975,'e)','units','normalized','fontsize',8);
    end
    if i==length(yrs)
        % yticklabels({'0';'';'1000'})
        % text(70,1000,num2str(yrs(i)),'fontsize',6);
        % text(90,1000,{num2str(yrs(i));'(May-Jun.)'},'fontsize',5,...
        %     'horizontalalignment','right');
        set(gca,'xcolor','k');
        xlabel('O$_2$ ($\mu$mol kg$^{-1}$)','Interpreter','latex','fontsize',8);
    end
end

mnth=month(M.H1.time);
yr=year(M.H1.time);
yrs=2017:2023;
pH=[];mH=[];
for i=1:length(yrs)
    if i==1
        axes('position', [0.78+0.06 0.87-((i-1)*0.1) 0.095 0.08],'TickDir','out','fontsize',8);
        ylim([0 1220]);
        yticks([0 500 1000]);
        grid on
        yticklabels([]);
        box off
        set(gca,'xcolor','none','ycolor','none');
        xlim([40 90]);
        title('Hak1','fontsize',8,'fontweight','normal')
         text(120,1000,{num2str(yrs(i));'(May-Oct.)'},'fontsize',7,...
                'horizontalalignment','right');
         text(0.025,0.975,'f)','units','normalized','fontsize',8);
         text(80,400,'No data','fontsize',7,'HorizontalAlignment','center');
         xlim([40 120]);
    else

        % if yrs(i)==2020
        %     keyboard
        % end

        msk=yr==yrs(i) & mnth>=5 & mnth<=10;
        msk2=M.H1.CTD133m.oxygen(msk)<61;
        pH(i)=sum(msk2)/sum(msk);
        mH(i)=median(M.H1.CTD133m.oxygen(msk),'omitnan');
        
        axes('position', [0.78+0.06 0.87-((i-1)*0.1) 0.095 0.08],'TickDir','out','fontsize',8);
        hold on
        tmp=M.H1.CTD133m.oxygen(msk);
        lqH(i)=quantile(tmp,0.25);
        uqH(i)=quantile(tmp,0.75);

        line([mH(i) mH(i)],[0 700],'color',lnes(2,:),'linewidth',1.75)
        line([quantile(tmp,0.25) quantile(tmp,0.25)],[0 700],'linestyle',...
            '-','color',lnes(2,:),'linewidth',0.75)
        line([quantile(tmp,0.75) quantile(tmp,0.75)],[0 700],'linestyle',...
            '-','color',lnes(2,:),'linewidth',0.75)

        histogram(M.H1.CTD133m.oxygen(msk),40:2.5:160,'facecolor',lnes(2,:),...
            'edgecolor','none');
        grid on;
        ylim([0 1220]);
        yticks([0 500 1000]);
        yticklabels([]);
        box off
        set(gca,'xcolor','none','ycolor','none');
        xlim([40 120]);
        xticks([60 100])

        if i<length(yrs)
            xticklabels([]);
        end

        if i==length(yrs)
            set(gca,'xcolor','k','fontsize',8);
            text(120,1000,{'2023';'(May-Jul.)'},'fontsize',7,...
                'horizontalalignment','right')
        else
            text(120,1000,{num2str(yrs(i));'(May-Oct.)'},'fontsize',7,...
                'horizontalalignment','right');
        end
    end

end

set(findall(gcf,'-property','fontname'),'fontname','Latin Modern Roman');

%%
export_fig /Users/samuelstevens/Dropbox/Hakai/gliders/oxygen/figures/moorings_R1.pdf -dpdf -nofontswap
% exportgraphics(gcf, '/Users/samuelstevens/Dropbox/Hakai/gliders/oxygen/figures/mooringsHR.pdf', ...
%     'ContentType', 'vector', 'Resolution', 300);

%% PP plots
figure('units','centimeters','outerposition',[0 0 9 11],'color','w');

t = tiledlayout(2,3,'TileSpacing','Compact','Padding','Compact');

msk=M.S2.time>datenum(2023,05,01) & M.S2.time<datenum(2023,08,01) & ...
    ~isnan(M.S2.CTD280m.oxygen) & ~isnan(M.S2.CTD280m.sigma_theta);
Sso=M.S2.CTD280m.oxygen(msk); Sst=M.S2.CTD280m.sigma_theta(msk);
Sao=M.S2.CTD280m.aou(msk);
k = boundary(double(Sso'),double(Sst'));

[N, edgesX, edgesY] = histcounts2(M.S2.CTD280m.oxygen,...
    M.S2.CTD280m.sigma_theta,20:5:250,25:0.05:27.05);
AXh1=nexttile;
set(gca,'TickDir','out'); 
hold on;
X = edgesX(1:end-1) + diff(edgesX)/2;
Y = edgesY(1:end-1) + diff(edgesY)/2;
p=pcolor(X, Y, N');p.FaceAlpha=0.8;shading flat; colormap(gca,cm1);
xlabel('O$_2$ ($\mu$mol kg$^{-1}$)','Interpreter','latex');
ylabel('\sigma_\theta (kg m^{-3})');
axis tight;
caxis([0 1300]);
ylim([25 27]);
box on
axis ij
yticks([25:0.5:27])
xlim([20 210]);
xticks(50:50:200);
set(gca,'XTickLabelRotation',0);
text(0.05,0.95,'a)','units','normalized');
set(gca, 'Layer', 'top','GridLineStyle', ':');
grid on;  
plot(Sso(k), Sst(k), '-', 'LineWidth', 1,'color',[0.9 00 0]);  % Draw boundary

[N, edgesX, edgesY] = histcounts2(M.S2.CTD280m.AOU,...
    M.S2.CTD280m.sigma_theta,50:5:300,25:0.05:27.05);
k = boundary(double(Sao'),double(Sst'));

AXh2=nexttile;
set(gca,'TickDir','out'); 
hold on;
X = edgesX(1:end-1) + diff(edgesX)/2;
Y = edgesY(1:end-1) + diff(edgesY)/2;
p=pcolor(X, Y, N');p.FaceAlpha=0.8;shading flat; colormap(gca,cm2);
xlabel('AOU ($\mu$mol kg$^{-1}$)','Interpreter','latex');
axis tight;
caxis([0 1300]);
ylim([25 27]); xlim([50 300]);
yticklabels([]);
box on
axis ij
yticks([25:0.5:27]);
set(gca,'XTickLabelRotation',0);
text(0.05,0.95,'b)','units','normalized');
set(gca, 'Layer', 'top','GridLineStyle', ':');
grid on;  
plot(Sao(k), Sst(k), '-', 'LineWidth', 0.5,'color',[0.9 00 0]);  % Draw boundary

[N, edgesX, edgesY] = histcounts2(M.S2.CTD280m.temperature,...
    M.S2.CTD280m.salinity,5:0.1:9,31.6:0.005:34.25);
nexttile;
set(gca,'TickDir','out'); 
hold on;
X = edgesX(1:end-1) + diff(edgesX)/2;
Y = edgesY(1:end-1) + diff(edgesY)/2;
p=pcolor(X, Y, N');p.FaceAlpha=0.8;shading flat; colormap(gca,cm2);
xlabel('AOU ($\mu$mol kg$^{-1}$)','Interpreter','latex');
axis tight;
caxis([0 200]);
ylim([32.5 34.25]); xlim([5.2 9]);
% yticklabels([]);
box on
yticks([25:0.5:27]);
set(gca,'XTickLabelRotation',0);
text(0.05,0.95,'b)','units','normalized');
set(gca, 'Layer', 'top','GridLineStyle', ':');
grid on;  
% plot(Sz(k), Sst(k), '-', 'LineWidth', 0.5,'color',[0.9 00 0]);  % Draw boundary


%%%%%%

[N, edgesX, edgesY] = histcounts2(M.H1.CTD133m.oxygen,...
    M.H1.CTD133m.sigma_theta,20:5:250,25:0.05:27.05);
Hso(isnan(Hso))=[];Hst(isnan(Hst))=[];Ha(isnan(Ha))=[]; 
Ht(isnan(Ht))=[];Hs(isnan(Hs))=[];

% Compute the boundary of the points
k = boundary(double(Hso'),double(Hst'));

% Modify the plotting commands for the oxygen vs sigma theta histogram
nexttile;
set(gca, 'TickDir', 'out'); hold on;
X = edgesX(1:end-1) + diff(edgesX)/2;
Y = edgesY(1:end-1) + diff(edgesY)/2;
p = pcolor(X, Y, N'); p.FaceAlpha = 0.8; shading flat; colormap(gca, cm1);
xlabel('O$_2$ ($\mu$mol kg$^{-1}$)', 'Interpreter', 'latex');
ylabel('\sigma_\theta (kg m^{-3})');
ylim([25 27]); xlim([20 210]);
plot(Hso(k), Hst(k), '-', 'LineWidth', 0.5,'color',[0.1 0.1 0.1]);  % Draw boundary
caxis([0 1300]);
box on;
axis ij;
yticks(25:0.5:27);
xticks(50:50:200);
set(gca, 'XTickLabelRotation', 0);
text(0.05,0.95,'c)','units','normalized');
set(gca, 'Layer', 'top','GridLineStyle', ':');
grid on;  

[N, edgesX, edgesY] = histcounts2(M.H1.CTD133m.AOU,...
    M.H1.CTD133m.sigma_theta,50:5:300,25:0.05:27.05);
nexttile;

% Compute the boundary of the points
k = boundary(double(Ha'), double(Hst'));
% Modify the plotting commands for the AOU vs sigma theta histogram
set(gca, 'TickDir', 'out'); hold on;
X = edgesX(1:end-1) + diff(edgesX)/2;
Y = edgesY(1:end-1) + diff(edgesY)/2;
p = pcolor(X, Y, N'); p.FaceAlpha = 0.8; shading flat; colormap(gca, cm2);
xlabel('AOU ($\mu$mol kg$^{-1}$)', 'Interpreter', 'latex');
ylim([25 27]); xlim([50 300]);
plot(Ha(k), Hst(k), '-', 'LineWidth', 0.5,'color',[0.1 0.1 0.1]);
caxis([0 1300]);
box on;
axis ij;
yticks(25:0.5:27);
yticklabels([]);
set(gca, 'XTickLabelRotation', 0);
text(0.05,0.95,'f)','units','normalized');
set(gca, 'Layer', 'top','GridLineStyle', ':');
grid on;

msk=month(M.H1.time)==7;
[N, edgesX, edgesY] = histcounts2(M.H1.CTD133m.temperature(msk),...
    M.H1.CTD133m.salinity(msk),5:0.1:9.5,31.6:0.005:34.25);
k = boundary(double(Ht'), double(Hs'));
nexttile;
set(gca,'TickDir','out'); 
hold on;
X = edgesX(1:end-1) + diff(edgesX)/2;
Y = edgesY(1:end-1) + diff(edgesY)/2;
p=pcolor(X, Y, N');p.FaceAlpha=0.8;shading flat; colormap(gca,cm2);
xlabel('AOU ($\mu$mol kg$^{-1}$)','Interpreter','latex');
axis tight;
% caxis([0 200]);
ylim([32.5 34.25]); xlim([5.2 9]);
% yticklabels([]);
box on
% yticks([25:0.5:27]);
set(gca,'XTickLabelRotation',0);
text(0.05,0.95,'b)','units','normalized');
set(gca, 'Layer', 'top','GridLineStyle', ':');
grid on;  
plot(Ht(k), Hs(k), '-', 'LineWidth', 0.5,'color',[0.9 00 0]);  % Draw boundary


axes(AXh1)
c=colorbar('position',[0.425,0.8,0.01,0.1],'fontsize',6);
c.Ticks=[0 1000];
title(c,'N','fontsize',6);

axes(AXh2)
c=colorbar('position',[0.85,0.8,0.01,0.1],'fontsize',6);
c.Ticks=[0 1000];
title(c,'N','fontsize',6);


set(findall(gcf, '-property', 'fontsize'), 'fontsize', 8);
set(findall(gcf, '-property', 'fontname'), 'fontname', 'Latin Modern Roman');

%%
% exportgraphics(gcf, '/Users/samuelstevens/Dropbox/Hakai/gliders/oxygen/figures/mooringsPPV2.pdf', ...
%     'ContentType', 'vector', 'Resolution', 300);
%% Calculate mean and STD for Oxygen and Sigma-Theta

% Extract month and year from the time variables
mnth = month(M.H1.time);
yr = year(M.H1.time);
yrs = 2017:2023;

% Initialize arrays to store results
pH = [];        % Proportion of low oxygen events at H1
mH = [];        % mean oxygen at H1
mH_MAD = [];    % std of oxygen at H1
mst = [];       % mean sigma-theta at H1
mst_std = [];   % std of sigma-theta at H1

pS = [];        % Proportion of low oxygen events at S2
mS = [];        % mean oxygen at S2
mS_std = [];    % std of oxygen at S2
mSst = [];      % mean sigma-theta at S2
mSst_std = [];  % std of sigma-theta at S2

% Loop over each year to calculate means and stds
for i = 1:length(yrs)
    % Create mask for the specified year and months (May to September)
    msk = yr == yrs(i);% & mnth>=5 & mnth<=10;
    
    % H1 Mooring Calculations
    % Proportion of oxygen measurements below 61 at H1
    msk2 = M.H1.CTD133m.oxygen(msk) < 61;
    pH(i) = sum(msk2) / sum(msk);
    
    % Oxygen data at H1 for the current mask
    data_H1_oxygen = M.H1.CTD133m.oxygen(msk);
    % mean oxygen at H1
    mH(i) = mean(data_H1_oxygen, 'omitnan');
    % std of oxygen at H1
    mH_std(i) = std(data_H1_oxygen, 'omitnan');
    
    % Sigma-theta data at H1 for the current mask
    data_H1_sigma_theta = M.H1.CTD133m.sigma_theta(msk);
    % mean sigma-theta at H1
    mst(i) = mean(data_H1_sigma_theta, 'omitnan');
    % std of sigma-theta at H1
    mst_std(i) = std(data_H1_sigma_theta, 'omitnan');
    
    % S2 Mooring Calculations
    % Proportion of oxygen measurements below 61 at S2
    msk2 = M.S2.CTD280m.oxygen(msk) < 61;
    pS(i) = sum(msk2) / sum(msk);
    
    % Oxygen data at S2 for the current mask
    data_S2_oxygen = M.S2.CTD280m.oxygen(msk);
    % mean oxygen at S2
    mS(i) = mean(data_S2_oxygen, 'omitnan');
    % std of oxygen at S2
    mS_std(i) = std(data_S2_oxygen, 'omitnan');
    
    % Sigma-theta data at S2 for the current mask
    data_S2_sigma_theta = M.S2.CTD280m.sigma_theta(msk);
    % mean sigma-theta at S2
    mSst(i) = mean(data_S2_sigma_theta, 'omitnan');
    % std of sigma-theta at S2
    mSst_std(i) = std(data_S2_sigma_theta, 'omitnan');
end

% Special handling for the S2_17 data in 2017
mnth_S2_17 = month(M.S2_17.time);
msk_S2_17 = year(M.S2_17.time) == 2017 & mnth_S2_17 >= 5 & mnth_S2_17 < 10;

% Oxygen data for S2_17 in 2017
data_S2_17_oxygen = M.S2_17.CTD280m.oxygen(msk_S2_17);
% mean oxygen at S2_17 in 2017
mS(1) = mean(data_S2_17_oxygen, 'omitnan');
% std of oxygen at S2_17 in 2017
mS_std(1) = std(data_S2_17_oxygen, 'omitnan');

% Sigma-theta data for S2_17 in 2017
data_S2_17_sigma_theta = M.S2_17.CTD280m.sigma_theta(msk_S2_17);
% mean sigma-theta at S2_17 in 2017
mSst(1) = mean(data_S2_17_sigma_theta, 'omitnan');
% std of sigma-theta at S2_17 in 2017
mSst_std(1) = std(data_S2_17_sigma_theta, 'omitnan');

% Save the results to a .mat file
save('/Users/samuelstevens/Dropbox/Hakai/gliders/data/mooringDistsMean.mat', ...
     'mS', 'mH', 'yrs', 'mH_std', 'mS_std', 'mst', 'mst_std', 'mSst', 'mSst_std','SclimWM');

%% Calculate Median and Median Absolute Deviation (MAD) for Oxygen and Sigma-Theta

% Extract month and year from the time variables
mnth = month(M.H1.time);
yr = year(M.H1.time);
yrs = 2017:2023;

% Initialize arrays to store results
pH = [];        % Proportion of low oxygen events at H1
mH = [];        % Median oxygen at H1
mH_MAD = [];    % MAD of oxygen at H1
mst = [];       % Median sigma-theta at H1
mst_MAD = [];   % MAD of sigma-theta at H1

pS = [];        % Proportion of low oxygen events at S2
mS = [];        % Median oxygen at S2
mS_MAD = [];    % MAD of oxygen at S2
mSst = [];      % Median sigma-theta at S2
mSst_MAD = [];  % MAD of sigma-theta at S2

% Loop over each year to calculate medians and MADs
for i = 1:length(yrs)
    % Create mask for the specified year and months (May to September)
    msk = yr == yrs(i);% & mnth>=5 & mnth<=10;
    
    % H1 Mooring Calculations
    % Proportion of oxygen measurements below 61 at H1
    msk2 = M.H1.CTD133m.oxygen(msk) < 61;
    pH(i) = sum(msk2) / sum(msk);
    
    % Oxygen data at H1 for the current mask
    data_H1_oxygen = M.H1.CTD133m.oxygen(msk);
    % Median oxygen at H1
    mH(i) = median(data_H1_oxygen, 'omitnan');
    % MAD of oxygen at H1
    mH_MAD(i) = median(abs(data_H1_oxygen - mH(i)), 'omitnan');
    
    % Sigma-theta data at H1 for the current mask
    data_H1_sigma_theta = M.H1.CTD133m.sigma_theta(msk);
    % Median sigma-theta at H1
    mst(i) = median(data_H1_sigma_theta, 'omitnan');
    % MAD of sigma-theta at H1
    mst_MAD(i) = median(abs(data_H1_sigma_theta - mst(i)), 'omitnan');
    
    % S2 Mooring Calculations
    % Proportion of oxygen measurements below 61 at S2
    msk2 = M.S2.CTD280m.oxygen(msk) < 61;
    pS(i) = sum(msk2) / sum(msk);
    
    % Oxygen data at S2 for the current mask
    data_S2_oxygen = M.S2.CTD280m.oxygen(msk);
    % Median oxygen at S2
    mS(i) = median(data_S2_oxygen, 'omitnan');
    % MAD of oxygen at S2
    mS_MAD(i) = median(abs(data_S2_oxygen - mS(i)), 'omitnan');
    
    % Sigma-theta data at S2 for the current mask
    data_S2_sigma_theta = M.S2.CTD280m.sigma_theta(msk);
    % Median sigma-theta at S2
    mSst(i) = median(data_S2_sigma_theta, 'omitnan');
    % MAD of sigma-theta at S2
    mSst_MAD(i) = median(abs(data_S2_sigma_theta - mSst(i)), 'omitnan');
end

% Special handling for the S2_17 data in 2017
mnth_S2_17 = month(M.S2_17.time);
msk_S2_17 = year(M.S2_17.time) == 2017 & mnth_S2_17 >= 5 & mnth_S2_17 < 10;

% Oxygen data for S2_17 in 2017
data_S2_17_oxygen = M.S2_17.CTD280m.oxygen(msk_S2_17);
% Median oxygen at S2_17 in 2017
mS(1) = median(data_S2_17_oxygen, 'omitnan');
% MAD of oxygen at S2_17 in 2017
mS_MAD(1) = median(abs(data_S2_17_oxygen - mS(1)), 'omitnan');

% Sigma-theta data for S2_17 in 2017
data_S2_17_sigma_theta = M.S2_17.CTD280m.sigma_theta(msk_S2_17);
% Median sigma-theta at S2_17 in 2017
mSst(1) = median(data_S2_17_sigma_theta, 'omitnan');
% MAD of sigma-theta at S2_17 in 2017
mSst_MAD(1) = median(abs(data_S2_17_sigma_theta - mSst(1)), 'omitnan');

% Save the results to a .mat file
save('/Users/samuelstevens/Dropbox/Hakai/gliders/data/mooringDists.mat', ...
     'mS', 'mH', 'yrs', 'mH_MAD', 'mS_MAD', 'mst', 'mst_MAD', 'mSst', 'mSst_MAD');



%% Calculate annual IQR
mS_aIQR = [];
mH_aIQR = [];

mnth = month(M.H1.time);
yr = year(M.H1.time);
yrs = 2017:2023;

for i = 1:length(yrs)
    msk = yr == yrs(i);
    mS_aIQR(i)=iqr(M.S2.CTD280m.oxygen(msk));
    mH_aIQR(i)=iqr(M.H1.CTD133m.oxygen(msk)); 
end
mean(mS_aIQR)
std(mS_aIQR)
mean(mH_aIQR)
std(mH_aIQR)

%% Define Constants
hypoxia_threshold = 61; % O₂ threshold for hypoxia
clc

% Scott2 Mooring (280m) - Hypoxic Days & Mean O₂ with Uncertainty

% Pre-2022
pre_mask = M.S2.time < datenum(2022,01,01) & ~isnan(M.S2.CTD280m.oxygen);
pre_hypoxic_mask = pre_mask & M.S2.CTD280m.oxygen < hypoxia_threshold;
pre_hypoxic_days = sum(pre_hypoxic_mask) / 24;
total_days_S2_pre = sum(pre_mask) / 24;
mean_O2_pre = mean(M.S2.CTD280m.oxygen(pre_hypoxic_mask), 'omitnan');
std_O2_pre = std(M.S2.CTD280m.oxygen(pre_hypoxic_mask), 'omitnan');
uncertainty_O2_pre = 2 * std_O2_pre; % ±2 standard deviations

% Post-2022
post_mask = M.S2.time >= datenum(2022,01,01) & ~isnan(M.S2.CTD280m.oxygen);
post_hypoxic_mask = post_mask & M.S2.CTD280m.oxygen < hypoxia_threshold;
post_hypoxic_days = sum(post_hypoxic_mask) / 24;
total_days_S2_post = sum(post_mask) / 24;
mean_O2_post = mean(M.S2.CTD280m.oxygen(post_hypoxic_mask), 'omitnan');
std_O2_post = std(M.S2.CTD280m.oxygen(post_hypoxic_mask), 'omitnan');
uncertainty_O2_post = 2 * std_O2_post; % ±2 standard deviations

% Print Scott2 results
fprintf(['Scott2: %d/%d days (%.1f%%) hypoxic before 2022, mean O₂ = %.1f ± %.1f µmol/kg\n'], ...
        round(pre_hypoxic_days), round(total_days_S2_pre), (pre_hypoxic_days/total_days_S2_pre)*100, mean_O2_pre, uncertainty_O2_pre);
fprintf(['Scott2: %d/%d days (%.1f%%) hypoxic since 2022, mean O₂ = %.1f ± %.1f µmol/kg\n\n'], ...
        round(post_hypoxic_days), round(total_days_S2_post), (post_hypoxic_days/total_days_S2_post)*100, mean_O2_post, uncertainty_O2_post);

% Hak1 Mooring (133m) - Hypoxic Days & Mean O₂ with Uncertainty

% Pre-2022
pre_mask = M.H1.time < datenum(2022,01,01) & ~isnan(M.H1.CTD133m.oxygen);
pre_hypoxic_mask = pre_mask & M.H1.CTD133m.oxygen < hypoxia_threshold;
pre_hypoxic_days = sum(pre_hypoxic_mask) / 24;
total_days_H1_pre = sum(pre_mask) / 24;
mean_O2_pre = mean(M.H1.CTD133m.oxygen(pre_hypoxic_mask), 'omitnan');
std_O2_pre = std(M.H1.CTD133m.oxygen(pre_hypoxic_mask), 'omitnan');
uncertainty_O2_pre = 2 * std_O2_pre; % ±2 standard deviations

% Post-2022
post_mask = M.H1.time >= datenum(2022,01,01) & ~isnan(M.H1.CTD133m.oxygen);
post_hypoxic_mask = post_mask & M.H1.CTD133m.oxygen < hypoxia_threshold;
post_hypoxic_days = sum(post_hypoxic_mask) / 24;
total_days_H1_post = sum(post_mask) / 24;
mean_O2_post = mean(M.H1.CTD133m.oxygen(post_hypoxic_mask), 'omitnan');
std_O2_post = std(M.H1.CTD133m.oxygen(post_hypoxic_mask), 'omitnan');
uncertainty_O2_post = 2 * std_O2_post; % ±2 standard deviations

% Print Hak1 results
fprintf(['Hak1: %d/%d days (%.1f%%) hypoxic before 2022, mean O₂ = %.1f ± %.1f µmol/kg\n'], ...
        round(pre_hypoxic_days), round(total_days_H1_pre), (pre_hypoxic_days/total_days_H1_pre)*100, mean_O2_pre, uncertainty_O2_pre);
fprintf(['Hak1: %d/%d days (%.1f%%) hypoxic since 2022, mean O₂ = %.1f ± %.1f µmol/kg\n'], ...
        round(post_hypoxic_days), round(total_days_H1_post), (post_hypoxic_days/total_days_H1_post)*100, mean_O2_post, uncertainty_O2_post);

%% Zoom on event

load HakaiPruthDockProvisional.mat
PD=HakaiPruthDockProvisional;
load SWSadcpH1.mat
%%
msk=M.H1.time>datenum(2023,06,01);

figure('color','w');
t = tiledlayout(8, 1, 'TileSpacing', 'Compact', 'Padding', 'Compact');
nexttile
rectangle('position',[datenum(2023,07,18) min(M.H1.CTD133m.temperature(msk))...
    14 range(M.H1.CTD133m.temperature(msk))],'edgecolor',...
    rgb_x('red'),'facecolor','r','facealpha',0.05 );
hold on
plot(M.H1.time(msk),M.H1.CTD133m.temperature(msk),'r');
grid on; axis tight;
axdate
ylabel('T (^oC)')

nexttile
rectangle('position',[datenum(2023,07,18) min(M.H1.CTD133m.salinity(msk))...
    14 range(M.H1.CTD133m.salinity(msk))],'edgecolor',...
    rgb_x('red'),'facecolor','r','facealpha',0.05 );
hold on
plot(M.H1.time(msk),M.H1.CTD133m.salinity(msk),'b');
grid on; axis tight;
axdate
ylabel('SA (g/kg)')

nexttile
rectangle('position',[datenum(2023,07,18) min(M.H1.CTD133m.sigma_theta(msk))...
    14 range(M.H1.CTD133m.sigma_theta(msk))],'edgecolor',...
    rgb_x('red'),'facecolor','r','facealpha',0.05 );
hold on
plot(M.H1.time(msk),M.H1.CTD133m.sigma_theta(msk),'m');
grid on; axis tight;
axdate
ylabel('\sigma_\theta (kg/m^3)')

nexttile
rectangle('position',[datenum(2023,07,18) min(M.H1.CTD133m.spice(msk))...
    14 range(M.H1.CTD133m.spice(msk))],'edgecolor',...
    rgb_x('red'),'facecolor','r','facealpha',0.05 );
hold on
plot(M.H1.time(msk),M.H1.CTD133m.spice(msk),'k');
grid on; axis tight;
axdate
ylabel('Spice')

nexttile
rectangle('position',[datenum(2023,07,18) min(M.H1.CTD133m.oxygen(msk))...
    14 range(M.H1.CTD133m.oxygen(msk))],'edgecolor',...
    rgb_x('red'),'facecolor','r','facealpha',0.05 );
hold on
plot(M.H1.time(msk),M.H1.CTD133m.oxygen(msk),'g');
grid on; axis tight;
axdate
ylabel('O2 (umol/kg)');

nexttile

PD.mtime=datenum(1970,01,01,00,00,00)+(PD.time/86400);
msk=PD.mtime>datenum(2023,06,01) & PD.mtime<datenum(2023,08,01);
rectangle('position',[datenum(2023,07,18) min(PD.tideheight_avg(msk))...
    14 range(PD.tideheight_avg(msk))],'edgecolor',...
    rgb_x('red'),'facecolor','r','facealpha',0.05 );
hold on
plot(PD.mtime(msk),PD.tideheight_avg(msk),'color',[0.4 0.4 0.4]);
grid on; axis tight;
axdate
ylabel('Tidal Height (m)');

lat=ncread('dt_global_twosat_phy_l4_20230719_vDT2021.nc',...
    'latitude');
lon=ncread('dt_global_twosat_phy_l4_20230719_vDT2021.nc',...
    'longitude');

[~,Lidx]=min(abs(lon-M.H1.longitude));
[~,LLidx]=min(abs(lat-M.H1.latitude));

d=dir('/Users/samuelstevens/Dropbox/Hakai/gliders/data/SWOT/ssh/*.nc');
for i=1:length(d)
    fname=d(i).name;

    ssh(i)=squeeze(ncread(fname,...
        'sla',[Lidx LLidx 1],...
        [1 1 1],[1 1 1]));

    u(i)=squeeze(ncread(fname,...
        'ugos',[Lidx LLidx 1],...
        [1 1 1],[1 1 1]));

    v(i)=squeeze(ncread(fname,...
        'vgos',[Lidx LLidx 1],...
        [1 1 1],[1 1 1])); 

    satT(i)=ncread(fname,...
        'time')+datenum(1950,01,01);
end

nexttile
rectangle('position',[datenum(2023,07,18) min(ssh)...
    14 range(ssh)],'edgecolor',...
    rgb_x('red'),'facecolor','r','facealpha',0.05 );
hold on
plot(satT,ssh,'k');
grid on; axis tight; xlim([datenum(2023,06,01) datenum(2023,08,01)]);
axdate
ylabel('SLA @ HAK1 (m)')

nexttile
hold on

% Scale factor to adjust the size of the vectors for better visibility
sf = 5; % Adjust this value as needed

% Starting positions for the vectors
x = satT;
y = zeros(size(satT)); % Since we are plotting over time, y can be zeros

% Plot the event rectangle
rectangle('Position', [datenum(2023,07,18), -0.1, 14, 0.2], 'EdgeColor', rgb_x('red'), 'FaceColor', 'r', 'FaceAlpha', 0.05);


HT=datenum(adcp{1, 1}.Time);
msk=HT>datenum(2023,06,01) & HT<datenum(2023,08,01);
HT=HT(msk);

subU=smoothdata(adcp{1, 1}.velu(msk,1),'movmean',30,'omitnan');
subV=smoothdata(adcp{1, 1}.velv(msk,1),'movmean',30,'omitnan');

% Plot the vector stick plot using quiver
q(1)=quiver(HT(1:5:end), zeros(size(HT(1:5:end))), subU(1:5:end), subV(1:5:end), 0, 'r','maxheadsize',10);

% Plot the vector stick plot using quiver
q(2)=quiver(x, y, u, v, 0, 'k','maxheadsize',10,'linewidth',1.5);

% Formatting the plot
grid on; axis tight;
axdate;
ylabel('Currents @ HAK1 (m/s)');
legend(q,'42m @ HAK1','0m GS','location','best','fontsize',6);

% Optional: Adjust y-axis limits if necessary
% ylim([min(min(u), min(v))-0.1, max(max(u*sf), max(v*sf))+0.1]);

% 
% keyboard
% % 
% exportgraphics(gcf, '/Users/samuelstevens/Dropbox/Hakai/gliders/oxygen/figures/extremeHyp.jpg', ...
%     'Resolution', 300);